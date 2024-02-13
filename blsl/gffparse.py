#!/usr/bin/env python3
"""
A simple parser for the GFF3 format.

Test with transcripts.gff3 from
http://www.broadinstitute.org/annotation/gebo/help/gff3.html.

Format specification source:
http://www.sequenceontology.org/gff3.shtml

Version 1.1: Python3 ready
Version 2.0: @kdm9's version. Add code for comparison & normalisation
"""
import argparse
from collections import namedtuple, Counter
import gzip
import urllib.request, urllib.parse, urllib.error
from sys import stderr, stdout
from copy import deepcopy

from tqdm.auto import tqdm

__author__  = "Uli Köhler"
__license__ = "Apache License v2.0"
__version__ = "2.1"

#Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

ontological_order = [
    "chromosome", "gene", "pseudogene", "pseudogenic_region", "mRNA", "transcript", "pseudogeneic_transcript",
    "exon", "pseudogenic_exon", "five_prime_UTR", "three_prime_UTR", "CDS", "pseudogenic_CDS",
]

def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""#
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        if "=" in attribute:
            key, value = attribute.split("=")
        else:
            key = attribute
            value = "True"
        if not key and not value:
            continue
        ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
    return ret

def parseGFF3(filename, return_as=dict):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
    
    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    openFunc = gzip.open if str(filename).endswith(".gz") else open
    with openFunc(filename, "rt") as infile:
        for line in infile:
            #if line.startswith("###"):
            #    ### Yield last gene if we move that here
            #    pass
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else urllib.parse.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.parse.unquote(parts[1]),
                "type": None if parts[2] == "." else urllib.parse.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.parse.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.parse.unquote(parts[7]),
                "attributes": parseGFFAttributes(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield return_as(**normalizedInfo)
            

def gff_heirarchy(filename, progress=None, make_missing_genes=False):
    last = {"l1": None, "l2": None, "l3": None}
    level_canonicaliser = {
        "transcript": "mRNA",
    }
    levels = {
        "repeat_region": "l1",
        "pseudogene": "l1",
        "pseuodgenic_region": "l1",
        "transposable_element_gene": "l1",
        "gene": "l1",
        "tRNA": "l2",
        "tmRNA": "l2",
        "mRNA": "l2",
        "exon": "l3",
        "CDS": "l3",
        "five_prime_UTR": "l3",
        "three_prime_UTR": "l3",
    }
    ignore = {
        "source",
        "stop_codon",
    }
    records = {}
    l2l1 = {}
    i = 0
    warned = set()
    recordsrc =  parseGFF3(filename, return_as=dict)
    if progress is not None:
        if progress == True:
            desc=None
        else:
            desc=progress
        recordsrc = tqdm(recordsrc, desc=desc)
    for i, record in enumerate(recordsrc):
        typ = record["type"]
        typ = level_canonicaliser.get(typ, typ)
        record["type"] = typ
        lvl = None
        if "RNA" in typ.upper():
            lvl = "l2"
        lvl = levels.get(typ, lvl)
        if lvl is None:
            if typ not in ignore:
                print(f"WARNING: {typ} is not a nice feature, skipping", file=stderr)
                ignore.add(typ)
            continue
        try:
            id = record["attributes"]["ID"]
        except:
            if "fakeid" not in warned:
                print(f"WARNING: no id for {record}, faked to fakeid{i}")
                warned.add("fakeid")
            id = f"fakeid{i}"
        if lvl == "l1":
            record["children"] = {}
            records[id] = record
            i = 0
        elif lvl == "l2":
            try:
                parent = record["attributes"]["Parent"]
            except KeyError:
                if make_missing_genes:
                    parent = f"fakegenefor_{id}"
                    record["attributes"]["Parent"] = parent
                    records[parent] = deepcopy(record)
                    records[parent]["type"] = "gene"
                    records[parent]["attributes"] = {"ID": parent}
                    records[parent]["children"] = {}
                else:
                    parent = None
            l2l1[id] = parent
            record["children"] = {}
            try:
                records[parent]["children"][id] = record
            except KeyError:
                print(f"L2 entry {id} parent {parent} not in records? {record}")
        else:
            try:
                parent = record["attributes"]["Parent"]
                top = l2l1[parent]
                if id in records[top]["children"][parent]["children"]:
                    i += 1
                    id = f"{id}_{i}"
                    record["attributes"]["ID"] = id
                records[top]["children"][parent]["children"][id] = record
            except KeyError:
                print(f"L3 entry {id} parent not in records? {record}")
    return records


def attr2line(attrs):
    from urllib.parse import quote as Q
    return ";".join("=".join((Q(k), Q(v))) for k, v in attrs.items())


def reformat_names(gene, geneid=None, changenames=True):
    def prefix_name(entry, sid):
        name = entry["attributes"].get("Name")
        if name:
            name = f"{sid}_{name}"
        else:
            name = sid
        entry["attributes"]["Name"] = name
    if not geneid:
        geneid = gene["attributes"]["ID"]
    gene["attributes"]["ID"] = geneid
    subids = Counter()
    for child in gene.get("children", {}).values():
        ct = child["type"]
        subids[ct] += 1
        ci = subids[ct]
        subid = f"{geneid}_{ct}{ci:02d}"
        child["attributes"]["ID"] = subid
        child["attributes"]["Parent"] = geneid
        if changenames:
            prefix_name(child, subid)
        #print("E", child["attributes"]["Name"])
        gcids = Counter()
        for gchild in child.get("children", {}).values():
            gt = gchild["type"]
            gcids[gt] += 1
            gi = gcids[gt]
            gcid = f"{geneid}_{gt}{gi:02d}"
            #print("F", gt, gcid)
            gchild["attributes"]["ID"] = gcid
            gchild["attributes"]["Parent"] = subid
            #print("G", gcid, gchild["attributes"].get("Name"))
            if changenames:
                prefix_name(gchild, gcid)
            #print("H", gcid, gchild["attributes"]["Name"])
    

def write_gff_line(line, file=None):
    attr = ";".join(f"{k}={v}" for k, v in line["attributes"].items())
    cols = [line.get(f, ".") for f in gffInfoFields[:-1]] + [attr]
    cols = ["." if x is None else x for x in cols]
    print(*cols, file=file, sep="\t")

def write_line(entry, file):
    x = [entry[field] if entry[field] is not None else "." for i, field in enumerate(gffInfoFields)]
    x[-1] = attr2line(x[-1])
    print(*x, sep="\t", file=file)

def write_gene(gene, geneid=None, file=None, changenames=False):
    if geneid or changenames:
        reformat_names(gene, geneid, changenames)
    write_line(gene, file)
    for child in gene.get("children", {}).values():
        write_line(child, file)
        for gchild in child.get("children", {}).values():
            write_line(gchild, file)       
    print("###", file=file)
    

def get_all_features(annot, ftype):
    feat = set()
    if annot["type"] == ftype:
        feat.add((annot["start"], annot["end"], annot["strand"], annot["type"]))
    for child in annot.get("children", {}).values():
        feat.update(get_all_features(child, ftype))
    return feat


def features_equal(a, b, ftype):
    a = set(get_all_features(a, ftype))
    b = set(get_all_features(b, ftype))
    return a == b


def feature_distance(a, b, ftype):
    def overlap(a, b):
        return a[0] <= b[0] <= a[1] or a[0] <= b[1] <= a[1]
    a = list(get_all_features(a, ftype))
    b = list(get_all_features(b, ftype))
    ovl = list()
    noovl = list()
    b_with_ovl = set()
    for ai in range(len(a)):
        for bi in range(len(b)):
            if overlap(a[ai], b[bi]):
                ovl.append((a[ai], b[bi]))
                b_with_ovl.add(b[bi])
                break
        else:
            # a does not overlap with an b
            noovl.append(a[ai])
    for bi in range(len(b)):
        if b[bi] not in b_with_ovl:
            noovl.append(b[bi])
    dist = 0
    for A, B in ovl:
        dist += abs(A[1] - B[1]) + abs(A[0] - B[0])
    for x in noovl:
        dist += x[1] - x[0]
    return dist


def gffparse_main(argv=None):
    """Format a GFF sanely"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="/dev/stdout",
            help="Output file")
    ap.add_argument("-t", "--tabular", action="store_true",
            help="Convert gff including attributes to simple TSV")
    ap.add_argument("-c", "--prefix-chrom", action="store_true",
            help="Prefix gene names with chromsome name (e.g. useful for concatenating augustus results)")
    ap.add_argument("input", help="Input GFF")
    args = ap.parse_args(argv)

    with open(args.output, "w") as fh:
        if args.tabular:
            records =  list(tqdm(parseGFF3(args.input, return_as=dict), desc="Parse GFF"))
            fields = set()
            for line in records:
                fields.update(line["attributes"].keys())
            print(*gffInfoFields[:-1], *fields, sep="\t", file=fh)
            for line in tqdm(records, desc="Tabularise"):
                cols = [line[x] for x in gffInfoFields[:-1]] + [line["attributes"].get(x) for x in fields]
                cols = [x if x is not None else "" for x in cols]
                print(*cols, sep="\t", file=fh)
        else:
            gff = gff_heirarchy(args.input, progress="Parse   ")
            for _, gene in tqdm(gff.items(), desc="Process "):
                if args.prefix_chrom:
                    newgid = f"{gene['seqid']}_{gene['attributes']['ID']}"
                    reformat_names(gene, geneid=newgid, changenames=False)
                write_gene(gene, file=fh)

if __name__ == "__main__":
    gffparse_main()
