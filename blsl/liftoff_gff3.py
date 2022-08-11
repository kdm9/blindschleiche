#!/usr/bin/env python3
import argparse
from collections import OrderedDict, defaultdict
from sys import stdin, stdout, stderr

reserved_keys = ["ID", "Name", "Alias", "Parent", "Target", "Gap",
        "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular"]

def fixcase(key):
    if key not in reserved_keys:
        return key.lower()
    else:
        return key

def liftoff_gff3_main(argv=None):
    """Obtain an actually-useful GFF3 from Liftoff by fixing basic GFF3 format errors"""
    ap = argparse.ArgumentParser("blsl liftoff-gff3")
    ap.add_argument("-o", "--output", type=argparse.FileType("w"), default="-",
            help="Output mostly-kosher GFF3 file")
    ap.add_argument("input", type=argparse.FileType("r"), default='-',
            help="Input liftoff 'gff' file.")
    args=ap.parse_args(argv)
    
    emptyattrcount = 0
    ids = defaultdict(set)
    replacement_suffixes = dict()
    replacements = dict()
    for lineno, line in enumerate(args.input):
        if lineno == 0 and not line.startswith("##gff-version"):
            print("##gff-version 3", file=args.output)
        if line.startswith("#"):
            args.output.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            print("WARNING: wrong number of columns on line", lineno, f"\n{line}", file=stderr)
        type = fields[2]
        kvp = OrderedDict([(fixcase(x[0]),x[1]) for x in [y.split("=") for y in fields[-1].split(";")]])
        if "" in kvp:  # sometime you get ID=***;=; type bullshit in there, clear it out
            emptyattrcount += 1
            del kvp[""]
        if "Parent" in kvp:
            parents = kvp["Parent"].split(",")
            for i, parent in enumerate(parents):
                if parent in replacements:
                    parents[i] = replacements[parent]
            parents = [p for p in parents if not p.endswith("-Protein")]  # remove protein parents as the aren't defined anywhere
            kvp["Parent"] = ",".join(parents)
        if "ID" in kvp:
            if type.lower() == "gene" and kvp["ID"] in ids[type] :
                # Many data types are *supposed* to have the same ID for a
                # block of discontinuous entries. We handle genes like this,
                # and then just add this number again to the exon/mRNA/CDS.
                oldname = kvp["ID"]
                newname = oldname
                for i in range(2, 9999999999):
                    if newname not in ids[type]:
                        break
                    newname = f"{oldname}_{i}"
                    replacement_suffixes[newname] = i
                ids[type].add(newname)
                kvp["ID"] = newname
                kvp["kdm_liftoff_fix_oldid"] = oldname
                print("WARNING: duplicate", type, "ID:", oldname, f"at {fields[0]}:{fields[3]}-{fields[4]} now called", newname, file=stderr)
                replacements[oldname] = newname
            elif "Parent" in kvp and kvp["Parent"] in replacement_suffixes:
                oldname = kvp["ID"]
                suf = replacement_suffixes[kvp['Parent']]
                newname = f"{oldname}_{suf}"
                kvp["ID"] = newname
                kvp["kdm_liftoff_fix_oldid"] = oldname
                print("WARNING: duplicate", type, "ID:", oldname, f"at {fields[0]}:{fields[3]}-{fields[4]} now called", newname, file=stderr)
                replacements[oldname] = newname
                replacement_suffixes[newname] = suf
                ids[type].add(kvp["ID"])
            else:
                ids[type].add(kvp["ID"])
        kvglue = ";".join(f"{k}={v}" for k, v in kvp.items())
        fields[-1] = kvglue
        print(*fields, sep="\t", file=args.output)
    if emptyattrcount > 0:
        print("Fixed empty attribute (;=;)", emptyattrcount, "times", file=stderr)
