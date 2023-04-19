#!/usr/bin/env python3
from sys import stdin, stdout, stderr, exit, argv
import re
import argparse

def main(argv):
    """Make a ncbi-style acc2taxid.map file for a uniref fasta"""

    prog = "blsl uniref-acc2taxid"
    ap = argparse.ArgumentParser(prog=prog)
    ap.add_argument("--strip-uniref", "-s", action="store_true",
            help="Strip UniRefXX_ prefixes from IDs (for Diamond, see https://github.com/bbuchfink/diamond/issues/639)")
    args = ap.parse_args(argv)
    
    if stdout.isatty() and stdin.isatty():
        print(f"USAGE: zcat uniref.fasta.gz | {prog} [-s] >uniref.acc2taxid.tsv")
        exit(0)

    #if args.strip_uniref:
    #    print("Stripping UniRefXX_ from each sequence ID (for Diamond)", file=stderr)
    print("accession.version", "taxid", sep="\t", file=stdout)
    for line in stdin:
        if not line.startswith(">"):
            continue
        seqid, *line = line.rstrip().lstrip(">").split(" ")
        fields = dict(x.split("=", 1) for x in line if "=" in x)
        if args.strip_uniref:
            # needed for Diamond, see https://github.com/bbuchfink/diamond/issues/639
            seqid = re.sub(r"UniRef\d\d_", "", seqid)
        print(seqid, fields.get("TaxID", 0), sep="\t", file=stdout)

if __name__ == "__main__":
    main(argv)
