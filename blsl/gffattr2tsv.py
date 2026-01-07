#!/usr/bin/env python3
import argparse
from sys import stdout
from .gffparse import parseGFF3, write_line
from tqdm import tqdm


def main(argv=None):
    """Extract the attributes column into a TSV"""
    ap = argparse.ArgumentParser("gfftagsane")
    ap.add_argument("-o", "--output", type=argparse.FileType("wt"), default=stdout,
            help="Output gff file")
    ap.add_argument("-t", "--type", type=str, action="append",
            help="Only emit records of this type (give multiple times for multiple types)")
    ap.add_argument("-x", "--extended", action="store_true",
            help="Only emit records of this type")
    ap.add_argument("-f", "--fields", default=["ID",], action="append", type=str,
            help="Attribute tags to keep (case sensitive).")
    ap.add_argument("--lax", action="store_true",
            help="Parse messy GFFs")
    ap.add_argument("input", nargs="+")

    args = ap.parse_args(argv)

    fname = []
    if len(args.input) > 1:
        fname.append("filename")

    cols = []
    if args.extended:
        cols.extend(["seqid", "source", "type", "start", "end", "score", "strand", "phase"])
    print(*cols, *args.fields, *fname, sep="\t")

    for filename in args.input:
        for record in tqdm(parseGFF3(filename, return_as=dict, lax=args.lax)):
            if args.type and record["type"] not in args.type:
                continue
            v = []
            if args.extended:
                v.extend([record[c] for c in cols])
            v.extend([record["attributes"].get(k, "") for k in args.fields])
            if len(args.input) > 1:
                v.append(filename)
            print(*v, sep="\t")

if __name__ == "__main__":
    main()
