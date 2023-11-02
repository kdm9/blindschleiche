#!/usr/bin/env python3
import argparse
from sys import stdout
from .gffparse import parseGFF3, write_line
from tqdm import tqdm

def main(argv=None):
    """Sanitise a messy gff attribute column to just simple tags """
    ap = argparse.ArgumentParser("gfftagsane")
    ap.add_argument("-o", "--output", type=argparse.FileType("wt"), default=stdout,
            help="Output gff file")
    ap.add_argument("-f", "--fields", default="ID,Name,Parent,locus_tag",
            help="Attribute tags to keep (case sensitive, give multiple times like -f ID -f tag2 -f tag3).")
    ap.add_argument("input")

    args = ap.parse_args(argv)
    args.fields = args.fields.split(",")

    for record in tqdm(parseGFF3(args.input, return_as=dict)):
        attrs = {k: record["attributes"][k] for k in args.fields if k in record["attributes"]}
        record["attributes"] = attrs
        write_line(record, args.output)

if __name__ == "__main__":
    main()
