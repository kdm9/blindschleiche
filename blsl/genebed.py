#!/usr/bin/env python3
import argparse
from sys import stdout
from .gffparse import parseGFF3, write_line
from tqdm import tqdm

def main(argv=None):
    """Extract a BED file of genes from a GFF"""
    ap = argparse.ArgumentParser("blsl genebed")
    ap.add_argument("-o", "--output", type=argparse.FileType("wt"), default=stdout,
            help="Output gff file")
    #ap.add_argument("-f", "--fields", default="ID,Name,Parent,locus_tag",
    #        help="Attribute tags to keep (case sensitive, give multiple times like -f ID -f tag2 -f tag3).")
    ap.add_argument("input")

    args = ap.parse_args(argv)
    #args.fields = args.fields.split(",")

    def N(x):
        if x is None:
            return "."
        return str(x)

    for record in tqdm(parseGFF3(args.input, return_as=dict)):
        if record["type"] != "gene":
            continue
        print(record["seqid"], record["start"], record["end"], record["attributes"]["ID"], N(record["score"]), N(record["strand"]), sep="\t", file=args.output)

if __name__ == "__main__":
    main()
