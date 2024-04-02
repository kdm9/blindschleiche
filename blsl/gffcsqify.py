#!/usr/bin/env python3
from blsl.gffparse import parseGFF3, write_gene, gffInfoFields, write_gff_line
import argparse


def main(argv=None):
    """Format a reasonably compliant GFF for use with bcftools csq"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="/dev/stdout",
            help="Output GFF compatible with bcftools csq")
    ap.add_argument("input", help="Input GFF")
    args = ap.parse_args(argv)

    with open(args.output, "w") as fh:
        for entry in parseGFF3(args.input):
            if entry["type"] == "gene":
                id = entry["attributes"]["ID"]
                entry["attributes"]["ID"] = f"gene:{id}"
                entry["attributes"]["biotype"] =  "protein_coding"  # FIXME properly parse this
            elif entry["type"] in ("mRNA", ):
                id = entry["attributes"]["ID"]
                parent = entry["attributes"]["Parent"]
                entry["attributes"]["ID"] = f"transcript:{id}"
                entry["attributes"]["Parent"] = f"gene:{parent}"
                entry["attributes"]["biotype"] =  "protein_coding"  # FIXME properly parse this
            elif entry["type"] in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
                parent = entry["attributes"]["Parent"]
                entry["attributes"]["Parent"] = f"transcript:{parent}"
                #if entry["type"] == "CDS" and entry["attributes"]["original_annotator"] == "Liftoff" and entry["phase"] is None:
                #    entry["phase"] = 0
            else:
                continue
            write_gff_line(entry, file=fh)


if __name__ == "__main__":
    main()
