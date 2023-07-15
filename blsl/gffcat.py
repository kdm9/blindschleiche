#!/usr/bin/env python3
import argparse
from sys import stderr, stdout, stdin

KNOWN_HEADERS = {
        "##gff-version",
        "##feature-ontology",
}


def gffcat_main(argv=None):
    """Concatenate GFF3 files, resepcting header lines and FASTA sections"""
    ap = argparse.ArgumentParser("gffcat")
    ap.add_argument("-o", "--output", type=argparse.FileType("wt"), default=stdout,
            help="Output gff file")
    ap.add_argument("inputs", nargs="+")
    args = ap.parse_args(argv)

    fastalines = []
    bodylines = []
    headerlines = {}

    for file in args.inputs:
        with open(file) as fh:
            in_fasta = False
            for line in fh:
                line = line.rstrip()
                if not line:
                    continue
                ll = line.lower().split()[0]
                if ll in KNOWN_HEADERS:
                    if ll in headerlines:
                        if headerlines[ll] != line:
                            print("WARN: header line with different value in differnt files. Only using value from first file", file=stderr)
                    else:
                        headerlines[ll] = line
                elif ll.startswith("##fasta"):
                    in_fasta = True
                else:
                    if in_fasta:
                        fastalines.append(line)
                    else:
                        bodylines.append(line)

    for hdr in KNOWN_HEADERS:
        if hdr in headerlines:
            print(headerlines[hdr], file=args.output)

    for line in bodylines:
        print(line, file=args.output)

    if fastalines:
        print("##fasta", file=args.output)

    for line in fastalines:
        print(line, file=args.output)
    
    print(f"DONE! {len(headerlines)} headers, {len(bodylines)} body lines, {len(fastalines)} FASTA lines", file=stderr)


if __name__ == "__main__":
    gffcat_main()
