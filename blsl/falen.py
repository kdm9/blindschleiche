#!/usr/bin/env python3
from sys import argv
import argparse

from Bio import SeqIO

def falen_main(argv):
    """Tabulate the lengths of sequences in a FASTA file"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-t", "--total", action="store_true",
            help="Output only a total length per file")
    ap.add_argument("fastas", help="FASTA file of contigs", nargs="+")
    args = ap.parse_args(argv)
    
    for infile in args.fastas:
        with open(infile) as fh:
            file_total = 0
            for seq in SeqIO.parse(fh, 'fasta'):
                file_total += len(seq.seq)
                if not args.total:
                    print(infile, seq.name, len(seq.seq), sep="\t")
            if args.total:
                print(infile, file_total, sep="\t")
                
if __name__ == "__main__":
    falen_main()
