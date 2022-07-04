#!/usr/bin/env python3
from sys import argv
import argparse

from Bio import SeqIO

def falen_main(argv):
    """Tabulate the lengths of sequences in a FASTA file"""
    ap = argparse.ArgumentParser()
    ap.add_argument("fastas", help="FASTA file of contigs", nargs="+")
    args = ap.parse_args(argv)
    
    for infile in args.fastas:
        with open(infile) as fh:
            for seq in SeqIO.parse(fh, 'fasta'):
                print(infile, seq.name, len(seq.seq), sep="\t")
                
if __name__ == "__main__":
    falen_main()
