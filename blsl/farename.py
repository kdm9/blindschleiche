#!/usr/bin/env python3
import argparse

from Bio import SeqIO

def farename_main(argv=None):
    """Rename sequences in a fasta file sequentially"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--prefix", type=str, default="",
            help="Prefix for each sequence ID (e.g. SEQ_ -> SEQ_1, SEQ_2 etc)")
    ap.add_argument("fastas", help="FASTA file of contigs", nargs="+")
    args = ap.parse_args(argv)
    
    for infile in args.fastas:
        i = 0
        with open(infile) as fh:
            for seq in SeqIO.parse(fh, 'fasta'):
                i += 1
                print(f">{args.prefix}{i} {seq.id} {seq.description}")
                print(seq.seq)

                
if __name__ == "__main__":
    farename_main()
