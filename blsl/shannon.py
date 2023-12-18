#!/usr/bin/env python3
from Bio import AlignIO
from tqdm import tqdm

import argparse
from sys import stdin, stdout, stderr
from math import log
from statistics import mean
from collections import Counter
warned = set()

class FreqCounter(Counter):
    def freqs(self):
        N = sum(self.values())
        return {b:v/N for b, v in self.most_common()}

def shannon_entropy(alncol):
    alncol = [b for b in alncol if b not in ("-", "*", "X")]
    aafreq = FreqCounter(alncol)
    entropy = -0
    for aa, freq in aafreq.freqs().items():
        entropy += freq*(log(freq,2))
    return len(alncol), -entropy

def shanent_aln(alignment, alnformat="fasta"):
    aln = AlignIO.read(alignment, alnformat)
    if len(aln) < 15 and "15_seqs" not in warned:
        print(f"WARN: {alignment} has <15 seqs, treat SE values with a grain of salt", file=stderr)
        warned.add("15_seqs")
    se = []
    for i in range(aln.get_alignment_length()):
        se.append(shannon_entropy(aln[:,i]))
    return aln, se

def main(argv=None):
    """Calculate Shannon's entropy (in bits) at each column of one or more alignments"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fofn", action="store_true",
            help="The first positional arg is a file of file names")
    ap.add_argument("-o", "--out-tsv", default=stdout, type=argparse.FileType("wt"),
            help="Output TSV file")
    ap.add_argument('alignments', nargs="+",
            help='One or more MSAs, in fasta format')
    args = ap.parse_args(argv)

    if args.fofn:
        fofn = args.alignments[0]
        if fofn == "-":
            fofn = "/dev/stdin"
        with open(fofn) as fh:
            args.alignments = [line.rstrip() for line in fh]

    print("alnfile", "nseqs",  "column", "shannon_entropy", sep="\t", file=args.out_tsv)
    for alnfile in tqdm(args.alignments, unit=" Alignment"):
        aln, se = shanent_aln(alnfile)
        for i, (l, s) in enumerate(se):
            print(alnfile, l, i+1, s, sep="\t", file=args.out_tsv)
        print(alnfile, len(aln), "*", mean([s for l, s in se]), sep="\t", file=args.out_tsv)

if __name__ == '__main__':
    main()
