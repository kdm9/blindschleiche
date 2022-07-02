# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from Bio import SeqIO
from shlex import quote
import argparse

def calc_n50(fastafile):
    lengths = []
    seqs = SeqIO.parse(fastafile, "fasta")
    for seq in seqs:
        lengths.append(len(seq.seq))
    lengths.sort()
    asmlen = sum(lengths)
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= asmlen/2:
            return sum(lengths), len(lengths), l

def n50_main(argv=None):
    """Calculate N50 and total length of a set of contigs"""
    ap = argparse.ArgumentParser()
    ap.add_argument("assembly", nargs="+")
    args = ap.parse_args(argv)
    
    for file in args.assembly:
        nbases, ncontig, n50 = calc_n50(file)
        print(quote(file), ":", sep="")
        print("  ncontigs:", ncontig)
        print("  n50:", n50)
        print("  total_length:", nbases)


if __name__ == "__main__":
    n50_main()
