# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from Bio import SeqIO
from shlex import quote
import argparse


full_help = """
Repeats or other problematic sequences are often masked from a reference
genome.  Ideally, this is done by converting each base to lower-case ("soft
masking"), though is sometimes done by replacing bases with N ("hard masking").

This tool reads a genome file, and for each sequence, outputs a zero-based bed
file of the masked regions. By default this means lower-case regions, although
if -N/--hardmasked is given, stretches of N are considered masked too.
"""

def mask2bed_main(argv=None):
    """The inverse of bedtools maskfasta: softmasked fasta -> unmasked fasta + mask.bed"""
    ap = argparse.ArgumentParser(epilog=full_help)
    ap.add_argument("-N", "--hardmasked", action="store_true",
            help="count N characters as a masked base (i.e. FASTA is hard-masked)")
    ap.add_argument("fasta", help="FASTA file of contigs")
    args = ap.parse_args(argv)
    
    seqs = SeqIO.parse(args.fasta, "fasta")
    for seq in seqs:
        last_start = 0
        last_masked = None
        for i, base in enumerate(str(seq.seq)):
            base_masked = base.islower() or (args.hardmasked and base.lower() == "n")
            if last_masked is None:
                last_masked = base_masked
                continue
            if last_masked != base_masked:
                if last_masked:
                    print(seq.name, last_start, i, sep="\t")
                else:
                    last_start = i
                last_masked = base_masked
        if last_masked:
            print(seq.name, last_start, i+1, sep="\t")

if __name__ == "__main__":
    mask2bed_main()
