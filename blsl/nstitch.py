# Copyright (c) 2023 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import gzip
import sys
import argparse
from ._utils import fqparse, rc

def nstitch_main(argv=None):
    """Combine R1 + R2 into single sequences, with an N in the middle"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-s", "--seq", default="N",
            help="Sequence used to glue together R1 + R2 of each pair")
    ap.add_argument("-q", "--qual", default="!",
            help="Quality score for glue sequences")
    ap.add_argument("-I", "--input-interleaved", action="store_true",
            help="Only read R1 which contains interleaved sequences")
    ap.add_argument("r1", help="R1 fastq.gz")
    ap.add_argument("r2", help="R2 fastq.gz", nargs="?")
    args = ap.parse_args(argv)
    
    if args.input_interleaved:
        r1 = gzip.open(args.r1, "rt")
        seqs = zip(fqparse(r1), fqparse(r1))
    else:
        r1 = gzip.open(args.r1, "rt")
        r2 = gzip.open(args.r2, "rt")
        seqs = zip(fqparse(r1), fqparse(r2))
    for s1, s2 in seqs:
        n1, n2 = s1[0].lstrip("@"), s2[0].lstrip("@")
        i1 = n1.split(" ", 1)[0]
        i2 = n2.split(" ", 1)[0]
        if i1 != i2:
            print(f"ERROR: read IDs don't match: {i1} {i2}", file=sys.stderr)
            sys.exit(1)
        newseq = f"{s1[1]}{args.seq}{rc(s2[1])}"
        newqual = f"{s1[3]}{args.qual*len(args.seq)}{s2[3]}"
        assert len(s1[1]) == len(s1[3])
        assert len(s2[1]) == len(s2[3])
        assert len(newseq) == len(newqual)
        print(f"@{n1} {n2}")
        print(newseq)
        print("+")
        print(newqual)
    r1.close()
    r2.close()

if __name__ == "__main__":
    nstich_main()
