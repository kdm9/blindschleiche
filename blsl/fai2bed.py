# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import argparse
import sys


def parsefai(fai):
    with open(fai) as fh:
        for l in fh:
            cname, clen, _, _, _ = l.split()
            clen = int(clen)
            yield cname, clen


def fai2bed_main(argv=None):
    """Create a BED file from a samtools fasta index (.fai)"""
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", "-w", type=int, default=None,
            help="Emit non-overlapping windows of SIZE bp (default: whole sequence)")
    ap.add_argument("fai", help="FASTA index file (.fai)")
    args = ap.parse_args(argv)

    for cname, clen in parsefai(args.fai):
        if args.window is None:
            print(cname, 0, clen, sep="\t")
        else:
            for start in range(0, clen, args.window):
                end = min(start + args.window, clen)
                print(cname, start, end, sep="\t")


if __name__ == "__main__":
    fai2bed_main()
