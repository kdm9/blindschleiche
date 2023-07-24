# Copyright (c) 2023 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import csv
import sys
import math
import argparse


def tabcat_main(argv=None):
    """Concatenate table (c/tsv) files, adding the filename as a column"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--delim", default=",",
            help="Column delimiter.")
    ap.add_argument("-o", "--output", default="-", type=argparse.FileType("wt"),
            help="Column delimiter.")
    ap.add_argument("-N", "--add-name", action="store_true",
            help="Add a column with the file name")
    ap.add_argument("-f", "--fofn", action="store_true",
            help="Each input is a file of filenames")
    ap.add_argument("xsvs", help="Input {C,T}SVs to concatenate", nargs="+")
    args = ap.parse_args(argv)

    class xsv(csv.unix_dialect):
        delimiter = args.delim
        quotechar = '"'
        doublequote = True
        skipinitialspace = True
        quoting=csv.QUOTE_MINIMAL

    if args.fofn:
        files = []
        for fofn in args.xsvs:
            with open(fofn) as fh:
                for file in fh:
                    files.append(file.rstrip())
    else:
        files = args.xsvs
    outcsv = None
    for file in files:
        with open(file) as fh:
            for rec in csv.DictReader(fh, dialect=xsv):
                rec["tabcat_filename"] = file
                if outcsv is None:
                    outcsv = csv.DictWriter(args.output, dialect=xsv, fieldnames=[k for k in rec])
                    outcsv.writeheader()
                outcsv.writerow(rec)

if __name__ == "__main__":
    tabcat_main()


