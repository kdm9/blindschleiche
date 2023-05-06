# Copyright (c) 2023 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import csv
import sys
import math
import argparse

class blast6(csv.unix_dialect):
    delimiter = "\t"
    quotechar = '"'
    doublequote = False
    escapechar = "\\"
    skipinitialspace = True

def output_best(recs, fields):
    nbegin = len(recs)
    mineval =None
    maxalnlen=None
    if nbegin < 1:
        return
    mineval = min(float(x["evalue"]) for x in recs)
    recs = [
        rec for rec in recs if float(rec["evalue"]) <= mineval
    ]
    maxalnlen = max(float(x["length"]) for x in recs)
    recs = [
        rec for rec in recs if float(rec["length"]) >= maxalnlen
    ]
    nend=len(recs)
    print(recs[0]["qseqid"], nbegin, nend, mineval, maxalnlen, file=sys.stderr)
    for rec in recs:
        print(*[rec[x] for x in fields], sep="\t")
    

def equalbestblast_main(argv=None):
    """Output only the best blast hits."""
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--outfmt", default="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            help="Blast outfmt headers (default is the default for blast -outfmt 6)")
    ap.add_argument("table", help="Table of blast hits")
    args = ap.parse_args(argv)
    

    fields = args.outfmt.split(" ")
    query = None
    hits = []
    with open(args.table) as fh:
        for rec in csv.DictReader(fh, dialect=blast6, fieldnames=fields):
            if query != rec["qseqid"]:
                if query is not None:
                    output_best(hits, fields)
                hits = []
                query = rec["qseqid"]
            hits.append(rec)
        if query is not None:
            output_best(hits, fields)

if __name__ == "__main__":
    equalbestblast_main()

