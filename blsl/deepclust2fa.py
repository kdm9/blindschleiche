# Copyright (c) 2023 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import csv
import sys
import argparse
from collections import defaultdict
from pathlib import Path
import gzip

from Bio.SeqIO import parse
from tqdm import tqdm


class tsv(csv.unix_dialect):
    delimiter = "\t"
    quotechar = '"'
    doublequote = False
    escapechar = "\\"
    skipinitialspace = True


def deepclust2fa_main(argv=None):
    """Split a .faa by the clusters diamond deepclust finds"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-z", "--min-cluster-size", type=int, default=1,
            help="Minimum cluster size")
    ap.add_argument("-b", "--buffer", action="store_true",
            help="Buffer sequences rather than writing (more ram, more speed)")
    ap.add_argument("-c", "--clusters", required=True,
            help="Diamond deepclust result")
    ap.add_argument("-f", "--faa", required=True,
            help="Fasta sequences (that you gave to diamond deepclust --db)")
    ap.add_argument("-o", "--outdir", required=True, type=Path,
            help="Output directory (files named by their centroids)")
    args = ap.parse_args(argv)

    seq2cent = {}
    cent2seq = defaultdict(list)
    with open(args.clusters) as fh:
        for rec in csv.reader(fh, dialect=tsv):
            if rec[0] == "centroid" and rec[1] == "member":
                continue
            c, s = rec[0:2]
            seq2cent[s] = c
            cent2seq[c].append(s)
    cent2seq_new = {}
    seq2cent_new = {}
    for cent, seqs in cent2seq.items():
        if len(seqs) >= args.min_cluster_size:
            cent2seq_new[cent] = seqs
            for seq in seqs:
                seq2cent_new[seq] = cent
    cent2seq = cent2seq_new
    seq2cent = seq2cent_new
    print(f"Read clusters, {len(cent2seq)} clusters with size >= {args.min_cluster_size}", file=sys.stderr)

    for cent in cent2seq:
        outf = args.outdir/ (c+".fa")
        with open(outf, "w") as fh:
            pass

    if args.faa.endswith(".gz"):
        fh = gzip.open(args.faa, "rt")
    else:
        fh = open(args.faa)

    buffers = defaultdict(list)
    for seq in tqdm(parse(fh, "fasta")):
        centroid = seq2cent.get(seq.id)
        if centroid is None:
            continue
        if args.buffer:
            buffers[centroid].append((seq.id, seq.seq))
        else:
            outf = args.outdir/ (centroid+".fa")
            with open(outf, "a") as ofh:
                print(f">{seq.id}", file=ofh)
                print(seq.seq, file=ofh)
    fh.close()

    for centroid, seqs in buffers.items():
        outf = args.outdir/ (centroid+".fa")
        with open(outf, "w") as ofh:
            for seq in seqs:
                print(f">{seq[0]}", file=ofh)
                print(seq[1], file=ofh)


if __name__ == "__main__":
    deepclust2fa_main()

