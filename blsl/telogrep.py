# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from Bio import SeqIO
import re
import argparse
from collections import Counter

def telogrep_main(argv=None):
    """Search contigs for known telomere repeats"""
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", help="FASTA file of contigs")
    ap.add_argument(
            "-w", "--window", type=int, default=500, metavar="BASES",
            help="Scan BASES (default 500) at the end of each contig for telomeres",
    )
    ap.add_argument(
            "-t", "--threshold", type=float, default=30, metavar="N",
            help="If >N% of the window is telomeric repeats, consider that a telomere",
    )
    ap.add_argument(
            "-o", "--omit-telomereless", action="store_true",
            help="Don't output for contigs with neither 5' or 3' telomeres"
    )
    ap.add_argument(
            "-T", "--tsv", action="store_true",
            help="Output as a TSV for easy reading into other tools."
    )
    args = ap.parse_args(argv)

    # These are taken from https://dx.doi.org/10.1093/gbe/evt019
    tel_forward = "C{2,4}T{1,2}A{1,3}"
    tel_reverse = "T{1,3}A{1,2}G{2,4}"

    window = args.window
    min_bp = window * args.threshold/100
    nrev = 0
    nfwd = 0
    nt2t = 0

    if args.tsv:
        print("contig", "contig_len", "n_bp_fwd", "teloseq_fwd", "n_bp_rev", "teloseq_rev", sep="\t")

    for record in SeqIO.parse(args.fasta, "fasta"):
        ctgid = str(record.description).split(" ")[0]
        sequence = str(record.seq).upper()

        fwds = re.findall(tel_forward, sequence[:window])
        fwdlens = [len(x) for x in fwds]
        fwdbp = sum(fwdlens)

        revs = re.findall(tel_reverse, sequence[-window:])
        revlens = [len(x) for x in revs]
        revbp = sum(revlens)

        n_t = 0
        fwdstr = ""
        mcfwd = ""
        if fwdbp > min_bp:
            mcfwd = Counter(fwds).most_common(1)[0][0]
            fwdstr = "{}bp of ({})".format(fwdbp, mcfwd)
            nfwd += 1
            n_t += 1
        else:
            fwdbp = 0

        revstr = ""
        mcrev = ""
        if revbp > min_bp:
            mcrev = Counter(revs).most_common(1)[0][0]
            revstr = "{}bp of ({})".format(revbp, mcrev)
            nrev += 1
            n_t += 1
        else:
            revbp = 0
        
        if n_t == 2:
            nt2t += 1

        if args.omit_telomereless and fwdstr == "" and revstr == "":
            continue
        
        if args.tsv:
            print(ctgid, len(sequence), fwdbp, mcfwd, revbp, mcrev, sep="\t")
        else:
            print('{:<20} ({:>9}bp) {:>20} ------- {:<20}'.format(ctgid, len(sequence), fwdstr, revstr))
    
    if not args.tsv:
        print()
        print('{:<20}               {:>20} ------- {:<20}  ({} T2T chroms)'.format(
                "TOTAL:", "{} 5' Telos".format(nfwd), "{} 3' Telos".format(nrev), nt2t))

if __name__ == "__main__":
    telogrep_main()
