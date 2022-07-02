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
    args = ap.parse_args(argv)

    tel_forward = "C{2,4}T{1,2}A{1,3}"
    tel_reverse = "T{1,3}A{1,2}G{2,4}"

    window = args.window
    min_bp = window * args.threshold/100
    nrev = 0
    nfwd = 0

    for record in SeqIO.parse(args.fasta, "fasta"):
        ctgid = str(record.description).split(" ")[0]
        sequence = str(record.seq).upper()

        fwds = re.findall(tel_forward, sequence[:window])
        fwdlens = [len(x) for x in fwds]
        fwdbp = sum(fwdlens)

        revs = re.findall(tel_reverse, sequence[-window:])
        revlens = [len(x) for x in revs]
        revbp = sum(revlens)

        fwdstr = ""
        if fwdbp > min_bp:
            mcfwd = Counter(fwds).most_common(1)[0][0]
            fwdstr = "{}bp of ({})".format(fwdbp, mcfwd)
            nfwd += 1

        revstr = ""
        if revbp > min_bp:
            mcrev = Counter(revs).most_common(1)[0][0]
            revstr = "{}bp of ({})".format(revbp, mcrev)
            nrev += 1
        
        if args.omit_telomereless and fwdstr == "" and revstr == "":
            continue
        
        print('{:<20}  {:>20} ------- {:<20}'.format(ctgid, fwdstr, revstr))
    
    print()
    print('{:<20}  {:>20} ------- {:<20}'.format("TOTAL:", "{} 5' Telos".format(nfwd),
                                                           "{} 3' Telos".format(nrev)))

if __name__ == "__main__":
    telogrep_main()
