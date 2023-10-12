import argparse

def parsefai(fai):
    with open(fai) as fh:
        for l in fh:
            cname, clen, _, _, _ = l.split()
            clen = int(clen)
            yield cname, clen


def make_regions(fas, window=1e6):
    window = int(window)
    if not fas.endswith(".fai"):
        fas = fas + ".fai"
    for cname, clen in parsefai(fas):
        for start in range(0, clen, window):
            wlen = min(clen - start, window)
            yield (cname, start, start+wlen, wlen)

def main(argv):
    """Make a bed/region file of genome windows"""
    ap = argparse.ArgumentParser()
    ap.add_argument("--windowsize", "-w", default=100000, type=int,
            help="Window size of each non-overlapping region")
    ap.add_argument("--bed", "-b", type=argparse.FileType("w"),
            help="Output bed file")
    ap.add_argument("--regions", "-r", type=argparse.FileType("w"), 
            help="Output region (chr:start-stop) file")
    ap.add_argument("fasta", help="Fasta reference (must have fai)")
    args = ap.parse_args(argv)

    for chrm, start, stop, leng in make_regions(args.fasta, window=args.windowsize):
        region = f"{chrm}:{start+1}-{stop}"
        if args.bed is not None:
            print(chrm, start, stop, region, file=args.bed, sep="\t")
        if args.regions is not None:
            print(region, file=args.regions, sep="\t")
