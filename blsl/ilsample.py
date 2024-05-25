from io import BytesIO
from pathlib import Path
from sys import stdin, stdout, stderr
from tqdm import tqdm
from collections import Counter
import random
import gzip

from ._utils import fqparse

def main(argv=None):
    """Sample a fraction of read pairs from an interleaved fastq file
    """
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="/dev/fd/1", type=argparse.FileType("wb"),
            help="Output file")
    ap.add_argument("-f", "--frac", default=0.1, type=float,
            help='Fraction of reads to sample')
    ap.add_argument("-i", "--input", type=Path, default="/dev/fd/0",
            help="Input fastq file")
    ap.add_argument("-s", "--singleend", action="store_const", const=4, default=8, dest="nlines",
            help="Single-end reads (default: paired end, interleaved)")
    args=ap.parse_args(argv)
    
    if args.input.suffix == ".gz":
        stream = gzip.open(args.input, "rb")
    else:
        stream = open(args.input, "rb")

    total = 0
    sampled = 0
    for pair in tqdm(fqparse(stream, args.nlines), unit="read pairs"):
        total += 1
        if random.random() < args.frac:
            sampled += 1
            args.out.write(b"".join(pair))

    print(f"Done! Sampled {sampled} of {total} pairs", file=stderr)

if __name__ == "__main__":
    main()
