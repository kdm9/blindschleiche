from io import BytesIO
from pathlib import Path
from sys import stdin
from tqdm import tqdm
from collections import Counter

from ._utils import fqpair, rc

def hamming_cloud(seq, distance=1):
    """Generate DNA seqs whose Hamming distance from seq is <= distance

    >>> sorted(hamming_cloud('GAT', 0))
    ['GAT']
    >>> sorted(hamming_cloud('GAT', 1))
    ['AAT', 'CAT', 'GAA', 'GAC', 'GAG', 'GAN', 'GCT', 'GGT', 'GNT', 'GTT', 'NAT', 'TAT']
    >>> max(hammdist(x, "GATTACA") for x in hamming_cloud('GAT', 2))
    2
    """
    from itertools import combinations, product
    seqs = set([seq, ])
    if distance == 0:
        return seqs
    elif distance > 0:
        seqs.update(hamming_cloud(seq, distance=distance-1))
    ACGTN = "ACGTN"
    for idx in combinations(range(len(seq)), distance):
        for rep in product(range(len(ACGTN) - 1), repeat=distance):
            mutated = list(seq)
            for p, r in zip(idx, rep):
                if mutated[p] == ACGTN[r]:
                    # we have N-1 replacements (hence -1 above), so we only
                    # cacluated the product above for N-1. So when we find the
                    # "replacement" that keeps the base the same, it's actually
                    # the replacement by the final letter in the alphabet.
                    r = -1
                mutated[p] = ACGTN[r]
            seqs.add(''.join(mutated))
    return seqs


class ByteBucket(object):
    """
    This bullshit is needed to avoid having too many files open (linux
    typically limits this to 1024). This thing buffers lines and only opens one
    file at a time to write output.

    And yes it needs pigz installed, if you're reading this that's why it
    crashed.
    """
    threads = 32
    ziplevel = 6
    def biff2pigz(self, bytes, outfile):
        from subprocess import Popen, PIPE
        with Popen(["pigz", "-n", "-p", str(self.threads), f"-{self.ziplevel}"], stdin=PIPE, stdout=outfile, bufsize=0) as proc:
            proc.stdin.write(bytes)

    def __init__(self, outfile, size=2**27, sname=""):
        self.outfile = str(outfile)
        self.buf = BytesIO()
        self.maxsize = size
        self.first = True
        self.sname = sname

    def writelines(self, lines):
        self.write("".join(lines).encode('ascii'))

    def write(self, somebytes):
        if self.first:
            # test write
            with open(self.outfile, "wb") as fh:
                pass
            self.first = False
        if isinstance(somebytes, str):
            somebytes.encode('ascii')
        self.buf.write(somebytes)

        if len(self.buf.getbuffer()) > self.maxsize:
            self.flush()

    def flush(self):
        if len(self.buf.getbuffer()) > 0:
            assert(self.first)
            with open(self.outfile, "ab", buffering=0) as fh:
                if self.outfile.endswith(".gz"):
                    self.biff2pigz(self.buf.getbuffer(), fh)
                else:
                    fh.write(self.buf.getbuffer())
                self.buf = BytesIO()

    def __del__(self):
        self.flush()

def gethdrdat(hdr):
    spot, idx = hdr.rstrip("\n").split(" ")
    idxpair = idx.split(":")[3]
    return spot, idx, idxpair

def fqp2idx(fqp):
    """fastq pair to index pair ACGT+ACGT"""
    s1, i1, ip1 = gethdrdat(fqp[0])
    s2, i2, ip2 = gethdrdat(fqp[4])
    assert(s1 == s2)
    assert(ip1 == ip2)
    return ip1

def make_sample_map(tablefile, outdir, fileext=".fq", distance=1):
    from csv import DictReader
    from itertools import product
    smap = {}
    with open(tablefile) as tabfh:
        for sample in DictReader(tabfh, dialect="excel-tab"):
            sname = sample["sample"]
            opath = outdir / f"{sname}.{fileext}"
            bucket = ByteBucket(opath, sname=sname)
            for i7, i5 in product(hamming_cloud(sample["i7"], distance=distance),
                                  hamming_cloud(rc(sample["i5"]), distance=distance)):
                i7plusi5 = f"{i7}+{i5}"
                if i7plusi5 in smap:
                    print("ERROR: duplicate i7+i5. reduce -m or cry into your beer.")
                    sys.exit(1)
                smap[i7plusi5] = bucket
    return smap


def main():
    """Demultiplex modern illumina reads from read headers.

    Supports *ONLY* out-of-read indicies (in @header lines e.g. 1:N:0:ACGT+TGCA), and only reads with i7 + i5.
    """

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outdir", required=True, type=Path,
            help="Output directory (must exist)")
    ap.add_argument("-c", "--justcount", action="store_true",
            help="Compress outputs with gzip.")
    ap.add_argument("-z", "--zip", nargs="?", default=None, const=6, type=int,
            help="Compress outputs with gzip.")
    ap.add_argument("-j", "--threads", default=8,
            help="Number of compression threads")
    ap.add_argument("-m", "--mismatch", default=1, type=int,
            help="Number tolerable mismatches (hamming distance).")
    ap.add_argument("-k", "--keyfile",
            help="Mapping of i7/i5 -> sample as tsv (with cols i7, i5, sample)")
    args=ap.parse_args()

    ByteBucket.threads = args.threads
    if args.zip is not None:
        ByteBucket.ziplevel = args.zip

    fileext = "il.fastq.gz" if args.zip is not None else "il.fastq"

    samps = make_sample_map(args.keyfile, args.outdir, fileext=fileext, distance=args.mismatch)
    print("set up sample map with", len(samps), "mappings of barcodes to files")

    stats = Counter()
    file_undef = ByteBucket(args.outdir/f"undefined.{fileext}", sname="undefined")
    file_undef.sname = "undefined"
    try:
        for i, pair in enumerate(tqdm(fqpair(stdin))):
            idxpair = fqp2idx(pair)
            ofile = samps.get(idxpair, file_undef)
            if not args.justcount:
                ofile.writelines(pair)
            stats[ofile.sname] += 1
    finally:
        with open(args.outdir / "stats.tsv", "w") as fh:
            print("sample", "n_pairs", file=fh, sep="\t")
            for samp, n in stats.most_common():
                print(samp, n, file=fh, sep="\t")

if __name__ == "__main__":
    main()
