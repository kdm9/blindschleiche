#!/usr/bin/env python3
import gzip
import io
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from concurrent.futures import as_completed, ProcessPoolExecutor
from sys import stdout
import argparse
import multiprocessing

from tqdm import tqdm

def head(raw, size=1_000_000):
    return io.BytesIO(raw.read(size))

@dataclass(order=True)
class FQStat:
    path: Path
    file_size: int
    estimated_nreads: float
    mean_read_len: float
    mean_record_size: float
    n_reads_sampled: int
    bytes_per_record: float


def estimate_fq_stats(fq, head_bytes=1_000_000):
    fq = Path(fq)
    with open(fq, "rb") as fh:
        buf = head(fh, size=head_bytes)
        bytes_read = len(buf.getbuffer())
        zfh = gzip.open(buf)
        n = 0
        recsizes = 0
        readlens = 0
        try:
            for hdr, seq, qhdr, qual in zip(zfh, zfh, zfh, zfh):
                n += 1
                recsizes += len(hdr) + len(seq) + len(qhdr) + len(qual)
                readlens += len(seq) -1 
        except (EOFError, gzip.BadGzipFile):
            pass
        fsize = fq.stat().st_size
        if n < 1:
            return FQStat(path=fq, file_size=fsize, estimated_nreads=0, mean_read_len=0, mean_record_size=0, n_reads_sampled=n, bytes_per_record=0)
        estim_reads = fsize / bytes_read * n
        return FQStat(path=fq, file_size=fsize, estimated_nreads=round(estim_reads),
                      mean_read_len=readlens/n, mean_record_size=recsizes/n,
                      n_reads_sampled=n, bytes_per_record=bytes_read/n)


def parse_fofn(file):
    with open(file) as fh:
        res = set()
        for fn in fh:
            res.add(fn.rstrip())
        return res


def main(argv=None):
    """Estimate stats from a fastq file based on the first kilobytes of the file, keeping high accuracy"""
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", "-o", type=argparse.FileType("w"), default=stdout,
        help="Output table")
    ap.add_argument("--threads", "-j", type=int, default=multiprocessing.cpu_count(),
        help="Parallel CPUs")
    ap.add_argument("--head", "-b", type=int, default=20_000,
        help="Inspect the first N bytes (default: 20kb, R^2 is still ~1.0!. Increase to improve accuracy slightly, above 1e6 is pointless.)")
    ap.add_argument("--fofn", "-f", action="store_true",
        help="Treat args as files of file names (one per line")
    ap.add_argument("fastqs", nargs="+")

    args = ap.parse_args(argv)

    if args.fofn:
        res = set()
        for fofn in args.fastqs:
            res.update(parse_fofn(fofn))
        args.fastqs = list(sorted(res))

    results = []
    with ProcessPoolExecutor(args.threads) as exc:
        jobs = set()
        for fq in args.fastqs:
            jobs.add(exc.submit(estimate_fq_stats, fq, head_bytes=args.head))
        for res in tqdm(as_completed(jobs)):
            results.append(res.result())

    print("path", "file_size", "estimated_n_reads", "read_length", "record_size", "bytes_per_record", sep="\t", file=args.out)
    for res in sorted(results):
        print(res.path, res.file_size, res.estimated_nreads, res.mean_read_len, res.mean_record_size, res.bytes_per_record, sep="\t", file=args.out)
