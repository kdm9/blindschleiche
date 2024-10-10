#!/usr/bin/env python3
import gzip
import io
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from concurrent.futures import as_completed, ProcessPoolExecutor
import argparse
import multiprocessing


def head(raw, size=1_000_000):
    return io.BytesIO(raw.read(size))

@dataclass
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
        bytes_read =  len(buf.getbuffer())
        zfh = gzip.open(buf)
        n = 0
        recsizes = 0
        readlens = 0
        try:
            for hdr, seq, qhdr, qual in zip(zfh, zfh, zfh, zfh):
                n += 1
                recsizes += len(hdr) + len(seq) + len(qhdr) + len(qual)
                readlens += len(seq) -1 
        except EOFError:
            pass
        fsize = fq.stat().st_size
        estim_reads = fsize / bytes_read * n
        return FQStat(path=fq, file_size=fsize, estimated_nreads=round(estim_reads),
                      mean_read_len=readlens/n, mean_record_size=recsizes/n,
                      n_reads_sampled=n, bytes_per_record=bytes_read/n)



def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", "-o", type=argparse.FileType("w"), default=stdout,
        help="Output table")
    ap.add_argument("--threads", "-j", type=int, default=multiprocessing.cpu_count(),
        help="Parallel CPUs")
    ap.add_argument("--head", "-b", type=int, default=1_000_000,
        help="Inspect the first N bytes")
    ap.add_argument("fastqs", nargs="+")

    args = ap.parse_args(argv)

    results = []
    with ProcessPoolExecutor(args.threads) as exc:
        jobs = set()
        for fq in args.fastqs:
            jobs.add(exc.submit(estimate_fq_stats, fq, head_bytes=args.head))
        for res in tqdm(as_completed(jobs)):
            results.append(res.result())

    print("path", "file_size", "estimated_n_reads", "read_length", "record_size", "bytes_per_record", sep="\t", file=args.out)
    for res in results:
        print(res.path, res.file_size, res.estimated_nreads, res.mean_read_len, res.mean_record_size, res.bytes_per_record, sep="\t", file=args.out)

