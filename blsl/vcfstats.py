#!/usr/bin/env python3
from tqdm import tqdm
from cyvcf2 import VCF
import pandas as pd

import sys
from sys import stdin, stdout, stderr
from subprocess import Popen, PIPE
import argparse
import math
from concurrent.futures import ProcessPoolExecutor, as_completed

def parallel_regions(vcf, cores=1):
    V = VCF(vcf)
    # 10 chunks per chrom per core
    chunks = len(V.seqlens)*10*cores
    for cname, clen in zip(V.seqnames, V.seqlens):
        chunk = int(math.ceil(clen/chunks))
        for start in range(0, clen, chunk):
            s = start+1
            e = start+chunk+1
            yield f"{cname}:{s}-{e}"

def variant2dict(v):
    for i, alt in enumerate(v.ALT):
        dat = {"CHROM": v.CHROM, "POS": v.POS, "REF": v.REF, "ALT": alt, "QUAL": v.QUAL}
        dat["call_rate"] = v.call_rate
        for key, val  in v.INFO:
            if isinstance(val, str) and "," in val:
                val = val.split(',')
            if isinstance(val, tuple) or isinstance(val, list):
                val = val[i]
            dat[f"INFO_{key}"] = val
        yield dat

def bcftools_info_with_tags(vbcf):
    res = []
    cmd=f"bcftools +fill-tags {vbcf} -Ou -- -d -t all,F_MISSING"
    with Popen(cmd, shell=True, stdout=PIPE) as proc:
        for v in tqdm(VCF(proc.stdout)):
            for r in variant2dict(v):
                res.append(r)
    return res

def one_chunk_stats(vcf, chunk, fill=True):
    cmd=f"bcftools view -r {chunk} {vcf} -Ou"
    if fill:
        cmd = f"{cmd} | bcftools +fill-tags - -Ou -- -d -t all,F_MISSING"
    res = []
    with Popen(cmd, shell=True, stdout=PIPE) as proc:
        for v in VCF(proc.stdout):
            for r in variant2dict(v):
                res.append(r)
    return res

def chunkwise_bcfools_stats(args):
    with ProcessPoolExecutor(args.threads) as exc:
        jobs = []
        for region in parallel_regions(args.vcf):
            jobs.append(exc.submit(one_chunk_stats, args.vcf, region, fill=args.fill_tags_first))
        for job in tqdm(as_completed(jobs), total=len(jobs), unit="chunk"):
            for res in job.result()
                if args.fields:
                    res = {k: dat[k] for k in fields if k in dat}
                print(json.dumps(res), file=args.output)


def main(argv=None):
    """Use bcftools to calculate various statistics, outputing an R-ready table"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default=stdout, type=argparse.FileType("wt"),
            help="Output jsonlines file")
    ap.add_argument("-f", "--fields", nargs="+",
            help="Use only these fields")
    ap.add_argument("-F", "--fill-tags-first", action="store_true",
            help="Use bcftools +fill-tags to pre-fill more fields (forces no parallelism)")
    ap.add_argument("-t", "--threads", default=1, type=int,
            help="Number of parallel threads")

    ap.add_argument("vcf")
    args = ap.parse_args(argv)

    chunkwise_bcfools_stats(args)


if __name__ == "__main__":
    main()
