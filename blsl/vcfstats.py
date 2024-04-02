#!/usr/bin/env python3
from tqdm import tqdm
from cyvcf2 import VCF

import sys
from sys import stdin, stdout, stderr
from subprocess import Popen, PIPE
import argparse
import math
import json
from concurrent.futures import ProcessPoolExecutor, as_completed

def parallel_regions(vcf, cores=1):
    V = VCF(vcf)
    # 10 chunks per chrom per core
    chunks = len(V.seqlens)*100*cores
    for cname, clen in zip(V.seqnames, V.seqlens):
        chunk = 10000#int(math.ceil(clen/chunks))
        for start in range(0, clen, chunk):
            s = start+1
            e = start+chunk+1
            yield f"{cname}:{s}-{e}"

def variant2dict(v, fields=None):
    for i, alt in enumerate(v.ALT):
        dat = {"CHROM": v.CHROM, "POS": v.POS, "REF": v.REF, "ALT": alt, "QUAL": v.QUAL}
        dat["call_rate"] = v.call_rate
        for key, val  in v.INFO:
            if isinstance(val, str) and "," in val:
                val = val.split(',')
            if isinstance(val, tuple) or isinstance(val, list):
                val = val[i]
            dat[f"INFO_{key}"] = val
        if fields:
            dat = {K:V for K, V in dat.items() if K in fields}
        yield json.dumps(dat)

def one_chunk_stats(vcf, chunk, fill=True, fields=None):
    cmd=f"bcftools view -r {chunk} {vcf} -Ou"
    if fill:
        cmd = f"{cmd} | bcftools +fill-tags - -Ou -- -d -t all,F_MISSING"
    res = []
    with Popen(cmd, shell=True, stdout=PIPE) as proc:
        vcf = VCF(proc.stdout)
        for v in vcf:
            for r in variant2dict(v, fields=fields):
                res.append(r)
        del vcf
    return res

def chunkwise_bcfools_stats(args):
    regions = list(parallel_regions(args.vcf))
    with tqdm(total=len(regions), unit="chunk") as pbar:
        for i in range(0, len(regions), 10000):
            to=min(len(regions), i+10000)
            with ProcessPoolExecutor(args.threads) as exc:
                jobs = (exc.submit(one_chunk_stats, args.vcf, region, fill=args.fill_tags_first, fields=args.fields) for region in regions[i:to])
                for job in as_completed(jobs):
                    pbar.update(1)
                    for res in job.result():
                        args.output.write(res)
                        args.output.write("\n")


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
