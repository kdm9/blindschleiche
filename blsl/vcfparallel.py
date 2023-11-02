#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from tqdm import tqdm
from cyvcf2 import VCF

import sys
from sys import stdin, stdout, stderr
import subprocess
import argparse
import math
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import uuid


def parallel_regions(vcf, cores=1):
    V = VCF(vcf)
    # 10 chunks per chrom per core
    chunks = len(V.seqlens)*10*cores
    for cname, clen in zip(V.seqnames, V.seqlens):
        chunk = int(max(min(clen, 1000), math.ceil(clen/chunks)))
        for start in range(0, clen, chunk):
            s = start+1
            e = min(clen, start+chunk+1)
            yield f"{cname}:{s:09d}-{e:09d}"


def one_chunk(vcf, chunk, temp_prefix, filter="", pipeline=""):
    cmd=f"bcftools view -r {chunk} {filter} -Ou {vcf}"
    if pipeline:
        cmd = f"{cmd} | {pipeline}"
    ofile = f"{temp_prefix}{chunk}.bcf"
    with open(ofile, "wb") as ofh, open(ofile + ".log", "wb") as log:
        subprocess.run(cmd, shell=True, stdout=ofh, stderr=log, check=True)
        subprocess.run(f"bcftools index -f {ofile}", shell=True, stderr=log, check=True)
    return ofile


def chunkwise_pipeline(args):
    filestomerge = []
    with ProcessPoolExecutor(args.threads) as exc:
        jobs = []
        for region in parallel_regions(args.vcf):
            jobs.append(exc.submit(one_chunk, args.vcf, region, args.temp_prefix, filter=args.filter, pipeline=args.commands))
        for job in tqdm(as_completed(jobs), total=len(jobs), unit="chunk"):
            ofile = job.result()
            filestomerge.append(ofile)
    return filestomerge


def merge_results(args, filestomerge):
    fofn = f"{args.temp_prefix}fofn.txt"
    with open(fofn, "w") as fh:
        for file in sorted(filestomerge):
            print(file, file=fh)
    cmd = f"bcftools concat --file-list {fofn} --rm-dup all --allow-overlaps --threads {args.threads} -O{args.outformat} -o {args.output}"
    proc = subprocess.run(cmd, shell=True, check=True)
    return proc.returncode


def main(argv=None):
    """Use bcftools to calculate various statistics, outputing an R-ready table"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", required=True,
            help="Output V/BCF file")
    ap.add_argument("-O", "--outformat", default="z8",
            help="Output file format (passed directly to bcftools -O, see man bcftools but z=vcf.gz and b=bcf)")
    ap.add_argument("-T", "--temp-prefix", default=None,
            help="Temporary output directory")
    ap.add_argument("-f", "--filter", default="", type=str,
            help="bcftools view arguments for variant filtering")
    ap.add_argument("-c", "--commands", default="", type=str,
            help="command(s) to operate. Must take uncompressed bcf on stdin and yield uncompresed bcf on stdout. Can use | and other shell features.")
    ap.add_argument("-t", "--threads", default=1, type=int,
            help="Number of parallel threads")
    ap.add_argument("vcf")
    args = ap.parse_args(argv)

    if args.temp_prefix is None:
        args.temp_prefix = tempfile.gettempdir() + "/bcffilter"
    tp = Path(args.temp_prefix)
    if tp.is_dir() and not args.temp_prefix.endswith("/"):
        args.temp_prefix += "/"
    args.temp_prefix += str(uuid.uuid4())
    print(f"Using temporary file prefix: {args.temp_prefix}. This WILL NOT be automatically cleaned!!", file=stderr)

    chunkfiles = chunkwise_pipeline(args)
    merge_results(args, chunkfiles)

if __name__ == "__main__":
    main()
