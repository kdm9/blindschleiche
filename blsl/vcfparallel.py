#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from tqdm import tqdm
from cyvcf2 import VCF

from sys import stdin, stdout, stderr
import shutil
import subprocess
import argparse
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import math
import uuid
import re


def parallel_regions(vcf, cores=1, npercore=10):
    V = VCF(vcf)
    # 10 chunks per chrom per core
    chunks = len(V.seqlens)*npercore*cores
    for cname, clen in zip(V.seqnames, V.seqlens):
        chunk = int(max(min(clen, 1000), math.ceil(clen/chunks)))
        for start in range(0, clen, chunk):
            s = start+1
            e = min(clen, start+chunk+1)
            yield f"{cname}:{s:09d}-{e:09d}"


def one_chunk(vcf, chunk, temp_prefix, filterarg="", pipeline="", index=True):
    ofile = f"{temp_prefix}{chunk}.bcf"
    cmd=f"bcftools view -r {chunk} {filterarg} -Ou {vcf}"
    if pipeline:
        cmd = f"{cmd} | {pipeline}"
    if index:
        cmd = f"{cmd} | bcftools view -Ob0 --write-index -o {ofile}"
    cmd = f"({cmd}) >{ofile} 2>{ofile}.log"
    subprocess.run(cmd, shell=True, check=True)
    return ofile


def chunkwise_pipeline(args):
    filestomerge = []
    with ProcessPoolExecutor(args.threads) as exc:
        jobs = []
        if args.regions is None:
            regions = parallel_regions(args.vcf)
        else:
            regions = set()
            for line in args.regions:
                chrom, start0, end, *_ = line.rstrip("\n").split("\t")
                start0 = int(start0)
                regions.add(f"{chrom}:{start0+1}-{end}")
            regions = list(sorted(regions))
        for region in regions:
            jobs.append(exc.submit(one_chunk, args.vcf, region, args.temp_prefix, filterarg=args.filter, pipeline=args.commands, index=not args.merge_with_cat))
        for job in tqdm(as_completed(jobs), total=len(jobs), unit="chunk", desc="chunkwise variant pipeline"):
            ofile = job.result()
            filestomerge.append(ofile)
    return filestomerge


def merge_one(files, prefix, threads=1, merge_type="fast"):
    fofn = f"{prefix}fofn.txt"
    output = f"{prefix}output.bcf"
    with open(fofn, "w") as fh:
        for file in sorted(files):
            print(file, file=fh)
    merge = "--allow-overlaps --rm-dup all" if merge_type == "slow" else ""
    cmd = f"(bcftools concat --file-list {fofn} {merge} --threads {threads} -Ob0 --write-index -o {output}) >{output}.log 2>&1"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        with open(output + ".log") as fh:
            print(fh.read())
        raise
    return output


def merge_results(args, filestomerge):
    if args.merge_with_cat:
        with open(args.output, "wb") as ofh:
            for file in tqdm(filestomerge, desc="merge with cat", unit="file"):
                with open(file, "rb") as fh:
                    shutil.copyfileobj(fh, ofh)
        return args.output
    filestomerge =list(filestomerge)
    if len(filestomerge) > 1000:
        n_groups = int(math.ceil(math.sqrt(len(filestomerge))))
        print(f"Grouped merge -> {n_groups} groups")
        groupsize = int(math.ceil(len(filestomerge) / n_groups))
        groups = [filestomerge[i:min(i+groupsize, len(filestomerge))] for i in range(0, len(filestomerge), groupsize)]
        final_merge = []
        with ProcessPoolExecutor(args.threads) as exc:
            jobs = []
            for i, files in enumerate(groups):
                jobs.append(exc.submit(merge_one, files, f"{args.temp_prefix}merge_group_{i}_", args.merge_type))
            for job in tqdm(as_completed(jobs), total=len(jobs), unit="group"):
                ofile = job.result()
                final_merge.append(ofile)
    else:
        final_merge = filestomerge

    fofn = f"{args.temp_prefix}final_fofn.txt"
    with open(fofn, "w") as fh:
        for file in sorted(final_merge):
            print(file, file=fh)
    index = "--write-index" if re.match(r"[zb]", args.outformat) else ""
    merge = "--allow-overlaps --rm-dup all" if args.merge_type == "slow" else ""
    print(f"Finally merging {len(final_merge)} files to {args.output}")
    cmd = f"bcftools concat --file-list {fofn} {merge} --threads {args.threads} -O{args.outformat} -o {args.output} {index} >{args.output}.log 2>&1"
    subprocess.run(cmd, shell=True, check=True)
    return args.output


def main(argv=None):
    """Use bcftools to calculate various statistics, outputing an R-ready table"""
    ap = argparse.ArgumentParser()
    ap.add_argument("-R", "--regions", type=argparse.FileType("rt"),
            help="Use regions defined in this BED file. NB! Will only do these regions and will not check for overlaps. Be careful. Cols should be chrom, start0, end, ... and tab separated (i.e. a BED file).")
    ap.add_argument("-o", "--output", required=True,
            help="Output V/BCF file")
    ap.add_argument("-O", "--outformat", default="z8",
            help="Output file format (passed directly to bcftools -O, see man bcftools but z=vcf.gz and b=bcf)")
    ap.add_argument("-T", "--temp-prefix", default=None,
            help="Temporary output directory")
    ap.add_argument("-S", "--slow-merge", default="fast", action="store_const", const="slow", dest="merge_type",
            help="Use slow bcftools merging (with --allow-overlaps and --remove-duplicates)")
    ap.add_argument("-f", "--filter", default="", type=str,
            help="bcftools view arguments for variant filtering")
    ap.add_argument("-c", "--commands", default="", type=str,
            help="command(s) to operate. Must take uncompressed bcf on stdin and yield bcf (i.e -Ob0) on stdout. Can use | and other shell features.")
    ap.add_argument("-M", "--merge-with-cat", action="store_true",
            help="command(s) produce a non-BCF file, so merge with 'cat'.")
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
