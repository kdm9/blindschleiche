#!/usr/bin/env python3
import re
import json
import select
from sys import stdin, stdout, stderr, exit
import sys

def main(argv=None):
    """INTERNAL: MPI Tübingen tool. Make a runlib-to-sample map table from ebio sra files"""
    samples = {}
    if not select.select([stdin,],[],[],2.0)[0]:
        print("This is an internal MPI Tübingen tool to make a table that mapps runs of a given library to sample names, e.g. for Acanthophis.", file=stderr)
        print(f"USAGE: find /ebio/abt6_sra/... -type f -name '*.fq.gz' | ebiosra2rl2s > run-lib2sample.tsv", file=stderr)
        exit(1)
    for line in stdin:
        m = re.match(r"(.*/illumina_.+_flowcell.+_SampleId(\w+)_(RunId\d+_LaneId\d+)/.+_(R\d)_\d+.fastq.gz)", line.strip())
        if m is not None:
            path, sample, run, read = m.groups()
            if sample not in samples:
                samples[sample] = {run: {read: path}}
            elif run not in samples[sample]:
                samples[sample][run] = {read: path}
            else:
                samples[sample][run][read] = path
    #json.dump(samples, stdout, indent=2)
    print("library", "run", "read1_uri", "read2_uri", file=stdout, sep='\t')
    for sample, dat in samples.items():
        for run, reads in dat.items():
            print(sample, run, reads["R1"], reads["R2"], file=stdout, sep='\t')


if __name__ == "__main__":
    main()
