#!/usr/bin/env python3
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(*args, **kwargs):
        return args[0]

import argparse
import re
from pathlib import Path
from sys import stderr


def make_replacer(reffa, mode, delim="~", outdelim="~", refdelim="~"):
    names = []
    if not reffa.endswith(".fai"):
        reffa += ".fai"
    with open(reffa) as fh:
        for line in fh:
            names.append(line.split()[0])
    fromto = {}
    for name in names:
        splitname = name.strip().split(refdelim)
        fromname = delim.join(splitname)
        if mode == "cat": 
            fromto[fromname] = outdelim.join(splitname)
        elif mode == "add":
            fromto[splitname[2]] = outdelim.join(splitname)
        elif mode == "rm":
            fromto[fromname] = splitname[2]
    return fromto


def do_tsv(infh, outfh, cols, replacer, coldelim="\t"):
    for line in tqdm(infh):
        if line.startswith("#"):
            outfh.write(line)
            continue
        fields = line.rstrip("\r\n").split(coldelim)
        for col in cols:
            fields[col] = replacer[fields[col]]
        print(*fields, sep=coldelim, file=outfh)


def do_fasta(infh, outfh, replacer):
    for line in infh:
        if not line.startswith(">"):
            outfh.write(line)
            continue
        name = line.lstrip(">").rstrip("\n\r").split(" ")[0]
        name = replacer[name]
        print(">", name, sep="", file=outfh)


fasta_exts = {".fa", ".fna", ".fasta"}
tsv_exts = {".gff", ".gff3", ".gtf", ".bed", ".bg", ".tsv"}


def main(argv=None):
    """Add, remove, or modify PanSN-style prefixes to contig/chromosome names in references"""
    extra="""
    Modes:

        cat: pass through, ensuring all names are PanSN-style. Useful to convert delimiters in PanSN names.

        rm: remove the PanSN specific prefix, e.g. at9900_1_chr1 -> chr1
        
        add: add a PanSN-specific prefix, e.g. chr1 -> at9900_1_chr1.
    """
    ap = argparse.ArgumentParser("blsl pansn-rename", epilog=extra)
    ap.add_argument("-C", "--colsep", default="\t",
        help="Delimiter between columns in c/tsv file. (default tab)")
    ap.add_argument("-d", "--delim", default="_",
        help="Delimiter between fields in PanSN names in <INPUT> (the # in indiv#1#chr1)")
    ap.add_argument("-D", "--out-delim",
        help="--delim for the output. Defaults to the same as --delim.")
    ap.add_argument("-r", "--reference-fasta", required=True,
        help="Reference fasta file with PanSN names.")
    ap.add_argument("-R", "--ref-delim",
        help="Delimiter between fields in PanSN names in the --reference-fasta file (default: same as --delim)")
    ap.add_argument("-c", "--colidx", action="append", type=int,
        help="Which column in TSV should be sub'd? (default depends on input "
             "filetype, normally 0, use more than once for multiple columns "
             "e.g. -c 0 -c 3 for 1st & 4th cols).")
    ap.add_argument("-m", "--mode", choices=["cat", "rm", "add"],
        help="What should I do? (see below)")
    ap.add_argument("-o", "--output", default="-", type=argparse.FileType("w"),
        help="output filename, default: stdout")
    ap.add_argument("input", type=argparse.FileType("r"))
    args = ap.parse_args(argv)

    if args.colidx is None:
        args.colidx = [0, ]
    args.colidx = list(sorted(set(args.colidx)))
    if args.out_delim is None:
        args.out_delim = args.delim
    if args.ref_delim is None:
        args.ref_delim = args.delim
        
    
    replacer = make_replacer(args.reference_fasta, args.mode, args.delim, args.out_delim, args.ref_delim)
    
    input_ftype = "tsv"
    try:
        input_filename = Path(args.input.name)
        if input_filename.suffix in fasta_exts:
            input_ftype = "fasta"
        elif input_filename.suffix in tsv_exts:
            input_ftype =  "tsv"
        else:
            pass
            #print("WARNING: unknown filetype", input_filename.suffix, ", assuming TSV-like format", file=stderr)
    except:
        print("WARNING: unable to determine input filetype, assuming TSV-like format", file=stderr)
    
    if input_ftype == "tsv":
        do_tsv(args.input, args.output, args.colidx, replacer, coldelim=args.colsep)
    else:
        do_fasta(args.input, args.output, replacer)


if __name__ == "__main__":
    main()
