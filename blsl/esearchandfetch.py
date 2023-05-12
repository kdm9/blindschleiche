#!/usr/bin/env python3
from Bio import Entrez
from tqdm import tqdm

import argparse
from time import sleep
from sys import stdout, stdin, stderr, exit
from pathlib import Path
import io

def esearch(*args, chunks=10000, **kwargs):
    kwargs["retmax"] = chunks
    ids = []
    n = None
    i = 0
    with tqdm(desc="esearch IDs") as progress:
        while n is None or i < n:
            this = Entrez.esearch(*args, **kwargs, retstart=i)
            for x in range(3):
                try:
                    sleep(0.5)
                    this = Entrez.read(this)
                    break
                except RuntimeError:
                    sleep(0.5)
                    pass
            else:
                this = Entrez.read(this)
            nret = len(this["IdList"])
            ids.extend(this["IdList"])
            i += nret
            n = int(this["Count"])
            progress.total = n
            progress.update(nret)
    return ids
            

def efetch(db="nucleotide", rettype="genbank", retmode="text", ids=None, chunks=2000):
    with tqdm(desc=f"efetch {rettype}", total=len(ids)) as progress:
        for start in range(0, len(ids), chunks):
            stop = min(start + chunks, len(ids))
            idstr = ",".join(ids[start:stop])
            res = Entrez.efetch(db=db, retmode=retmode, rettype=rettype, id=idstr)
            if start != 0 and rettype=="runinfo":
                # Skip header on non-first one
                hdr = next(res)
            for line in res:
                if isinstance(line, bytes):
                    line=line.decode('utf-8')
                yield line
            progress.update(stop-start)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("-e", "--email", type=str,
            help="Email for Entrez API")
    ap.add_argument("-c", "--chunksize", type=int, default=1000,
            help="Number of records to download per chunk")
    ap.add_argument("-o", "--out", type=argparse.FileType('w'), default="-",
            help="Genbank output")
    ap.add_argument("-d", "--db", required=True,
            help="db api param")
    ap.add_argument("-r", "--rettype", required=True,
            help="rettype api param")
    ap.add_argument("-m", "--retmode", default="text",
            help="API return mode")
    ap.add_argument("-i", "--id-file", type=Path,
            help="List of IDs (give with --term to save a list of IDs to this file, or without to use this list of IDs)")
    ap.add_argument("-t", "--term", type=str)
    args = ap.parse_args(argv)

    Entrez.email = args.email

    if args.term:
        ids = esearch(term=args.term, db=args.db)
        if args.id_file:
            with args.id_file.open("w") as fh:
                for id in ids:
                    print(id, file=fh)
    elif args.id_file:
        ids = []
        with args.id_file.open("r") as fh:
            for line in fh:
                ids.append(line.rstrip())
    else:
        print("ERROR: must give either --term to search for, or a list of ids in --id-file", file=stderr)
        exit(1)

    for line in efetch(ids=ids, db=args.db, rettype=args.rettype, retmode=args.retmode, chunks=1000):
        args.out.write(line)


if __name__ == "__main__":
    main()
