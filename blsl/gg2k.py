# Copyright (c) 2023 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import csv
import sys
import math
import argparse
from collections import defaultdict

def recursive_defaultdict():
    x = defaultdict(recursive_defaultdict)
    x["__count__"] = 0
    return x

def recursive_taxtable(lineage, level=0):
    yield (level, lineage["__name__"], lineage["__count__"], lineage["__line__"])
    for k, v in reversed(sorted(filter(lambda x: x[0] not in ["__count__", "__name__", "__line__"],
                                       lineage.items()),
                                key=lambda x: x[1]["__count__"])):
        yield from recursive_taxtable(v, level=level+1)

def insert_lineage(tree, lin, max_level=None, n=1):
    X = tree
    X["__count__"] += n
    X["__name__"] = "root"
    X["__line__"] = "NA"
    for i, t in enumerate(lin):
        X = X[t]
        X["__count__"] += n
        X["__name__"] = t
        X["__line__"] = ";".join(lin)
        if max_level is not None and i >= max_level:
            break

def gg2k_main(argv=None):
    """Summarise a table with GreenGenes-style lineages into a kraken-style report."""
    ap = argparse.ArgumentParser()
    ap.add_argument("-l", "--lca-column", type=lambda x: int(x) - 1,
            help="Take the lowest common ancestor per value in this column")
    ap.add_argument("-f", "--column", required=True, type=lambda x: int(x) - 1,
            help="Which column contains lineages?")
    ap.add_argument("-d", "--column-delim", default="\t",
            help="Quality score for glue sequences")
    ap.add_argument("-D", "--lineage-delim", default=";",
            help="Delimiter between nodes in a taxon lineage")
    ap.add_argument("-L", "--max-level", default=None, type=lambda x: int(x) - 1,
            help="Maximum level to descend to.")
    ap.add_argument("-V", "--vsearch-size", default=None, type=lambda x: int(x) - 1,
            help="Parse vsearch size=NNNN headers from this column, use NNNN as a count rather than # hits/reads.")
    ap.add_argument("table", help="Table of taxon hits")
    args = ap.parse_args(argv)
    
    
    class OurDialect(csv.unix_dialect):
        delimiter = args.column_delim
        quotechar = '"'
        doublequote = True
        skipinitialspace = True

    lineage = recursive_defaultdict()
    with open(args.table) as fh:
        current_id = None
        lca = []
        for rec in csv.reader(fh, dialect=OurDialect):
            linstr = rec[args.column]
            lin = linstr.split(args.lineage_delim)
            n = 1
            if args.lca_column is not None:
                if current_id != rec[args.lca_column]:
                    if args.vsearch_size is not None and current_id is not None:
                        hdr = current_id
                        fields = dict(x.split("=") for x in hdr.split(";") if "=" in x)
                        n = int(fields["size"])
                    insert_lineage(lineage, lca, args.max_level, n=n)
                    lca = lin
                    current_id = rec[args.lca_column]
                for i, tax in enumerate(lin):
                    if i >= len(lca) or lca[i] != tax:
                        break
                lca = lca[:i+1]
            else:
                if args.vsearch_size is not None:
                    hdr = rec[args.vsearch_size]
                    fields = dict(x.split("=") for x in hdr.split(";") if "=" in x)
                    n = int(fields["size"])
                insert_lineage(lineage, lin, args.max_level, n=n)
        if lca:
            if args.vsearch_size is not None:
                hdr = current_id
                fields = dict(x.split("=") for x in hdr.split(";") if "=" in x)
                n = int(fields["size"])
            insert_lineage(lineage, lca, args.max_level, n=n)

    n = int(math.ceil(math.log10(lineage["__count__"])))
    taxtbl = list(recursive_taxtable(lineage))
    max_name = max(len(t[1]) for t in taxtbl)
    max_level = max(t[0] for t in taxtbl)
    for level, name, count, linstr in taxtbl:
        linstr = linstr.split(args.lineage_delim, level)
        linstr = args.lineage_delim.join(linstr[:-1])
        print(name.ljust(max_name + (max_level-level)*2).rjust(max_name + max_level * 2), str(count).rjust(n), level, linstr, sep="|")

if __name__ == "__main__":
    gg2k_main()

