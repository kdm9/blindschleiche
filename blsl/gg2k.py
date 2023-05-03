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


def gg2k_main(argv=None):
    """Summarise a table with GreenGenes-style lineages into a kraken-style report."""
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--column", required=True, type=lambda x: int(x) - 1,
            help="Which column contains lineages?")
    ap.add_argument("-d", "--column-delim", default="\t",
            help="Quality score for glue sequences")
    ap.add_argument("-D", "--lineage-delim", default=";",
            help="Delimiter between nodes in a taxon lineage")
    ap.add_argument("-L", "--max-level", default=None, type=lambda x: int(x) - 1,
            help="Maximum level to descend to.")
    ap.add_argument("table", help="Table of taxon hits")
    args = ap.parse_args(argv)
    
    
    class OurDialect(csv.unix_dialect):
        delimiter = args.column_delim
        quotechar = '"'
        doublequote = True
        skipinitialspace = True


    lineage = recursive_defaultdict()
    with open(args.table) as fh:
        for rec in csv.reader(fh, dialect=OurDialect):
            linstr = rec[args.column]
            lin = linstr.split(args.lineage_delim)
            X = lineage
            X["__count__"] += 1
            X["__name__"] = "root"
            X["__line__"] = "NA"
            for i, t in enumerate(lin):
                X = X[t]
                X["__count__"] += 1
                X["__name__"] = t
                X["__line__"] = linstr
                if args.max_level is not None and i >= args.max_level:
                    break

    n = int(math.ceil(math.log10(lineage["__count__"])))
    for level, name, count, linstr in recursive_taxtable(lineage):
        linstr = linstr.split(args.lineage_delim, level)
        linstr = args.lineage_delim.join(linstr[:-1])
        print(str(count).rjust(n + level * 2), name.ljust(15), level, linstr, sep="|")

if __name__ == "__main__":
    gg2k_main()

