# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import argparse
import random
import time
import os
from tqdm import tqdm
try:
    import orjson as json
except ImportError:
    import json

def main(argv=None):
    """Parse jsonlines into a C/TSV"""
    ap = argparse.ArgumentParser()
    ap.add_argument("--fraction", "-f", type=float,
            help="Subsample this fraction of rows (0-1)")
    ap.add_argument("--seed", type=float, default=time.time() + os.getpid(),
            help="Seed the RNG with this seed (seeded adequately by default)")
    ap.add_argument("--header-same", "-s", action="store_true",
            help="Assume all rows have the same header, just use the first row's header")
    ap.add_argument("--out-csv", "-o", type=argparse.FileType("wt"), default="/dev/stdout",
            help="Output file")
    ap.add_argument("in_json")

    args=ap.parse_args(argv)
    header = set()
    first = True

    rand = random.Random(args.seed)
    with open(args.in_json) as fh:
        for line in tqdm(fh, desc="Read Headers", unit="lines"):
            if first or args.fraction is None or rand.random() < args.fraction:
                dat = json.loads(line)
                for key in dat:
                    header.add(key)
                if first and args.header_same:
                    break
                first = False

    try:
        sep = "\t" if args.out_csv.filename.endswith("tsv") else ","
    except:
        sep = ","
    rand = random.Random(args.seed)
    print(*header, sep=sep, file=args.out_csv)
    with open(args.in_json) as fh:
        for line in tqdm(fh, desc="Read Data", unit="lines"):
            if args.fraction is None or rand.random() < args.fraction:
                dat = json.loads(line)
                outline = [dat.get(key, "") for key in header]
                print(*outline, sep=sep, file=args.out_csv)
                

if __name__ == "__main__":
    main()
