
from tqdm import tqdm

import re
import argparse
from sys import stdin, stdout, stderr
from pathlib import Path
from collections import Counter, defaultdict
import json


def as_rtype(rtype, val):
    try:
        return rtype(val)
    except ValueError:
        pass
    return None

class PatternMatcher:
    def match(self, file):
        if hasattr(self, "firstline"):
            with open(file) as fh:
                line = next(fh)
                return re.match(self.firstline, line) is not None
        if hasattr(self, "anyline"):
            with open(file) as fh:
                for line in fh:
                    if re.match(self.anyline, line) is not None:
                        return True
                return False

    
    def read_patterns(self, file):
        def extract(line):
            for name, (rtype, regex) in self.patterns.items():
                m = re.match(regex, line)
                if m is not None:
                    return name, as_rtype(rtype, m[1])
            return None, None
        res = {}
        with open(file) as fh:
            for line in fh:
                key, val = extract(line)
                if key is not None:
                    res[key] = val
        return res

    def read(self, file):
        res = {
                "filename": str(file),
                "tool": self.name,
        }
        if hasattr(self, "patterns"):
            res.update(self.read_patterns(file))
        return res


class RSEMLog(PatternMatcher):
    name="RSEM"
    anyline = r"\d+ reads; of these:"
    patterns = {
            "total_reads": (int, r"(\d+) reads; of these:"),
            "paired_reads": (int, r"\s*(\d+) \(.+%\) were paired; of these:"),
            "unaligned_reads": (int, r"\s*(\d+) \(.+%\) aligned concordantly 0 times"),
            "mapped_reads": (int, r"\s*(\d+) \(.+%\) aligned concordantly exactly 1 time"),
            "multimapped_reads": (int, r"\s*(\d+) \(.+%\) aligned concordantly >1 times"),
            "overall_alignment_rate": (float, r"\s*([\d.]+)% overall alignment rate"),
    }


class AdapterRemoval(PatternMatcher):
    name="AdapterRemoval"
    anyline = r"AdapterRemoval ver. .+"
    patterns = {
            "total_reads": (int, r"Total number of read pairs: (\d+)"),
            "reads_with_adapter": (int, r"Number of reads with adapters[1]: (\d+)"),
            "reads_kept": (int, r"Number of retained reads: (\d+)"),
            "bases_kept": (int, r"Number of retained nucleotides: (\d+)"),
            "avergage_length_postqc": (float, r"Average length of retained reads: ([\d.]+)"),
    }



def main(argv=None):
    """What if MultiQC was in the tidyverse? (and much worse)"""
    tools = {RSEMLog(), AdapterRemoval()}
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", type=argparse.FileType("wt"), default=stdout,
            help="Output file (jsonl)")
    ap.add_argument("files", nargs="+", type=Path)

    args = ap.parse_args(argv)

    results = defaultdict(list)
    toolcounts = Counter()

    for infile in tqdm(args.files):
        for proc in tools:
            if proc.match(infile):
                results[proc.name].append(proc.read(infile))
                toolcounts[proc.name] += 1

    for tool in results:
        n = toolcounts[tool]
        print(f"Found {n} {tool} files", file=stderr)
        for res in results[tool]:
            json.dump(res, args.out)
            args.out.write("\n")

if __name__ == "__main__":
    main()
