#!/usr/bin/python3
from shutil import get_terminal_size
from collections import Counter
from sys import stdin, stdout, stderr, exit
from math import ceil, log10
import select

def main(argv=None):
    """Make a summary histogram of git-annex-list output"""
    if not select.select([stdin], [], [], 1)[0]:
        print("USAGE:  git annex list ... | blsl galhist", file=stderr)
        print("This tool makes a histogram of file occurence patterns to summarise the output of git annex list. Just pipe git anex list in on stdin.", file=stderr)
        exit(1)

    twidth = 80
    try:
        if stdout.isatty():
            twidth = get_terminal_size()[0]
    except:
        pass
        
    header = []
    counts = Counter()
    for line in stdin:
        if not header:
            header.append(line.rstrip())
        elif line.startswith("|"):
            header.append(line.rstrip())
        else:
            key, file = line.split(" ", maxsplit=1)
            counts[key]+=1 

    for l in header:
        print(l)
    _, mcv = counts.most_common(1)[0]
    nw = int(ceil(log10(mcv)))
    remaining = twidth - len(header) - 1 - nw - 2 

    for k, v in counts.most_common():
        vbarw = min(round(v/mcv*remaining), remaining)
        vbar = ("#" * vbarw).ljust(remaining, " ")
        npad = str(v).rjust(nw, " ")
        print(f"{k} {vbar} ({npad})")
