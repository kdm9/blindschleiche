import gzip

from ._utils import fqparse

def seqnum(n, c):
    if n.endswith("/1"):
        return 1
    if n.endswith("/2"):
        return 2
    if c:
        if c[0].startswith("1:"):
            return 1
        if c[0].startswith("2:"):
            return 2
    raise ValueError(f"Unknown pair from fastq header: '{n} {''.join(c)}'")

def main(argv=None):
    """Add an old-style /1 /2 pair indicator to paired-end fastq files"""

    import argparse
    ap = argparse.ArgumentParser("blsl pairslash")
    ap.add_argument("readfile", nargs='?', default="/dev/fd/0")
    args=ap.parse_args(argv)


    if args.readfile.endswith(".gz"):
        reads = gzip.open(args.readfile, "rt")
    else:
        reads = open(args.readfile, "r")

    for s in fqparse(reads):
        n, *c = s[0].split(" ", 1)
        if n.endswith("/1") or n.endswith("/2"):
            print(s[0])
        else:
            i = seqnum(n, c)
            print(f"{n}/{i} {''.join(c)}")
        print(s[1], s[2], s[3], sep="\n")

if __name__ == "__main__":
    main()

