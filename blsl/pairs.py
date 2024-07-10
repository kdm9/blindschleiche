from tqdm.auto import tqdm
from collections import defaultdict
import re
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


def pairslash(header):
    name, *comment = header.split(" ", 1)
    if name.endswith("/1") or name.endswith("/2"):
        return header
    i = seqnum(name, comment)
    return f"{name}/{i} {comment[0]}"


def read_fileset(readfiles):
    res = defaultdict(dict)
    for readfile in readfiles:
        stem, read = re.match(r"(.+?)[_.]?(R[12]|SE|IL|INTERLEAVED)?(_001)?\.(FQ|FASTQ)(\.GZ)?", readfile.upper()).groups()[:2]
        if read == "INTERLEAVED":
            read = "IL"
        if read in res[stem]:
            raise ValueError("Two files with same stem and read index: {readfiles}")
        res[stem][read] = readfile
    return dict(res)


def output_read(args, read):
    if args.pairslash:
        read[0] = pairslash(read[0])
    print("\n".join(read), file=args.output)


def handle_fileset(args, fileset):
    if "R1" in fileset and "R2" in fileset:
        for r1, r2 in tqdm(zip(fqparse(fileset["R1"]), fqparse(fileset["R2"]))):
            output_read(args, r1)
            output_read(args, r2)
    if "IL" in fileset:
        for r in tqdm(fqparse(fileset("IL"))):
            output_read(args, r)
    if "SE" in fileset:
        for r in tqdm(fqparse(fileset("SE"))):
            output_read(args, r)
            

def main(argv=None):
    """Handle paired-end reads, with various transformations."""
    import argparse
    ap = argparse.ArgumentParser("blsl pairs")
    ap.add_argument("--output", "-o", default="/dev/stdout", 
        help="Output file (default: stdout)")
    ap.add_argument("--pairslash", "-s", action="store_true",
        help="Add /1 /2 to paired reads")
    ap.add_argument("readfile", nargs='+', default="/dev/stdin")
    args=ap.parse_args(argv)

    for stem, fileset in read_fileset(args.readfile).items():
        handle_fileset(args, fileset)

if __name__ == "__main__":
    main()

