from tqdm.auto import tqdm
from collections import defaultdict
import re
from sys import stdout, stderr, stdin
from ._utils import fqparse



def parse_header(header):
    name, *comment = header.split(" ", 1)
    m = re.match(r"@([^ ]+)/([12])", name)
    if m:
        return header, m[1], m[2], ""
    if comment:
        if comment[0].startswith("1:"):
            return header, name, "1", comment[0]
        if comment[0].startswith("2:"):
            return header, name, "2", comment[0]
    raise ValueError(f"Unable to parse fastq header '{header}'")

def pairslash(header):
    _, name, i, comment = parse_header(header)
    if comment:
        return f"{name}/{i} {comment}"
    return f"{name}/{i}"

def check_pair(r1, r2):
    _, n1, p1, c1 = parse_header(r1[0])
    _, n2, p2, c2 = parse_header(r2[0])
    if n1 != n2:
        raise ValueError(f"Bad fastq pair: {r1}, {r2}")

def read_fileset(readfiles):
    res = defaultdict(dict)
    for readfile in readfiles:
        stem, read = re.match(r"(.+?)[_.]?(R[12]|SE|IL|INTERLEAVED)?(_001)?\.(FQ|FASTQ)(\.GZ)?", readfile, re.I).groups()[:2]
        if read is not None:
            read = read.upper()
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


def handle_fileset(args, stem, fileset):
    if "R1" in fileset and "R2" in fileset:
        r1s = fqparse(fileset["R1"])
        r2s = fqparse(fileset["R2"])
        for r1, r2 in tqdm(zip(r1s, r2s), desc=f"{stem} R1/R2"):
            check_pair(r1, r2)
            output_read(args, r1)
            output_read(args, r2)
        if not (next(r1s, None) is None and next(r2s, None) is None):
            raise ValueError(f"Uneven number of R1/R2 reads in {fileset['R1']} and {fileset['R2']}")
    if "IL" in fileset:
        for r in tqdm(fqparse(fileset["IL"]), desc=f"{stem} IL"):
            output_read(args, r)
    if "SE" in fileset:
        for r in tqdm(fqparse(fileset["SE"]), desc=f"{stem} SE"):
            output_read(args, r)
    if None in fileset:
        for r in tqdm(fqparse(fileset[None]), desc=stem):
            output_read(args, r)
            

def main(argv=None):
    """Handle paired-end reads, with various transformations."""
    import argparse
    ap = argparse.ArgumentParser("blsl pairs")
    ap.add_argument("--output", "-o", default=stdout, type=argparse.FileType("w"),
        help="Output file (default: stdout)")
    ap.add_argument("--pairslash", "-s", action="store_true",
        help="Add /1 /2 to paired reads")
    ap.add_argument("readfile", nargs='+', default="/dev/stdin")
    args=ap.parse_args(argv)

    for stem, fileset in read_fileset(args.readfile).items():
        handle_fileset(args, stem, fileset)

if __name__ == "__main__":
    main()

