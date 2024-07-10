from sys import stderr
from .pairs import main as pairs_main

def main(argv=None):
    """Add an old-style /1 /2 pair indicator to paired-end fastq files (DEPRECATED, use blsl pairs)"""
    print("WARNING: blsl pairslash is deprecated, replace use with blsl pairs", file=stderr)
    pairs_main(argv)

if __name__ == "__main__":
    main()

