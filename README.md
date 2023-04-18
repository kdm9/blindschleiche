# blindschleiche

Misc sequence tools in python. These tools are small things which I have at some point needed to write because I couldn't find a solution I liked. This is in no way a comprehensive toolkit. It is a companion to [seqhax](https://github.com/kdm9/seqhax), execpt that seqhax is written in C/C++ and generally contains tools to handle very large datasets where performance is somewhat important. This is all in python for ease of development, and so typically these tools perform less data- or compute-intensive tasks.

## Install

```
pip install blindschleiche
# or for the current main branch:
# pip install git+https://github.com/kdm9/blindschleiche.git
```

## Usage

```
USAGE: blsl <subtool> [options...]


Where <subtool> is one of:

  genigvjs:             Generate a simple IGV.js visualisation of some bioinf files.
  telogrep:             Search contigs for known telomere repeats
  n50:                  Calculate N50 and total length of a set of contigs
  falen:                Tabulate the lengths of sequences in a FASTA file
  mask2bed:             The inverse of bedtools maskfasta: softmasked fasta -> unmasked fasta + mask.bed
  liftoff-gff3:         Obtain an actually-useful GFF3 from Liftoff by fixing basic GFF3 format errors
  pansn-rename:         Add, remove, or modify PanSN-style prefixes to contig/chromosome names in references
  ildemux:              Demultiplex modern illumina reads from read headers.
  ilsample:             Sample a fraction of read pairs from an interleaved fastq file
  regionbed:            Make a bed/region file of genome windows
  uniref-acc2taxid:     Make a ncbi-style acc2taxid.map file for a uniref fasta
  help:                 Print this help message


Use blsl subtool --help to get help about a specific tool
```

## Why Blindschleiche

1) [They're awesome animals](https://www.google.com/search?q=blindschleiche&tbm=isch)
2) Their English name is Slow Worm, which is appropriate for this set of low-performance tools in Python.
3) All tools implemented in python must be named with a snake pun, and they're kinda a snake (not really, they're legless lizards)

