# blindschleiche

A collection of bioinformatics / sequence utilities needed for my research, and hopefully useful for yours.

[![DOI](https://zenodo.org/badge/509693094.svg)](https://zenodo.org/doi/10.5281/zenodo.10049825)

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

  deepclust2fa:         Split a .faa by the clusters diamond deepclust finds
  ebiosra2rl2s:         INTERNAL: MPI Tübingen tool. Make a runlib-to-sample map table from ebio sra files
  equalbestblast:       Output only the best blast hits.
  esearchandfetch:      Use the Entrez API to search for and download something. A CLI companion to the NCBI search box
  falen:                Tabulate the lengths of sequences in a FASTA file
  farename:             Rename sequences in a fasta file sequentially
  galhist:              Make a summary histogram of git-annex-list output
  genigvjs:             Generate a simple IGV.js visualisation of some bioinf files.
  gffcat:               Concatenate GFF3 files, resepcting header lines and FASTA sections
  gg2k:                 Summarise a table with GreenGenes-style lineages into a kraken-style report.
  ildemux:              Demultiplex modern illumina reads from read headers.
  ilsample:             Sample a fraction of read pairs from an interleaved fastq file
  liftoff-gff3:         Obtain an actually-useful GFF3 from Liftoff by fixing basic GFF3 format errors
  mask2bed:             The inverse of bedtools maskfasta: softmasked fasta -> unmasked fasta + mask.bed
  n50:                  Calculate N50 and total length of a set of contigs
  nstitch:              Combine R1 + R2 into single sequences, with an N in the middle
  pairslash:            Add an old-style /1 /2 pair indicator to paired-end fastq files
  pansn-rename:         Add, remove, or modify PanSN-style prefixes to contig/chromosome names in references
  regionbed:            Make a bed/region file of genome windows
  shannon-entropy:      Calculate Shannon's entropy (in bits) at each column of one or more alignments
  tabcat:               Concatenate table (c/tsv) files, adding the filename as a column
  telogrep:             Search contigs for known telomere repeats
  uniref-acc2taxid:     Make a ncbi-style acc2taxid.map file for a uniref fasta
  help:                 Print this help message


Use blsl subtool --help to get help about a specific tool
```

## Why the name Blindschleiche?

1) [They're awesome animals](https://www.google.com/search?q=blindschleiche&tbm=isch)
2) Their English name is Slow Worm, which is appropriate for this set of low-performance tools in Python.
3) All tools implemented in Python must be named with a snake pun, and they're kinda a snake (not really, they're legless lizards)

