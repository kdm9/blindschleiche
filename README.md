# blindschleiche

Misc sequence tools in python. These tools a small things which I have at some point needed to write because I couldn't find a solution I liked. This is in no way a comprehensive toolkit. It is a companion to [seqhax](https://github.com/kdm9/seqhax), execpt that seqhax is written in C/C++ and generally contains tools to handle very large datasets where performance is somewhat important. This is all in python for ease of development, and so typically these tools perform less data- or compute-intensive tasks.

## Install

```
pip install git+https://github.com/kdm9/blindschleiche.git
```

## Usage

```
USAGE: blsl <subtool> [options...]


Where <subtool> is one of:

  telogrep:       Search contigs for known telomere repeats
  n50:            Calculate N50 and total length of a set of contigs


Use `blsl subtool --help` to get help about a specific tool
```
