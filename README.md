# kreeq
k-mer analysis quality evaluation tool

**kreeq** is a single fast kmer-based QV computation tool for genome assemblies. It computes both the kmers and the QV on the fly. The kmer counting capabilities are inherited from (https://github.com/vgl-hub/kcount) kcount 

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/kreeq.git --recursive` and `make -j` in the `kreeq` folder.

## Usage

```
kreeq validate -f input.[fasta|fastq][.gz] -r reads.fastq[.gz] [-k 21]
```

It accepts multiple read files as input, separated by space. To check out all options and flags use `kreeq -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
kreeq validate -f testFiles/random1.fasta testFiles/random1.fastq
```

## How to cite

This tool is part of the **gfastar** tool suite. If you use **kreeq** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
