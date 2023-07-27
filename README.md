# kreeq
k-mer analysis quality evaluation tool

**kreeq** is a single fast kmer-based QV computation tool for genome assemblies. It computes both the kmers and the QV on the fly. The kmer counting capabilities are inherited from [kcount](https://github.com/vgl-hub/kcount).

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/kreeq.git --recursive` and `make -j` in the `kreeq` folder.

## Usage

```
kreeq validate -f input.[fasta|fastq][.gz] -r reads1.fastq[.gz] reads2.fastq[.gz] [...] [-k 21]
```

It accepts multiple read files as input, separated by space. To check out all options and flags use `kreeq -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
kreeq validate -f testFiles/random1.fasta -r testFiles/random1.fastq
```

Importantly, the kreeq database can only be computed once on the read set, and reused for multiple analyses to save runtime:

```
kreeq validate -r testFiles/random1.fastq -o db.kreeq
kreeq validate -f testFiles/random1.fasta -d db.kreeq
```

Similarly, kreeq databases can be generated separately for multiple inputs and combined, with increased performance in HPC environments:

```
kreeq validate -r testFiles/random1.fastq -o random1.kreeq
kreeq validate -r testFiles/random2.fastq -o random2.kreeq

time kreeq union -d random1.kreeq random2.kreeq -o union.kreeq
time kreeq validate -f testFiles/random1.fasta -d union.kreeq
```

## How to cite

This tool is part of the **gfastar** tool suite. If you use **kreeq** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
