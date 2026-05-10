# SeqKit

## Install seqkit

```bash
# Activate the enviornment
conda activate ngs1

# Install seqkit
conda install -y seqkit
```

## Create a directory for the tutorial

```
mkdir -p ~/workdir/seqkit/ && cd ~/workdir/seqkit/
```

## Download sample data

```bash
wget https://www.mirbase.org/download/hairpin.fa
```

<hr>

## 1- Read and print

### Read from file

`seqkit seq hairpin.fa`

### Read from stdin

`cat hairpin.fa | seqkit seq`

## 2- Sequence Types

> By default, seqkit seq automatically detect the sequence type

```bash
echo -e ">seq\nAACGATCGAATATCGCGTAGCTACGTACGC" | seqkit stats
echo -e ">seq\nACGUN ACGUN" | seqkit stats
echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit stats
echo -e "@read\nACTGCN\n+\n@IICCG" | seqkit stats
seqkit stat hairpin.fa
```

### Force to read specific file type by `-t` or `--seq-type`

```bash
echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit seq -t dna
```

## 3- Print names

### Full name

`seqkit seq hairpin.fa -n`

### Only IDs

`seqkit seq hairpin.fa -n -i`

## 4- Only print seq (global flag -w defines the output line width, 0 for no wrap)

`seqkit seq hairpin.fa -s -w 0`


## 5- Reverse complement sequence

`seqkit seq hairpin.fa -r -p`

## 6- Remove gaps and to lower/upper case

`echo -e ">seq\nACGT-ACTGC-ACC" | seqkit seq -g -u`

## 7- Convert RNA to DNA

`echo -e ">seq\nUCAUAUGCUUGUCUCAAAGAUUA" | seqkit seq --rna2dna`

## 8- Filter by sequence length

```bash
cat hairpin.fa | seqkit seq | seqkit stats
cat hairpin.fa | seqkit seq -m 100 | seqkit stats
cat hairpin.fa | seqkit seq -m 100 -M 1000 | seqkit stats
```

---

## GREP

### 1- Extract human hairpins (i.e. sequences with name starting with hsa)

`cat hairpin.fa | seqkit grep -r -p ^hsa`

### 2- Extract human and mice hairpins.

`cat hairpin.fa | seqkit grep -r -p ^hsa -p ^mmu`

### 3- Extract sequences containing AGGCG

`cat hairpin.fa | seqkit grep -s -i -p aggcg`

### 4- Extract sequences starting with AGGCG

`cat hairpin.fa | seqkit grep -s -r -i -p ^aggcg`

---

## Sorting

### 1- sort by ID

`echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet`

### 2- sort by ID, ignoring case.

`echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -i`

### 3- sort by seq, ignoring case.

`echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -s -i`

### 4- sort by sequence length.

`echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAAnnn\n>seq3\nacgt" | seqkit sort --quiet -l`


---

## shuffle

### General Usage

`seqkit shuffle hairpin.fa`

### Shuffle with a two-pass mode (For very large files)

`seqkit shuffle -2 hairpin.fa`

---
