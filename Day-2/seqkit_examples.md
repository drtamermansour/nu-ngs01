# Seqit


## 1- Read and print

### Read from file

`seqkit seq hairpin.fa.gz`

### Read from stdin

`zcat hairpin.fa.gz | seqkit seq`

## 2- Sequence Types

> By default, seqkit seq automatically detect the sequence type

```bash
echo -e ">seq\nacgtryswkmbdhvACGTRYSWKMBDHV" | seqkit stats
echo -e ">seq\nACGUN ACGUN" | seqkit stats
echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit stats
echo -e "@read\nACTGCN\n+\n@IICCG" | seqkit stats
seqkit stat hairpin.fa.gz
```

### Force to read specific file type by `-t` or `--seq-type`

```bash
echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit seq -t dna
```

## 3- Print names

### Full name

`seqkit seq hairpin.fa.gz -n`

### Only IDs

`seqkit seq hairpin.fa.gz -n -i`

## 4- Only print seq (global flag -w defines the output line width, 0 for no wrap)

`seqkit seq hairpin.fa.gz -s -w 0`


## 5- Reverse comlement sequence

`seqkit seq hairpin.fa.gz -r -p`

## 6- Remove gaps and to lower/upper case

`echo -e ">seq\nACGT-ACTGC-ACC" | seqkit seq -g -u`

## 7- Convert RNA to DNA

`echo -e ">seq\nUCAUAUGCUUGUCUCAAAGAUUA" | seqkit seq --rna2dna`

## 8- Filter by sequence length

```bash
zcat hairpin.fa.gz | seqkit seq | seqkit stats
zcat hairpin.fa.gz | seqkit seq -m 100 | seqkit stats
zcat hairpin.fa.gz | seqkit seq -m 100 -M 1000 | seqkit stats
```

---

## GREP

### 1- Extract human hairpins (i.e. sequences with name starting with hsa)

`zcat hairpin.fa.gz | seqkit grep -r -p ^hsa`

### 2- Remove human and mice hairpins.

`zcat hairpin.fa.gz | seqkit grep -r -p ^hsa -p ^mmu -v`

### 3- Extract sequences containing AGGCG

`zcat hairpin.fa.gz | seqkit grep -s -i -p aggcg`

### 4- Extract sequences starting with AGGCG

`zcat hairpin.fa.gz | seqkit grep -s -r -i -p ^aggcg`

### 5- zcat hairpin.fa.gz | seqkit grep -s -d -i -p TTSAA

`zcat hairpin.fa.gz | seqkit grep -s -d -i -p TTSAA`
`zcat hairpin.fa.gz | seqkit grep -s -r -i -p TT[CG]AA`

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
