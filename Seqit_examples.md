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

## 5- Convert multi-line FASTQ to 4-line FASTQ

`seqkit seq reads_1.fq.gz -w 0`

## 6- Reverse comlement sequence

`seqkit seq hairpin.fa.gz -r -p`

## 7- Remove gaps and to lower/upper case

`echo -e ">seq\nACGT-ACTGC-ACC" | seqkit seq -g -u`

## 8- Convert RNA to DNA

`echo -e ">seq\nUCAUAUGCUUGUCUCAAAGAUUA" | seqkit seq --rna2dna`

## 9- Filter by sequence length

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



















