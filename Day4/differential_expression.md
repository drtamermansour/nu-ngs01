# Differential Expression

## What is the common outcome of an RNA-Seq analysis?
The goal of most RNA-Seq analyses is to find genes or transcripts that change across experimental conditions.
This change is called **differential expression**. By finding these genes and transcripts, we can infer the functional characteristics of the different conditions.

## Why using RNA in differential expression?
Unlike DNA, which is static, the mRNA abundances change over time.
You will need to ensure not only that you observed a change but that this change is correlated with the experimental conditions.
Typically this is achieved by measuring the same state multiple times.

## What is a "spike-in" control?
The goal of the spike-in control is to determine how well we can measure and reproduce data with known (expected) properties.
A common product called the "ERCC ExFold RNA Spike-In Control Mix" can be added in different mixtures.
This spike-in consists of 92 transcripts that are present in known
concentrations across a wide abundance range (from very few copies to many copies).

---

## Prepare the data

### Download

```sh
wget -c https://transfer.sh/ExXNk/data.tar.gz
tar xvzf data.tar.gz
```

### Does this data have "spike-in" control?
Yes there are two mixes: ERCC Mix 1 and ERCC Mix2. The spike-in consists of 92 transcripts that are present in known concentrations across a wide abundance range (from very few copies to many copies).
[More info](http://tools.thermofisher.com/content/sfs/manuals/cms_086340.pdf)

### Description

The data consists of two commercially available RNA samples:

1. Universal Human Reference **(UHR)** is total RNA isolated from a diverse set of 10 cancer cell lines. [more info](https://www.chem-agilent.com/pdf/strata/740000.pdf)
2. Human Brain Reference **(HBR)** is total RNA isolated from the brains of 23 Caucasians, male and female, of varying age but mostly 60-80 years old. [more info](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/sp_6052.pdf)

And the Human Chromosome 22 will be used as a reference.

#### The data was produced in three replicates for each condition.

1. UHR + ERCC Mix1, Replicate 1, **HBR_1**
2. UHR + ERCC Mix1, Replicate 2, **HBR_2**
3. UHR + ERCC Mix1, Replicate 3, **HBR_3**
4. HBR + ERCC Mix2, Replicate 1, **UHR_1**
5. HBR + ERCC Mix2, Replicate 2, **UHR_2**
6. HBR + ERCC Mix2, Replicate 3, **UHR_3**

`tree data` to see the folders structure.

---

## Setup enviornemnt

```bash
conda create -n df_exp
conda activate df_exp
conda install hisat2
conda install kallisto
conda install samtools
conta install subread
```


---

## Analyzing control samples

### Alignment: Hisat2

#### Step 1 (Indexing)

```bash
REF_HUMAN=data/refs/22.fa
IDX_HUMAN=data/refs/22.fa

REF_ERCC=data/refs/ERCC92.fa
IDX_ERCC=data/refs/ERCC92.fa

GTF_ERCC=data/refs/ERCC92.gtf
GTF_HUMAN=data/refs/22.gtf

hisat2-build $REF_ERCC $IDX_ERCC
hisat2-build $REF_HUMAN $IDX_HUMAN
```

#### Step 2 (Alignment)

```bash
IDX=data/refs/ERCC92.fa
RUNLOG=runlog.txt

for SAMPLE in HBR UHR;
do
    for REPLICATE in 1 2 3;
    do
        R1=data/reads/${SAMPLE}_${REPLICATE}_R1.fq
        R2=data/reads/${SAMPLE}_${REPLICATE}_R2.fq
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        hisat2 $IDX -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
    done
done
```

#### Step 3 (Quantifying)

```bash
GTF=data/refs/ERCC92.gtf

# Generate the counts.
featureCounts -a $GTF -g gene_name -o counts.txt  bam/HBR*.bam  bam/UHR*.bam

# Simplify the file to keep only the count columns.
cat counts.txt | cut -f 1,7-12 > simple_counts.txt
```

> Head

| Geneid     | bam/HBR_1.bam | bam/HBR_2.bam | bam/HBR_3.bam | bam/UHR_1.bam | bam/UHR_2.bam | bam/UHR_3.bam | 
|------------|---------------|---------------|---------------|---------------|---------------|---------------| 
| ERCC-00002 | 37892         | 47258         | 42234         | 39986         | 25978         | 33998         | 
| ERCC-00003 | 2904          | 3170          | 3038          | 3488          | 2202          | 2680          | 
| ERCC-00004 | 910           | 1078          | 996           | 9200          | 6678          | 7396          | 
| ERCC-00009 | 638           | 778           | 708           | 1384          | 954           | 1108          | 
| ERCC-00012 | 0             | 0             | 0             | 2             | 0             | 0             | 
| ERCC-00013 | 0             | 0             | 0             | 4             | 4             | 0             | 
| ERCC-00014 | 26            | 20            | 8             | 20            | 4             | 16            | 
| ERCC-00016 | 0             | 0             | 0             | 0             | 0             | 0             | 
| ERCC-00017 | 0             | 0             | 0             | 0             | 0             | 2             | 

---
