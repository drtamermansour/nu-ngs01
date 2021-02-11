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

### Download the data (if it is not already downloaded. You should have done this for denovo assembly)

```sh
mkdir -p ~/workdir/sample_data && cd ~/workdir/sample_data

if [ ! -f HBR_UHR_ERCC_ds_5pc.tar ];then 
  wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
  tar -xvf HBR_UHR_ERCC_ds_5pc.tar
else
  echo "the TAR file already exist"
fi

# Save the path for your samples in a variable to use later
READS_DIR=~/workdir/sample_data/

```

### Does this data have "spike-in" control?
Yes there are two mixes: ERCC Mix 1 and ERCC Mix2. The spike-in consists of 92 transcripts that are present in known concentrations across a wide abundance range (from very few copies to many copies).
[More info](http://tools.thermofisher.com/content/sfs/manuals/cms_086340.pdf)

### Description

The data consists of 3 replicates from each of two commercially available RNA samples:

1. Universal Human Reference **(UHR)** is total RNA isolated from a diverse set of 10 cancer cell lines. [more info](https://www.chem-agilent.com/pdf/strata/740000.pdf)
2. Human Brain Reference **(HBR)** is total RNA isolated from the brains of 23 Caucasians, male and female, of varying age but mostly 60-80 years old. [more info](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/sp_6052.pdf)


#### The data was produced in three replicates for each condition.

1. UHR + ERCC Mix1, Replicate 1, **HBR_1**
2. UHR + ERCC Mix1, Replicate 2, **HBR_2**
3. UHR + ERCC Mix1, Replicate 3, **HBR_3**
4. HBR + ERCC Mix2, Replicate 1, **UHR_1**
5. HBR + ERCC Mix2, Replicate 2, **UHR_2**
6. HBR + ERCC Mix2, Replicate 3, **UHR_3**

---

## Setup enviornemnt for alignment and quantification

```bash
# For alignment and quantification: we will try 2 options: 
# A) Genome-based alignment by Hisat2 then quantification by featureCounts 
# B) Transcriptome-based psudo-alignment and quantification by kallisto


conda activate ngs1
# A) Genome-based alignment by Hisat2 then quantification by featureCounts

# Install Hisat2 and samtools 
conda install -c bioconda -y hisat2
conda install -y samtools
# Install subread, we will use featureCount : a software program developed for counting reads to genomic features such as genes, exons, promoters and genomic bins.
conda install subread
# Download the reference and its GTF annotation (if they are not already downloaded. You should have done this for reference-based assembly)
cd ~/workdir/sample_data
if [ ! -f chr22_with_ERCC92.fa ];then 
  wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
else
  echo "the chr22_with_ERCC92.fa file already exist"
fi

if [ ! -f chr22_with_ERCC92.gtf ];then 
  wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf
else
  echo "the chr22_with_ERCC92.gtf file already exist"
fi
# Save the path for the reference and annotation in variables to use later
REF=~/workdir/sample_data/chr22_with_ERCC92.fa
GTF=~/workdir/sample_data/chr22_with_ERCC92.gtf

# B) Transcriptome-based psudo-alignment and quantification by kallisto

# Install Kallisto
conda install -c bioconda -y kallisto

# Kallisto needs a reference transcriptome. Transform the GTF/REF into fasta file
conda install -c bioconda gffread
gffread chr22_with_ERCC92.gtf -g chr22_with_ERCC92.fa -w chr22_with_ERCC92_transcripts.fasta 
```

## Setup enviornemnt for differential expression & visualization

```bash
# For differential expression, we will use DESeq R package and for visualization, we will use gplots package. 
conda install r
conda install -y bioconductor-deseq r-gplots
https://raw.githubusercontent.com/drtamermansour/nu-ngs01/master/Day-6/deseq1.r
https://raw.githubusercontent.com/drtamermansour/nu-ngs01/master/Day-6/draw-heatmap.r

```


---


## A) Genome-based alignment pipeline  (Hisat2 for alignmnet, featureCounts for quantification and DESeq for DE) 

```bash
## Hisat2 Indexing
## You should have generated this index already before (if not, see how to do this in the reference-based assembly lecture)
hisatIndex=~/workdir/hisat_align/hisatIndex/chr22_with_ERCC92


## Alignment
mkdir -p ~/workdir/diff_exp && cd ~/workdir/diff_exp/
mkdir -p bams
for SAMPLE in HBR UHR;
do
    for REPLICATE in 1 2 3;
    do
        R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*read1.fastq.gz
        R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*read2.fastq.gz
        BAM=bams/${SAMPLE}_${REPLICATE}.bam

        hisat2 $hisatIndex -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
    done
done
```

> You can visualize BAM files using the [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/download)


```bash

## Quantification
featureCounts -a $GTF -g gene_name -o counts.txt  bams/HBR*.bam  bams/UHR*.bam

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


```bash
## Differential expression by DESeq1
cat simple_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.tsv
```

> Head DESeq1 Output

| id         | baseMean         | baseMeanA        | baseMeanB        | foldChange       | log2FoldChange   | pval                 | padj                 | 
|------------|------------------|------------------|------------------|------------------|------------------|----------------------|----------------------| 
| ERCC-00130 | 29681.8244237545 | 10455.9218232761 | 48907.7270242329 | 4.67751460376823 | 2.22574215774208 | 1.16729711209977e-88 | 9.10491747437823e-87 | 
| ERCC-00108 | 808.597670575459 | 264.877838024487 | 1352.31750312643 | 5.10543846632202 | 2.35203486825767 | 2.40956154792562e-62 | 9.39729003690993e-61 | 
| ERCC-00136 | 1898.3382995277  | 615.744918976546 | 3180.93168007886 | 5.16598932779828 | 2.36904466305553 | 2.80841619396469e-58 | 7.30188210430821e-57 | 
| ERCC-00116 | 952.57953992746  | 337.704944218003 | 1567.45413563692 | 4.64149004174798 | 2.21458802318734 | 1.72224091670521e-45 | 3.35836978757517e-44 | 



#### DeSEQ1 Output header description

- `id`: Gene or transcript name that the differential expression is computed for,
- `baseMean`: The average normalized value across all samples,
- `baseMeanA, baseMeanB`: The average normalized gene expression for each condition,
- `foldChange`: The ratio baseMeanB/baseMeanA ,
- `log2FoldChange`: log2 transform of foldChange . When we apply a 2-based logarithm the values
become symmetrical around 0. A log2 fold change of 1 means a doubling of the expression level, a log2 fold change of -1 shows show a halving of the expression level.
- `pval`: The probability that this effect is observed by chance,
- `padj`: The adjusted probability that this effect is observed by chance.

#### View only rows with pval < 0.05

```bash
cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > hisat_output.pdf
```

---

## B) Transcriptome-based psudo-alignment pipeline  (Kallisto for psudo-alignmnet & quantification and DESeq for DE) 

## What is Kallisto?
Kallisto is a software package for quantifying transcript abundances.
The tool perform a pseudoalignment of reads against a transcriptome, In pseudoalignment, the program tries to identify for each read the target that it originates from but not where in the target it aligns.
This makes the algorithm much faster than a 'real' alignment algorithm.

## Automate the same expirement with Kallisto

```bash

## Kallisto indexing
mkdir -p ~/workdir/kallisto_align/kallistoIndex && cd ~/workdir/kallisto_align/kallistoIndex
ln -s ~/workdir/sample_data/chr22_with_ERCC92_transcripts.fasta .
kallisto index -i chr22_with_ERCC92_transcripts.idx -k 25 chr22_with_ERCC92_transcripts.fasta


## Psudo-alignment and quantification
kallistoIndex=~/workdir/kallisto_align/kallistoIndex/chr22_with_ERCC92_transcripts.idx 
mkdir -p ~/workdir/diff_exp && cd ~/workdir/diff_exp/
mkdir -p quants
for SAMPLE in HBR UHR;
do
    for REPLICATE in 1 2 3;
    do
        # Build the name of the files (Paired End).
        R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*read1.fastq.gz
        R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*read2.fastq.gz

        echo "*** Running kallisto on ${SAMPLE}_${REPLICATE}"
        OUT=quants/${SAMPLE}_${REPLICATE}
        kallisto quant -i $kallistoIndex -o $OUT -b 100 $R1 $R2 
    done
done

echo "*** Created counts for ERCC control samples."
for quant_file in quants/*/abundance.tsv;do
  id=$(echo $quant_file | awk -F"/" '{print $2}'); 
  echo $id
  quant_file_simplified=$(dirname $quant_file)/abundance_simplified.tsv
  echo "target_id $id" | tr ' ' '\t' > $quant_file_simplified
  tail -n+2 $quant_file | cut -f 1,4 >> $quant_file_simplified
done
paste quants/H*/abundance_simplified.tsv quants/U*/abundance_simplified.tsv | cut -f 1,2,4,6,8,10,12  > ercc_kallisto_counts.tsv

## Differential expression by DESeq1
echo "*** Running the DESeq1 and producing the final result: ercc_kallisto_deseq1.tsv"
cat ercc_kallisto_counts.tsv | Rscript deseq1.r 3x3 > ercc_kallisto_deseq1.tsv 

# Filter pval < 0.05
cat ercc_kallisto_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_kallisto_deseq1.tsv

```
> Head kallisto DeSeq1 output

| id         | baseMean         | baseMeanA        | baseMeanB        | foldChange       | log2FoldChange   | pval                 | padj                 | 
|------------|------------------|------------------|------------------|------------------|------------------|----------------------|----------------------| 
| ERCC-00130 | 14840.7854784434 | 5227.93407242622 | 24453.6368844606 | 4.67749526786055 | 2.22573619391768 | 8.73963105304599e-87 | 6.81691222137587e-85 | 
| ERCC-00108 | 404.298508194945 | 132.438264826674 | 676.158751563215 | 5.10546368489592 | 2.35204199450587 | 9.72746450192469e-59 | 3.79371115575063e-57 | 
| ERCC-00136 | 949.168332376237 | 307.870824713046 | 1590.46584003943 | 5.16601675888527 | 2.36905232365751 | 2.5211305539385e-57  | 6.55493944024009e-56 | 
| ERCC-00116 | 476.289329027567 | 168.851590236675 | 783.727067818459 | 4.64151428316386 | 2.21459555802611 | 1.93683785478913e-44 | 3.77683381683879e-43 |


---

## Visualize


```bash

cat filtered_kallisto_deseq1.tsv | Rscript draw-heatmap.r > kallisto_output.pdf

```









