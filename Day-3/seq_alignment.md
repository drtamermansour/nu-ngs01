Sequence Alignment
==================


## Download reference file
```
mkdir -p ~/workdir/sample_data
cd ~/workdir/sample_data

# Use one of the 2 links for download
# Link 1
wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz

# OR
# Link 2
wget https://github.com/drtamermansour/nu-ngs01/raw/master/Day-3/sample_data/dog_chr5.fa.gz
---

gunzip dog_chr5.fa.gz

```

BLAST Alignment
===============

## install [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs)
```
conda activate ngs1
conda install -c bioconda -y blast
```

## index your genome
```
mkdir -p ~/workdir/blast_align/blastIndex && cd ~/workdir/blast_align/blastIndex
ln -s ~/workdir/sample_data/dog_chr5.fa .
makeblastdb -in dog_chr5.fa -dbtype nucl
```

## sequence alignment
```
cd ~/workdir/blast_align
R1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.fq.gz"
R2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.fq.gz"
zcat $R1 $R2 | paste - - - - | awk 'BEGIN{OFS="\n";}{print ">"$1,$2}' > sample.fa
/usr/bin/time -v blastn -db blastIndex/dog_chr5.fa -query sample.fa -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send" -out BD143_TGACCA_L005.blastOut

## try the xml format
head -n500 sample.fa | sed 's/^>@/>/' > sample_subset.fa
/usr/bin/time -v blastn -db blastIndex/dog_chr5.fa -query sample_subset.fa -outfmt 5 -out BD143_TGACCA_L005.xml
```

## convert to SAM
```
conda install -c conda-forge -y biopython
wget https://gist.githubusercontent.com/ozagordi/099bdb796507da8d9426/raw/6ca66616fd545fbb15d94b079e46a7c55edb54c0/blast2sam.py
python blast2sam.py BD143_TGACCA_L005.xml > BD143_TGACCA_L005.sam

# install Samtools
conda install -y samtools

# convert SAM file to BAM
samtools view -S -b BD143_TGACCA_L005.sam -o BD143_TGACCA_L005.bam

# Sorting the BAM file
samtools sort BD143_TGACCA_L005.bam -o sorted_BD143_TGACCA_L005.bam

# Indexing the BAM file
samtools index sorted_BD143_TGACCA_L005.bam


```

BWA Alignment
=============

## install [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)
```
conda activate ngs1
conda install -c bioconda -y bwa
```

## index your genome

```
mkdir -p ~/workdir/bwa_align/bwaIndex && cd ~/workdir/bwa_align/bwaIndex
ln -s ~/workdir/sample_data/dog_chr5.fa .
bwa index -a bwtsw dog_chr5.fa
```

## sequence alignment

```
cd ~/workdir/bwa_align
R1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.fq.gz"
R2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.fq.gz"
/usr/bin/time -v bwa mem bwaIndex/dog_chr5.fa $R1 $R2 > BD143_TGACCA_L005.sam
```


Kallisto
========

## Download reference & sample files

```bash
cd ~/workdir/sample_data

# Download human transcriptome
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz
# OR
wget -c https://github.com/drtamermansour/nu-ngs01/raw/master/Day-3/sample_data/gencode.v29.pc_transcripts.fa.gz

gunzip gencode.v29.pc_transcripts.fa.gz

wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

tar -xvf HBR_UHR_ERCC_ds_5pc.tar

# Download the Transcriptome Annotation File
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
# OR
wget -c https://github.com/drtamermansour/nu-ngs01/raw/master/Day-3/sample_data/gencode.v29.annotation.gtf.gz

gunzip gencode.v29.annotation.gtf.gz

# Select the transcripts of Chr22
cd ~/workdir/sample_data/
conda activate ngs1
READS=$(grep "^chr22" gencode.v29.annotation.gtf | awk -F'\t' '{print $9}' | awk -F';' '{print $1}' | awk -F' ' '{print $2}' | awk -F'"' '{print $2}' | sort | uniq)

for value in $READS
    do  
        echo "Processing: $value"
        seqkit grep -r -p ${value} gencode.v29.pc_transcripts.fa | awk -F'|' '{print $1}' >> gencode.v29.pc_transcripts.chr22.simplified.fa
    done
```


## install [Kallisto](https://pachterlab.github.io/kallisto/)

> For generating BAM files, must use Kallisto's latest release.

```bash
conda install -c bioconda -y kallisto
```

###  Run Indexing
```
mkdir -p ~/workdir/kallisto_align/kallistoIndex && cd ~/workdir/kallisto_align/kallistoIndex
ln -s ~/workdir/sample_data/gencode.v29.pc_transcripts.chr22.simplified.fa .
kallisto index -i human_pc.idx -k 25 gencode.v29.pc_transcripts.chr22.simplified.fa
```

### Run PseudoAlignment to generate PseudoBAM file

```
cd ~/workdir/kallisto_align
R1="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz"
R2="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"
kallisto quant -i kallistoIndex/human_pc.idx -o human_pc_bam $R1 $R2 --pseudobam

# install Samtools
conda install -y samtools

# Sorting the BAM file
samtools sort human_pc_bam/pseudoalignments.bam -o human_pc_bam/sorted_pseudoalignments.bam

# Indexing the BAM file
samtools index human_pc_bam/sorted_pseudoalignments.bam

```

### Run PseudoAlignment to generate PseudoBAM file with Coordinates

```
ln -s ~/workdir/sample_data/gencode.v29.annotation.gtf .
kallisto quant -i kallistoIndex/human_pc.idx -o human_pc_bam_gtf $R1 $R2 --genomebam --gtf gencode.v29.annotation.gtf

# Sorting and indexing done automatically in the last step

```

#### Visualization of human_pc_bam_gtf/*.bam

![](https://github.com/mr-eyes/nu-ngs01/blob/master/Day-3/pseudobam-GTF.png)

### IGV

[IGV Tutorial](https://bioinformatics-ca.github.io/resources/IGV_Tutorial.pdf)


[Open IGV in the web browser](https://igv.org/app/)


### OR run IGV locally (unrecommended)
[Download](http://software.broadinstitute.org/software/igv/download)

```
cd ~
wget http://data.broadinstitute.org/igv/projects/downloads/2.5/IGV_Linux_2.5.0.zip
unzip IGV_Linux_2.5.0.zip
sudo echo 'export IGV=$HOME/IGV_Linux_2.5.0/igv.sh' >> ~/.bashrc
source ~/.bashrc
source activate ngs1
bash $IGV -g ../sample_data/gencode.v29.pc_transcripts.chr22.simplified.fa human_pc_bam/sorted_pseudoalignments.bam
```
