Sequence Alignment
==================


## Download reference file
```
cd ~/workdir/sample_data
wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz

OR

wget https://transfer.sh/7g7zo/dog.fa.gz -O dog_chr5.fa.gz

gunzip dog_chr5.fa.gz
```

BLAST Alignment
===============

## install [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs)
```
source activate ngs1
conda install -c bioconda blast 
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
head -n500 sample.fa > sample_subset.fa
/usr/bin/time -v blastn -db blastIndex/dog_chr5.fa -query sample_subset.fa -outfmt 5 -out BD143_TGACCA_L005.xml
```

## convert to SAM
```
conda install -c conda-forge biopython 
wget https://gist.githubusercontent.com/ozagordi/099bdb796507da8d9426/raw/6ca66616fd545fbb15d94b079e46a7c55edb54c0/blast2sam.py
python blast2sam.py BD143_TGACCA_L005.xml > BD143_TGACCA_L005.sam

```

BWA Alignment
=============

## install [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)
```
source activate ngs1
conda install -c bioconda bwa 
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
```
cd ~/workdir/sample_data
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz
gunzip gencode.v29.pc_transcripts.fa.gz
cat gencode.v29.pc_transcripts.fa | awk -F'|' '{print $1}' > gencode.v29.pc_transcripts.fa.simplified.fa

wget https://transfer.sh/IbpI7/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar

# Download the Transcriptome Annotation File

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz

# subset the data

head -n 36189 gencode.v29.pc_transcripts.fa.simplified.fa > gencode.v29.pc_transcripts.1000.fa

```


## install [Kallisto](https://pachterlab.github.io/kallisto/)

> For generating BAM files, must use Kallisto's latest release.

```
##### Download
wget -c https://github.com/pachterlab/kallisto/releases/download/v0.45.1/kallisto_linux-v0.45.1.tar.gz

##### Extract and move to usr/local/bin

tar xvzf kallisto*gz
sudo cp kallisto_linux-v0.45.1/kallisto /usr/local/bin/
```

###  Run Indexing
```
mkdir -p ~/workdir/kallisto_align/kallistoIndex && cd ~/workdir/kallisto_align/kallistoIndex
ln -s ~/workdir/sample_data/gencode.v29.pc_transcripts.1000.fa.
kallisto index -i human_pc.idx -k 25 gencode.v29.pc_transcripts.1000.fa
```

### Run PseudoAlignment to generate PseudoBAM file

```
cd ~/workdir/kallisto_align
R1="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz"
R2="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"
kallisto quant -i kallistoIndex/human_pc.idx -o human_pc_bam $R1 $R2 --pseudobam

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


[Download](http://software.broadinstitute.org/software/igv/download)

