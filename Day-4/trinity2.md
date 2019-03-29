# Transcriptome Assembly Quality Assessment

## Assessing the Read Content of the Transcriptome Assembly

### 1. Index Building + Alignment

```bash
conda activate ngs1
conda install bwa

cd ~/workdir/trinity

# Build Index
bowtie2-build trinity_out_dir/Trinity.fasta Trinity.fa

# Run Alignment

R1="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz"
R2="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"

bowtie2 -p 1 -q --no-unal -k 20 -x Trinity.fa -1 $R1 -2 $R2 2> align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam

# View align_stats.txt
cat 2>&1 align_stats.txt

```

### 2. Visualize read support using IGV

```bash

# Sorting the BAM file
samtools sort bowtie2.bam -o bowtie2.coordSorted.bam

# Indexing the BAM file
samtools index bowtie2.coordSorted.bam

# Index Trinity.fasta file
samtools faidx Trinity.fasta

# Index Trinity.fasta
samtools faidx trinity_out_dir/Trinity.fasta

# Visualize
# Set $IGV to the igv.sh file
bash $IGV -g trinity_out_dir/Trinity.fasta  bowtie2.coordSorted.bam
```

## Counting Full Length Trinity Transcripts

