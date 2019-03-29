## Assessing the Read Content of the Transcriptome Assembly


### 1. Index Building + Alignment

```bash
cd ~/workdir/trinity

# Build Index
bowtie2-build trinity_out_dir/Trinity.fasta Trinity.fa

# Run Alignment
bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fa -1 $R1 -2 $R2 2> align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam

# View align_stats.txt
cat 2>&1 align_stats.txt

```

### 3. Visualize read support using IGV

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