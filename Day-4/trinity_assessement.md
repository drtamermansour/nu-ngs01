# Transcriptome Assembly Quality Assessment

## Assessing the Read Content of the Transcriptome Assembly

### 1. Index Building + Alignment

```bash
source activate ngs1

mkdir -p ~/workdir/trinity/trinity_out_dir/bowtie2_assessment/index  && cd ~/workdir/trinity/trinity_out_dir/bowtie2_assessment 

# Build Index
bowtie2-build ~/workdir/trinity/trinity_out_dir/Trinity.fasta index/Trinity.fa

# Run Alignment

R1="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz"
R2="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"

bowtie2 -p 1 -q --no-unal -k 20 -x index/Trinity.fa -1 $R1 -2 $R2 2> align_stats.txt| samtools view -Sb -o bowtie2.bam

# View align_stats.txt
cat align_stats.txt

# calc the alignments per transcript
samtools view bowtie2.bam | awk '{print $3}' | sort | uniq -c | sort -nr > alignment_per_transcript.count 

#### Some help
-q                 query input files are FASTQ .fq/.fastq (default)
-p 				         number of threads
--no-unal          suppress SAM records for unaligned reads
2> 				         pipe the stderr
```


### 2. Visualize read support using IGV

```bash

# Sorting the BAM file
samtools sort bowtie2.bam -o bowtie2.coordSorted.bam

# Indexing the BAM file
samtools index bowtie2.coordSorted.bam

# Index Trinity.fasta
samtools faidx ~/workdir/trinity/trinity_out_dir/Trinity.fasta

# Visualize
bash $IGV -g ~/workdir/trinity/trinity_out_dir/Trinity.fasta  bowtie2.coordSorted.bam
```

