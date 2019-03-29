### Download Data
```bash
cd ~/workdir/sample_data
wget https://transfer.sh/z8zqO/chr22_with_ERCC92.fa.gz
gunzip chr22_with_ERCC92.fa.gz
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf
#wget https://transfer.sh/IbpI7/HBR_UHR_ERCC_ds_5pc.tar
#tar -xvf HBR_UHR_ERCC_ds_5pc.tar
```

### install [Hisat](https://ccb.jhu.edu/software/hisat2/index.shtml)
```
source activate ngs1
conda install -c bioconda hisat2 
```

### Indexing
```
mkdir -p ~/workdir/hisat_align/hisatIndex && cd ~/workdir/hisat_align/hisatIndex
ln -s ~/workdir/sample_data/chr22_with_ERCC92.fa .
hisat2_extract_splice_sites.py ~/workdir/sample_data/chr22_with_ERCC92.gtf > splicesites.tsv
hisat2_extract_exons.py ~/workdir/sample_data/chr22_with_ERCC92.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_with_ERCC92
```

### Align
```
cd ~/workdir/hisat_align
R1="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz"
R2="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"
hisat2 -p 1 -x hisatIndex/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S UHR_Rep1.sam

## --dta/--downstream-transcriptome-assembly
## Report alignments tailored for transcript assemblers including StringTie. 
## With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. 
## This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

##--rna-strandness <string>
## Specify strand-specific information: the default is unstranded.
## For single-end reads, use F or R for forward and reverse 
## For paired-end reads, use either FR or RF. (RF:fr-firststrand corresponds & FR:fr-secondstrand)
## With this option, every read alignment will have an XS attribute tag: '+' means a read belongs to a transcript on '+' strand of genome. '-' means a read belongs to a transcript on '-' strand of genome.
```
Check the Alignment summary!

```
118571 reads; of these:
  118571 (100.00%) were paired; of these:
    487 (0.41%) aligned concordantly 0 times
    117360 (98.98%) aligned concordantly exactly 1 time
    724 (0.61%) aligned concordantly >1 times
    ----
    487 pairs aligned concordantly 0 times; of these:
      331 (67.97%) aligned discordantly 1 time
    ----
    156 pairs aligned 0 times concordantly or discordantly; of these:
      312 mates make up the pairs; of these:
        189 (60.58%) aligned 0 times
        113 (36.22%) aligned exactly 1 time
        10 (3.21%) aligned >1 times
99.92% overall alignment rate
```
### Prepare the SAM file for assembly
```
# install Samtools
source activate ngs1
conda install samtools
# convert the SAM file into BAM file 
samtools view -bS UHR_Rep1.sam > UHR_Rep1.bam
#convert the BAM file to a sorted BAM file. 
samtools sort UHR_Rep1.bam -o UHR_Rep1.sorted.bam
```

### install stringtie
```
source activate ngs1
conda install stringtie
```

### Assembly without known annotations 
```
stringtie UHR_Rep1.sorted.bam --rf -l ref_free -o ref_free.gtf
## how many transcript do you have?
cat ref_free.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l
```

### Assembly with known previous annotations 
```
stringtie UHR_Rep1.sorted.bam --rf -l ref_sup -G ~/workdir/sample_data/chr22_with_ERCC92.gtf -o ref_sup.gtf 
## how many transcript do you have?
cat ref_sup.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l
```


<hr>

