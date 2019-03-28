### Install [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

```
source activate ngs1
conda install -c bioconda trinity 
```

### [Run Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)
```
mkdir -p ~/workdir/trinity && cd ~/workdir/trinity

R1="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz"
R2="$HOME/workdir/sample_data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"

Trinity --seqType fq --max_memory 2G  \
        --left $R1 \
        --right $R2 \
        --SS_lib_type RF \
        --CPU 1 \
        --quality_trimming_params "ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 MINLEN:25"
```

### Explaining the identifiers
```
>TRINITY_DN1000|c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]
 AATCTTTTTTGGTATTGGCAGTACTGTGCTCTGGGTAGTGATTAGGGCAAAAGAAGACAC
 ACAATAAAGAACCAGGTGTTAGACGTCAGCAAGTCAAGGCCTTGGTTCTCAGCAGACAGA
 AGACAGCCCTTCTCAATCCTCATCCCTTCCCTGAACAGACATGTCTTCTGCAAGCTTCTC
 CAAGTCAGTTGTTCACAGGAACATCATCAGAATAAATTTGAAATTATGATTAGTATCTGA
 TAAAGCA
 
The accession encodes the Trinity 'gene' and 'isoform' information. 
In the example above, the accession 'TRINITY_DN1000|c115_g5_i1' indicates: 
a) Trinity read cluster 'TRINITY_DN1000|c115', 
b) gene 'g5', and 
c) isoform 'i1'. 
Because a given run of trinity involves many clusters of reads, each of which are assembled separately, and because the 'gene' numberings are unique within a given processed read cluster, the 'gene' identifier should be considered an aggregate of the read cluster and corresponding gene identifier, which in this case would be 'TRINITY_DN1000|c115_g5'.
```

### Explore the assembly 
```
TrinityStats.pl trinity_out_dir/Trinity.fasta
```
