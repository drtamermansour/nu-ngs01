### Download the sequencing data (if it is not already downloaded)

```bash
mkdir -p ~/workdir/sample_data && cd ~/workdir/sample_data

if [ ! -f HBR_UHR_ERCC_ds_5pc.tar ];then 
  wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
  tar -xvf HBR_UHR_ERCC_ds_5pc.tar
else
  echo "the TAR file already exist"
fi
```

### Install [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

```bash
conda activate ngs1
conda install -c bioconda -y trinity
conda install -c bioconda -y bowtie  ## bowtie is a short read aligner that triniy uses 
```

### [Run Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)

```bash
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

### What is SS_lib_type?

```
--SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.

Paired reads:
    *  RF: first read (/1) of fragment pair is sequenced as anti-sense (reverse(R)), and second read (/2) is in the sense strand (forward(F)); typical of the dUTP/UDG sequencing method.
    *  FR: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) is in the antisense strand (reverse)
    
Unpaired (single) reads:
    *  F: the single read is in the sense (forward) orientation
    *  R: the single read is in the antisense (reverse) orientation
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

```bash
TrinityStats.pl trinity_out_dir/Trinity.fasta
```
