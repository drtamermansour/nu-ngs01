## HISAT2


### Download Data

```bash
mkdir -p hisat/reference
cd hisat/reference/
wget https://transfer.sh/z8zqO/chr22_with_ERCC92.fa.gz
wget https://transfer.sh/HS0na/chr22_with_ERCC92.gtf
gunzip *gz
```
### Annotation
```
cd hisat/reference/
less -p start_codon -S chr22_with_ERCC92.gtf
grep ENST00000342247 chr22_with_ERCC92.gtf | less -p "exon\s" -S
```
### Indexing
```
~/hisat2-2.1.0/hisat2_extract_splice_sites.py ~/hisat/reference/chr22_with_ERCC92.gtf > ~/hisat/reference/splicesites.tsv
~/hisat2-2.1.0/hisat2_extract_exons.py ~/hisat/reference/chr22_with_ERCC92.gtf > ~/hisat/reference/exons.tsv
~/hisat2-2.1.0/hisat2-build -p 8 --ss ~/hisat/reference/splicesites.tsv --exon ~/hisat/reference/exons.tsv ~/hisat/reference/chr22_with_ERCC92.fa ~/hisat/reference/chr22_with_ERCC92
```
### Downloading Data
```
mkdir -p ~/hisat/data
cd ~/hisat/data
wget https://transfer.sh/IbpI7/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
ls
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -P "^\@HWI" | wc -l
```
### Align
```
mkdir -p ~/hisat/align
#hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
~/hisat2-2.1.0/hisat2 -p 8 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x ~/hisat/reference/chr22_with_ERCC92 --dta --rna-strandness RF -1 ~/hisat/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 ~/hisat/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep1.sam
```
