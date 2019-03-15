## Download Data

```bash

mkdir -p ~/kallisto/{reference,samples}
cd ~/kallisto/

# Reference

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.pc_transcripts.fa.gz

gunzip gencode.vM20.pc_transcripts.fa.gz
mv *fa reference/

# RNA-Seq Samples
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
mv *tar samples/
cd samples/
tar xvf HBR_UHR_ERCC_ds_5pc.tar
cd ..

```

## Kallisto

### Install

```bash

##### Download
wget -c https://github.com/pachterlab/kallisto/releases/download/v0.45.1/kallisto_linux-v0.45.1.tar.gz

##### Extract and move to usr/local/bin

tar xvzf kallisto*gz
cp sudo cp kallisto_linux-v0.45.1/kallisto /usr/local/bin/

```

###  Run Indexing

`kallisto index -i human_pc.idx -k 25 reference/gencode.vM20.pc_transcripts.fa`

### Run Alignment
`kallisto quant -i human_pc.idx -o XX samples/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz samples/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --pseudobam`

---

## HISAT2


### Download Data

```bash
mkdir -p hisat/reference
cd hisat/reference/
wget https://transfer.sh/z8zqO/chr22_with_ERCC92.fa.gz
gunzip *gz
cd ..
```






