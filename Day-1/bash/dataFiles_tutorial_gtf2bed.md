# Conversion of GTF to Bed file 

## Install [bioconda](https://docs.conda.io/en/latest/miniconda.html) if you don't have it already

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
## restart the terminal
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -y --name ngs1 python=3.6
```

## Install the tools

```
source activate ngs1
conda install -c bioconda ucsc-gtftogenepred
mkdir ~/workdir/scripts && cd ~/workdir/scripts
wget https://raw.githubusercontent.com/drtamermansour/horse_trans/master/scripts/genePredToBed
chmod 755 genePredToBed
```

## Download, unzip and explore the GTF file
```
mkdir ~/workdir/gtfToBed && cd ~/workdir/gtfToBed
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz
```

## Run the conversion code

```
gtfToGenePred gencode.v29.annotation.gtf gencode.v29.annotation.gpred
cat gencode.v29.annotation.gpred | ~/workdir/scripts/genePredToBed > gencode.v29.annotation.bed
```

## Explore the files


