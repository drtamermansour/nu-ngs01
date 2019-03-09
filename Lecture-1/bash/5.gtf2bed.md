# Conversion of GTF to Bed file 

## Install the tools

```
conda install -c bioconda ucsc-gtftogenepred
mkdir ~/workdir/scripts && cd ~/workdir/scripts
wget https://raw.githubusercontent.com/drtamermansour/horse_trans/master/scripts/genePredToBed
chmod 755 genePredToBed
```

## Run the conversion code

```
mkdir ~/workdir/gtfToBed && cd ~/workdir/gtfToBed
gtfToGenePred gencode.v29.annotation.gtf gencode.v29.annotation.gpred
cat gencode.v29.annotation.gpred | ~/workdir/scripts/genePredToBed > gencode.v29.annotation.bed
```

## Explore the files


