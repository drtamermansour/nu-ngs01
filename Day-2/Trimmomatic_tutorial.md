Trimmomatic Tutorial
====================

Install the software
```
source activate ngs1
conda install -c bioconda trimmomatic 
```

Run Trimmomatic
```
mkdir ~/workdir/trimmed && cd ~/workdir/trimmed 
f1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.fq.gz"
f2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.fq.gz"
newf1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.trim.fq.gz"
newf2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.trim.fq.gz"
newf1U="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.se.trim.fq.gz"
newf2U="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.se.trim.fq.gz"

adap="/home/ngs-01/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36
```
