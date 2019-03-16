Sequence Alignment
==================


## Download reference file
```
cd ~/workdir/sample_data
wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz

OR

wget https://transfer.sh/7g7zo/dog.fa.gz -O dog_chr5.fa.gz

gunzip dog_chr5.fa.gz
```

BLAST Alignment
===============

## install [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs)
```
source activate ngs1
conda install -c bioconda blast 
```

## index your genome
```
mkdir -p ~/workdir/blast_align/blastIndex && cd ~/workdir/blast_align/blastIndex
ln -s ~/workdir/sample_data/dog_chr5.fa .
makeblastdb -in dog_chr5.fa -dbtype nucl
```

## sequence alignment
```
cd ~/workdir/blast_align
R1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.fq.gz"
R2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.fq.gz"
zcat $R1 $R2 | paste - - - - | awk 'BEGIN{OFS="\n";}{print ">"$1,$2}' > sample.fa
/usr/bin/time -v blastn -db blastIndex/dog_chr5.fa -query sample.fa -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send" -out BD143_TGACCA_L005.blastOut

## try the xml format
head -n500 sample.fa > sample_subset.fa
/usr/bin/time -v blastn -db blastIndex/dog_chr5.fa -query sample_subset.fa -outfmt 5 -out BD143_TGACCA_L005.xml
```

## convert to SAM
```
conda install -c conda-forge biopython 
wget https://gist.githubusercontent.com/ozagordi/099bdb796507da8d9426/raw/6ca66616fd545fbb15d94b079e46a7c55edb54c0/blast2sam.py
python blast2sam.py BD143_TGACCA_L005.xml > BD143_TGACCA_L005.sam

```


BWA Alignment
=============

## install [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)
```
source activate ngs1
conda install -c bioconda bwa 
```

## index your genome
```
mkdir -p ~/workdir/bwa_align/bwaIndex && cd ~/workdir/bwa_align/bwaIndex
ln -s ~/workdir/sample_data/dog_chr5.fa .
bwa index -a bwtsw dog_chr5.fa
```

## sequence alignment
```
cd ~/workdir/bwa_align
R1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.fq.gz"
R2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.fq.gz"
/usr/bin/time -v bwa mem bwaIndex/dog_chr5.fa $R1 $R2 > BD143_TGACCA_L005.sam
```
