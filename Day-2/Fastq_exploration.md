# FASTA/Q manipulations

## Download example data
```
mkdir -p ~/workdir/sample_data && cd ~/workdir/sample_data
wget https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/sample_data/SP1.fq
wget https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/sample_data/hairpin.fa.gz
```

## Install [SeqKit](https://bioinf.shenwei.me/seqkit/)
```
conda activate ngs1
conda install -c bioconda seqkit
```

## Quick overview
```
cd ~/workdir/
## What seqkit can do?
seqkit -h > seqkit.help  ## compare to the software homepage (https://bioinf.shenwei.me/seqkit/) 
## some stats
seqkit stat sample_data/*
## GC content: seqkit fx2tab coverts FASTA/Q to 3-column tabular format (1st: name/ID, 2nd: sequence, 3rd: quality), and can also provide various information in new columns, including sequence length, GC content/GC skew, alphabet
seqkit fx2tab --name --only-id --gc sample_data/hairpin.fa.gz
```

## More detailed tutorials

*  https://bioinf.shenwei.me/seqkit/usage/
*  https://bioinf.shenwei.me/seqkit/tutorial/
