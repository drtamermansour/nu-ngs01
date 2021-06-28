## Building resources on HPC

### Pfam
```
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
subset=$(cat Pfam-A.hmm | grep -n "//" | head -n400 | tail -n1 | cut -d":" -f1)
head -n "$subset" Pfam-A.hmm > Pfam-A.sample.hmm
gzip Pfam-A.sample.hmm
```

### uniprot_sprot
```
wget http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
conda activate ngs1
conda install -c bioconda trinotate
prefix="Trinotate"
EMBL_swissprot_parser.pl uniprot_sprot.dat.gz $prefix
mv uniprot_sprot.dat.gz.pep uniprot_sprot.pep
grep _HUMAN $prefix.UniprotIndex > $prefix.UniprotIndex.human
cat $prefix.UniprotIndex.human | cut -f1 | uniq | paste - - - - - - - - - - | cut -f1 | head -n1000 > human_sample
cat human_sample | grep -A1 --no-group-separator -Ff - uniprot_sprot.pep > uniprot_sprot.sample.pep
cat human_sample | grep -A1 --no-group-separator -Ff - $prefix.UniprotIndex > $prefix.sample.UniprotIndex
gzip uniprot_sprot.sample.pep
gzip $prefix.sample.UniprotIndex
```
