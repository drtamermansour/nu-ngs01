# Functional Annotation
### Install
##### Trinotate automatic functional annotation of transcriptomes
##### TransDecoder for predicting coding regions in transcripts

``` 
source activate ngs1
conda install TransDecoder
conda install Trinotate
```
### Identification of likely protein-coding regions in transcripts

```
cd ~/workdir/trinity/trinity_out_dir
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta
ls -1 |grep transdecoder
```
### Sequence homology searches

```
mkdir -p ~/workdir/trinotate && cd ~/workdir/trinotate
blastx -db ../trinity/trinity_out_dir/data/mini_sprot.pep \
         -query Trinity.fasta -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > swissprot.blastx.outfmt6
          
blastp -query Trinity.fasta.transdecoder.pep \
         -db ../trinity/trinity_out_dir/data/mini_sprot.pep -num_threads 2 \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
          > swissprot.blastp.outfmt6

hmmpress ../data/Pfam-A.hmm

hmmscan --cpu 2 --domtblout TrinotatePFAM.out \
          ../trinity/trinity_out_dir/data/Pfam-A.hmm \
          Trinity.fasta.transdecoder.pep
          
```

### Preparing and Generating a Trinotate Annotation Report

```
Trinotate Trinotate.sqlite init \
     --gene_trans_map ../trinity/trinity_out_dir/Trinity.fasta.gene_trans_map \
     --transcript_fasta ../trinity/trinity_out_dir/Trinity.fasta \
     --transdecoder_pep Trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite \
       LOAD_swissprot_blastx swissprot.blastx.outfmt6
       
Trinotate Trinotate.sqlite \
       LOAD_swissprot_blastp swissprot.blastp.outfmt6
 
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
```
### Generate the Trinotate Annotation Report

```
$TRINOTATE_HOME/Trinotate Trinotate.sqlite report > Trinotate.xls
less Trinotate.xls
```
