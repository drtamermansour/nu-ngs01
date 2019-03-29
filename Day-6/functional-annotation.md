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
### Download  databases
```
mkdir -p ~/workdir/trinity/data && cd ~/workdir/trinity/data
wget https://get.station307.com/KP5ncqK7CUb/SWISSPROT-Hmm.tar.xz
tar -xvf SWISSPROT-Hmm.tar.xz
~/workdir/trinotate && cd ~/workdir/trinotate
wget https://get.station307.com/jYiGiAN78Qd/Trinotate.sqlite.tar.xz
tar -xvf Trinotate.sqlite.tar.xz
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

### Identifiers explanation

```
#gene_id        transcript_id   sprot_Top_BLASTX_hit    TrEMBL_Top_BLASTX_hit   RNAMMER prot_id prot_coords     sprot_Top_BLASTP_hit    TrEMBL_Top_BLASTP_hit   Pfam    SignalP TmHMM   eggnog  gene_ontology_blast     gene_ontology_pfam      transcript      peptide
TRINITY_DN144_c0_g1     TRINITY_DN144_c0_g1_i1  PUT4_YEAST^PUT4_YEAST^Q:1-198,H:425-490^74.24%ID^E:4e-29^RecName: Full=Proline-specific permease;^Eukaryota; Fungi; Dikarya; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces        .       .       .       .
   .       .       .       .       .       COG0833^permease        GO:0016021^cellular_component^integral component of membrane`GO:0005886^cellular_component^plasma membrane`GO:0015193^molecular_function^L-proline transmembrane transporter activity`GO:0015175^molecular_function^neutral amino acid transmembrane transporter activity`GO:0015812^biological_process^gamma-aminobutyric acid transport`GO:0015804^biological_process^neutral amino acid transport`GO:0035524^biological_process^proline transmembrane transport`GO:0015824^biological_process^proline transport      .       .       .
TRINITY_DN179_c0_g1     TRINITY_DN179_c0_g1_i1  ASNS1_YEAST^ASNS1_YEAST^Q:1-168,H:158-213^82.14%ID^E:5e-30^RecName: Full=Asparagine synthetase [glutamine-hydrolyzing] 1;^Eukaryota; Fungi; Dikarya; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces        .
   .       .       .       .       .       .       .       .       COG0367^asparagine synthetase   GO:0004066^molecular_function^asparagine synthase (glutamine-hydrolyzing) activity`GO:0005524^molecular_function^ATP binding`GO:0006529^biological_process^asparagine biosynthetic process`GO:0006541^biological_process^glutamine metabolic process`GO:0070981^biological_process^L-asparagine biosynthetic process    .       .       .
TRINITY_DN159_c0_g1     TRINITY_DN159_c0_g1_i1  ENO2_CANGA^ENO2_CANGA^Q:2-523,H:128-301^100%ID^E:4e-126^RecName: Full=Enolase 2;^Eukaryota; Fungi; Dikarya; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; Saccharomycetaceae; Nakaseomyces; Nakaseomyces/Candida clade      .       .       TRINITY_DN159_c0_g1_i1|m.1      2-523[+]        ENO2_CANGA^ENO2_CANGA^Q:1-174,H:128-301^100%ID^E:3e-126^RecName: Full=Enolase 2;^Eukaryota; Fungi; Dikarya; Ascomycota; Saccharomycotina; Saccharomycetes; Saccharomycetales; Saccharomycetaceae; Nakaseomyces; Nakaseomyces/Candida clade      .       PF00113.17^Enolase_C^Enolase, C-terminal TIM barrel domain^18-174^E:9.2e-79     .       .       .       GO:0005829^cellular_component^cytosol`GO:0000015^cellular_component^phosphopyruvate hydratase complex`GO:0000287^molecular_function^magnesium ion binding`GO:0004634^molecular_function^phosphopyruvate hydratase activity`GO:0006096^biological_process^glycolytic process     GO:0000287^molecular_function^magnesium ion binding`GO:0004634^molecular_function^phosphopyruvate hydratase activity`GO:0006096^biological_process^glycolytic process`GO:0000015^cellular_component^phosphopyruvate hydratase complex   .       .
...
```
