# Functional Annotation

### Install
``` 
conda activate ngs1
conda install -c bioconda transdecoder  # TransDecoder for predicting coding regions in transcripts
conda install -c bioconda trinotate     # Trinotate automatic functional annotation of transcriptomes
conda install -c bioconda hmmer         # for searching sequence databases for sequence homologs
```

### Coding Region Identification using [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki)
```
cd ~/workdir/trinity/trinity_out_dir
# Extract the long open reading frames
TransDecoder.LongOrfs -t Trinity.fasta &> longOrf.log
# Predict the more likely coding regions
TransDecoder.Predict -t Trinity.fasta &> predict.log
```

Output files 
```
## Trinity.fasta.transdecoder_dir
longest_orfs.pep   : all ORFs meeting the minimum length criteria, regardless of coding potential.
longest_orfs.gff3  : positions of all ORFs as found in the target transcripts
longest_orfs.cds   : the nucleotide coding sequence for all detected ORFs
longest_orfs.cds.top_500_longest   : the top 500 longest ORFs, used for training a Markov model for coding sequences.
hexamer.scores                     : log likelihood score for each k-mer  (coding/random)
longest_orfs.cds.scores            : the log likelihood sum scores for each ORF across each of the 6 reading frames
longest_orfs.cds.scores.selected   : the accessions of the ORFs that were selected based on the scoring criteria (described at top)

## in your current working directory
transcripts.fasta.transdecoder.pep : peptide sequences for the final candidate ORFs; all shorter candidates within longer ORFs were removed.
transcripts.fasta.transdecoder.cds  : nucleotide sequences for coding regions of the final candidate ORFs
transcripts.fasta.transdecoder.gff3 : positions within the target transcripts of the final selected ORFs
transcripts.fasta.transdecoder.bed  : bed-formatted file describing ORF positions, best for viewing using GenomeView or IGV.
```

### Sequence homology searches
## a) NCBI BLAST+ aganist SwissProt database (The UniProt Knowledgebase which include the Manually annotated proteins)
```
mkdir -p ~/workdir/trinotate && cd ~/workdir/trinotate

## Download SwissProt database (we are using small version for demonstration)
mkdir -p ~/workdir/databases && cd ~/workdir/databases
wget https://github.com/trinityrnaseq/KrumlovTrinityWorkshopJan2018/raw/master/data/mini_sprot.pep.gz
gunzip mini_sprot.pep.gz

## create the blast DB
makeblastdb -in mini_sprot.pep -dbtype prot

## Run the blast aanalysis
blastx -db ~/workdir/databases/mini_sprot.pep \
         -query ~/workdir/trinity/trinity_out_dir/Trinity.fasta  \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 1 \
          > swissprot.blastx.outfmt6
          
blastp -db ~/workdir/databases/mini_sprot.pep \
         -query ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 1 \
          > swissprot.blastp.outfmt6

## Assessment
analyze_blastPlus_topHit_coverage.pl swissprot.blastx.outfmt6 \
                                     ~/workdir/trinity/trinity_out_dir/Trinity.fasta \
                                     ~/workdir/databases/mini_sprot.pep | column -t

analyze_blastPlus_topHit_coverage.pl swissprot.blastp.outfmt6 \
                                     ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep \
                                     ~/workdir/databases/mini_sprot.pep | column -t
```

Here is the output for swissprot.blastx.outfmt6
```
#hit_pct_cov_bin  count_in_bin  >bin_below
100               4             4
90                3             7
80                1             8
70                2             10
60                2             12
50                2             14
40                3             17
30                1             18
20                4             22
10                3             25
```

## b) HMMER/PFAM Protein Domain Identification (http://hmmer.org/)
```
hmmpress ~/workdir/databases/Pfam-A.hmm
hmmscan --cpu 1 --domtblout TrinotatePFAM.out \
          ~/workdir/databases/Pfam-A.hmm \
          ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep      
```

### Trinotate boilerplate sqlite database 
```
# Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate  ## This step will download several data resources including the latest version of swissprot, pfam, and other companion resources, create and populate a Trinotate boilerplate sqlite database (Trinotate.sqlite)
# Today we will download the file ready for use
wget https://abdelrahmanma.com/NGS01/Trinotate.sqlite.tar.xz
tar -xvf Trinotate.sqlite.tar.xz
```

### Preparing and Generating a Trinotate Annotation Report
```
Trinotate Trinotate.sqlite init \
     --gene_trans_map ~/workdir/trinity/trinity_out_dir/Trinity.fasta.gene_trans_map \
     --transcript_fasta ~/workdir/trinity/trinity_out_dir/Trinity.fasta \
     --transdecoder_pep ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep

Trinotate Trinotate.sqlite \
       LOAD_swissprot_blastx swissprot.blastx.outfmt6
       
Trinotate Trinotate.sqlite \
       LOAD_swissprot_blastp swissprot.blastp.outfmt6
 
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
```

### Generate the Trinotate Annotation Report
```
Trinotate Trinotate.sqlite report > Trinotate.xls
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
