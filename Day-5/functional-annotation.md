# Functional Annotation

### Resources:

*   [Backgroud](https://github.com/Trinotate/Trinotate.github.io/wiki)
*   [Installation and annotation](https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required)
*   [Loading results into a Trinotate SQLite Database](https://github.com/Trinotate/Trinotate.github.io/wiki/Loading-generated-results-into-a-Trinotate-SQLite-Database-and-Looking-the-Output-Annotation-Report)

### Install 
``` 
conda activate ngs1
conda install -c bioconda transdecoder  # TransDecoder for predicting coding regions in transcripts
conda install -c bioconda trinotate     # Trinotate automatic functional annotation of transcriptomes
conda install -c bioconda hmmer         # for searching sequence databases for sequence homologs
```
### Download annotation databases 
```
## The analysis requires download of several annotation databases 
## Trinotate comes with a script to download recent version of all required data resources including the latest version of swissprot, pfam, and other companion resources, create and populate a Trinotate boilerplate sqlite database (Trinotate.sqlite) 
## WARINING: DO NOT RUN THIS command. It will take too much time and space ## Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate  

## For this tutorial, we will skip this step and install smaller versions of these resources 
## Download SwissProt database (we are using small version for demonstration)
mkdir -p ~/workdir/databases && cd ~/workdir/databases
wget https://data.cyverse.org/dav-anon/iplant/home/drtamermansour/Trinotate_sample_DB/uniprot_sprot.sample.pep.gz
gunzip uniprot_sprot.sample.pep.gz

## Download Pfam database (we are using small version for demonstration)
wget https://data.cyverse.org/dav-anon/iplant/home/drtamermansour/Trinotate_sample_DB/Pfam-A.sample.hmm.gz
gunzip Pfam-A.sample.hmm.gz

## Download and prepare Trinotate boilerplate sqlite database 
wget https://data.cyverse.org/dav-anon/iplant/home/drtamermansour/Trinotate_sample_DB/Trinotate.sample.UniprotIndex.gz
gunzip Trinotate.sample.UniprotIndex.gz
prefix="Trinotate"
EMBL_dat_to_Trinotate_sqlite_resourceDB.pl --sqlite "$prefix".sqlite --create
EMBL_dat_to_Trinotate_sqlite_resourceDB.pl --sqlite "$prefix".sqlite --uniprot_index $prefix.sample.UniprotIndex
```

### Coding Region Identification using [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki)

Check how does TransDecoder identify likely coding sequences on the its [wiki page](https://github.com/TransDecoder/TransDecoder/wiki) 
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

## Create the blast DB
makeblastdb -in ~/workdir/databases/uniprot_sprot.sample.pep -dbtype prot

## Run the blast aanalysis
blastx -db ~/workdir/databases/uniprot_sprot.sample.pep \
         -query ~/workdir/trinity/trinity_out_dir/Trinity.fasta  \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 1 \
          > swissprot.blastx.outfmt6
          
blastp -db ~/workdir/databases/uniprot_sprot.sample.pep \
         -query ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep \
         -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 1 \
          > swissprot.blastp.outfmt6

## Assessment
analyze_blastPlus_topHit_coverage.pl swissprot.blastx.outfmt6 \
                                     ~/workdir/trinity/trinity_out_dir/Trinity.fasta \
                                     ~/workdir/databases/uniprot_sprot.sample.pep | column -t

analyze_blastPlus_topHit_coverage.pl swissprot.blastp.outfmt6 \
                                     ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep \
                                     ~/workdir/databases/uniprot_sprot.sample.pep | column -t
```

Here is the output for swissprot.blastx.outfmt6
```
#hit_pct_cov_bin  count_in_bin  >bin_below
100               6             6
90                5             11
80                3             14
70                4             18
60                10            28
50                8             36
40                13            49
30                14            63
20                18            81
10                14            95
```

The above table lists bins of percent length coverage of the best matching protein sequence along with counts of proteins found within that bin. For example, 6 proteins are matched by 90-100% of their length. There are 5 matched by 80-90% of their length. The third column provides a running total, indicating that 11 transcripts match more than 80% of their length, and 14 transcripts match more than 70% of their length, etc.

## b) HMMER/PFAM Protein Domain Identification (http://hmmer.org/)

Instead of searching aganist all known sequences, HMMER compares sequences against all known domain families using a prebuilt dataset of statistical models e.g. Pfam profiles
```
## Run the HMMER aanalysis
hmmpress ~/workdir/databases/Pfam-A.sample.hmm
hmmscan --cpu 1 --domtblout TrinotatePFAM.out \
          ~/workdir/databases/Pfam-A.sample.hmm \
          ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep      
```


### Preparing and Generating a Trinotate Annotation Report
```
# Generate a map for genes and their corresponding transcripts 
get_Trinity_gene_to_trans_map=$(find /home/ngs/miniconda3/envs/ngs1 -name get_Trinity_gene_to_trans_map.pl)
$get_Trinity_gene_to_trans_map ~/workdir/trinity/trinity_out_dir/Trinity.fasta >  Trinity.fasta.gene_trans_map

# Load these info into the Trinotate sqlite database
Trinotate ~/workdir/databases/Trinotate.sqlite init \
     --gene_trans_map Trinity.fasta.gene_trans_map \
     --transcript_fasta ~/workdir/trinity/trinity_out_dir/Trinity.fasta \
     --transdecoder_pep ~/workdir/trinity/trinity_out_dir/Trinity.fasta.transdecoder.pep

# Loading BLAST homologies
Trinotate ~/workdir/databases/Trinotate.sqlite \
       LOAD_swissprot_blastx swissprot.blastx.outfmt6
       
Trinotate ~/workdir/databases/Trinotate.sqlite \
       LOAD_swissprot_blastp swissprot.blastp.outfmt6

# Load Pfam domain entries
Trinotate ~/workdir/databases/Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
```

### Generate the Trinotate Annotation Report
```
Trinotate ~/workdir/databases/Trinotate.sqlite report > Trinotate.xls
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
