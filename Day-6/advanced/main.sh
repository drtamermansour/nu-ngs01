## RNA-seq workflow: gene-level exploratory analysis and differential expression
## Source1: http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html (10/16/2019)
## Source2: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html (10/27/2020)
###########################################
## Software setup
conda create -n geneDE
conda activate geneDE
conda install salmon
conda install r

##########################################
## Experimental data
## Human Airway Smooth Muscle Transcriptome Changes in 4 different cell lines in Response to Asthma Medications(Dexamethasone, Albuterol, or both)
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
## You can use different ways to download this dataset like SRA-toolkit but I am using SRA explorer(https://sra-explorer.info/) to get download links by searching using the project ID "SRP033351". Select the 16 samples that show up and add them to collection. Now go to the cart and select the bash script for downloading the Fastq data. Note: All samples have _1.fastq.gz and _2.fastq.gz files for PE data but some samples have a 3rd file for unpaired reads.
## NOTE: SRR1039513,SRR1039514,SRR1039515,SRR1039520,SRR1039521,SRR1039522,SRR1039523 have ercc_mix (based on the SRA run selector metadata)
## This tutorial is using PE-reads from 8 samples only [SRR1039508,SRR1039509,SRR1039512,SRR1039513,SRR1039516,SRR1039517,SRR1039520,SRR1039521] to compare ubtreated vs Dexamethasone treatment in the 4 cell lines 
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513.fastq.gz -o SRR1039513_GSM1275867_N052611_Dex_Homo_sapiens_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513_1.fastq.gz -o SRR1039513_GSM1275867_N052611_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513_2.fastq.gz -o SRR1039513_GSM1275867_N052611_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512_1.fastq.gz -o SRR1039512_GSM1275866_N052611_untreated_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512_2.fastq.gz -o SRR1039512_GSM1275866_N052611_untreated_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509_1.fastq.gz -o SRR1039509_GSM1275863_N61311_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509_2.fastq.gz -o SRR1039509_GSM1275863_N61311_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039510/SRR1039510_1.fastq.gz -o SRR1039510_GSM1275864_N61311_Alb_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039510/SRR1039510_2.fastq.gz -o SRR1039510_GSM1275864_N61311_Alb_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz -o SRR1039508_GSM1275862_N61311_untreated_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz -o SRR1039508_GSM1275862_N61311_untreated_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/005/SRR1039515/SRR1039515.fastq.gz -o SRR1039515_GSM1275869_N052611_Alb_Dex_Homo_sapiens_RNA-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/005/SRR1039515/SRR1039515_1.fastq.gz -o SRR1039515_GSM1275869_N052611_Alb_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/005/SRR1039515/SRR1039515_2.fastq.gz -o SRR1039515_GSM1275869_N052611_Alb_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516.fastq.gz -o SRR1039516_GSM1275870_N080611_untreated_Homo_sapiens_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516_1.fastq.gz -o SRR1039516_GSM1275870_N080611_untreated_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516_2.fastq.gz -o SRR1039516_GSM1275870_N080611_untreated_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039519/SRR1039519.fastq.gz -o SRR1039519_GSM1275873_N080611_Alb_Dex_Homo_sapiens_RNA-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039519/SRR1039519_1.fastq.gz -o SRR1039519_GSM1275873_N080611_Alb_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039519/SRR1039519_2.fastq.gz -o SRR1039519_GSM1275873_N080611_Alb_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039511/SRR1039511_1.fastq.gz -o SRR1039511_GSM1275865_N61311_Alb_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039511/SRR1039511_2.fastq.gz -o SRR1039511_GSM1275865_N61311_Alb_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039518/SRR1039518_1.fastq.gz -o SRR1039518_GSM1275872_N080611_Alb_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039518/SRR1039518_2.fastq.gz -o SRR1039518_GSM1275872_N080611_Alb_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/007/SRR1039517/SRR1039517_1.fastq.gz -o SRR1039517_GSM1275871_N080611_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/007/SRR1039517/SRR1039517_2.fastq.gz -o SRR1039517_GSM1275871_N080611_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/004/SRR1039514/SRR1039514_1.fastq.gz -o SRR1039514_GSM1275868_N052611_Alb_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/004/SRR1039514/SRR1039514_2.fastq.gz -o SRR1039514_GSM1275868_N052611_Alb_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520.fastq.gz -o SRR1039520_GSM1275874_N061011_untreated_Homo_sapiens_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520_1.fastq.gz -o SRR1039520_GSM1275874_N061011_untreated_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520_2.fastq.gz -o SRR1039520_GSM1275874_N061011_untreated_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039522/SRR1039522.fastq.gz -o SRR1039522_GSM1275876_N061011_Alb_Homo_sapiens_RNA-Seq.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039522/SRR1039522_1.fastq.gz -o SRR1039522_GSM1275876_N061011_Alb_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039522/SRR1039522_2.fastq.gz -o SRR1039522_GSM1275876_N061011_Alb_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039521/SRR1039521.fastq.gz -o SRR1039521_GSM1275875_N061011_Dex_Homo_sapiens_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039521/SRR1039521_1.fastq.gz -o SRR1039521_GSM1275875_N061011_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039521/SRR1039521_2.fastq.gz -o SRR1039521_GSM1275875_N061011_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039523/SRR1039523_1.fastq.gz -o SRR1039523_GSM1275877_N061011_Alb_Dex_Homo_sapiens_RNA-Seq_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039523/SRR1039523_2.fastq.gz -o SRR1039523_GSM1275877_N061011_Alb_Dex_Homo_sapiens_RNA-Seq_2.fastq.gz 

########################
## Introduction:
## We start by abundance estimation for raw counts in each sample (We will try a couple ways).followed by 
## importing counts with an experiment design into one of a variety of Bioconductor packages for exploration and differential expression of the count data
## We will use DESeq2 for DE. DESeq has several import functions to work easily with different quantification and annotation programs (Check source1; DESeq2 import functions)
## Therefore, we will see how DESeq2 will import counts from different tools of abundance estimation 

########################
## Abundance estimation 
## You can align reads to the genome and then count the number of reads that are consistent with gene models (see the old NGS1 tutorial && the section of "Count matrix input" from Source2)
## Alternatively, transcript-level quantification data (e.g. by  Salmon, kallisto or RSEM) followed by aggregation to the gene-level (by  tximport or tximeta) is more recommended. Advantages: (1)this approach corrects for any potential changes in gene length across samples (e.g. from differential isoform usage) (Trapnell et al. 2013); (2) some of these methods are substantially faster and require less memory and disk usage compared to alignment-based methods; and (3) it is possible to avoid discarding those fragments that can align to multiple genes with homologous sequence (Robert and Watson 2015).  

## Abundance estimation using Salmon (https://salmon.readthedocs.io/en/latest/index.html)
# 1. Download the reference transcriptome (I will get the most recent version of Gencode transcriptome)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz

# 2. create Salmon index
salmon index --gencode -t gencode.v36.transcripts.fa.gz -i gencode.v36_salmonIndex ## the --gencode option is only to strip the extra information from the transcript names

# 3. run the alignment
for sample in {SRR1039508,SRR1039509,SRR1039512,SRR1039513,SRR1039516,SRR1039517,SRR1039520,SRR1039521};do
  echo $sample;
  R1=${sample}_*_1.fastq.gz
  salmon quant -i gencode.v36_salmonIndex --libType A -p 4 \
         --gcBias --numGibbsSamples 20 -o quants/$sample \
         -1 ${sample}_*_1.fastq.gz -2 ${sample}_*_2.fastq.gz
  ## --libType A: allow Salmon to automatically infer the lib type  (https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype)
  ## --numGibbsSamples and --numBootstraps options are mutually exclusive. Either one generate a background distribution to estimate the variance in abundance estimates
done

# 4. Import and summarize transcript-level estimates for transcript- and gene-level analysis
#    DESeq has several import functions to work easily with different quantification and annotation programs (Check source1; DESeq2 import functions)

# 4A. Expermint annotation
#     Typically, we have a table with detailed information for each of our samples.
#     For your own project, you might create such a comma-separated value (CSV) file using a text editor or spreadsheet software such as Excel.
#     In this tutorial I copied this file "sample_table.csv" from the airway package (an R package that summarizes RNA-seq & metadata of this study). 
#     You can see how to download and explore this package in explore_airway.R

# 4B. Transcript-to-gene mapping info: Package like "tximeta" can automatically do this for commonly used transcriptomes (GENCODE, Ensembl, RefSeq for human and mouse). 
#     However, other packages like "tximport" and other transcriptomes will reuire you to generate this file. 
#     tximport requires a two-column CSV linking transcript id (column 1) to gene id (column 2). the column names are not relevant, but this column order must be used
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.metadata.EntrezGene.gz
gunzip gencode.v36.metadata.EntrezGene.gz

# 4C. run the Import and summarization using tximport or tximeta R packages as described in DE.R
######################
