#### Load and prep the metadata object
## Note1: This table should have the samples names, conditions to be compared, and batch info. 
## Note2: Columns that connect the quantification files to the samples rows
# for tximeta package, the metadata table has to have 2 columns: names (to carry the samples name) and files (for file paths)
# This is not a requirement for the tximport package which takes a vector of counting files paths while the DESeqDataSetFromTximport function of DESeq2 takes the metadata object to merge with tximport object. However,if the order in the files vector is not the same order in coldata table OR the names of the files vectors are matching rownames of the coldata, the final columns names will be wrong.  
# Therefore, the suggested design here should prevent accidental mistakes with either packages     

# 1. load the metadata file (select the column of samples names to be the row.names)
coldata <- read.csv("sample_table.csv", row.names=1, stringsAsFactors=FALSE)
# have a look on the metadata
coldata
# 2. change the condition and batch columns to factors
coldata$dex = factor(coldata$dex)
coldata$cell = factor(coldata$cell)
# 3. Add column for the file paths of the quantification files
dir=getwd()
coldata$files <- file.path(dir, "quants", rownames(coldata), "quant.sf") ## This way the order of the files will be the same order of rows 
# Check if these paths exist
file.exists(coldata$files)
# 4. Add extra column for names (to stasify the requirement of tximeta package.
coldata$names <- rownames(coldata)
# let us have another look on the metadata object
coldata
########################################################################################################

#### Import and summarize transcript-level estimates for transcript- and gene-level analysis using tximport or tximeta
#### Note: These packages just read and sum the raw and TPM counts. It does not do any adjustments  
## tximport package: 
# http://bioconductor.org/packages/release/bioc/html/tximport.html
# Tutorial: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# tximport Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. 
# tximport produces a "list" that can be imported into DESeq2 by DESeqDataSetFromTximport(). 
# 1. install and load the tximport package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
if (!requireNamespace("tximport", quietly = TRUE))
    BiocManager::install("tximport")
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("readr")
if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")
# 2. Generate tx2gene table: a two-column data.frame linking transcript id (column 1) to gene id (column 2). the column names are not relevant, but this column order must be used. 
# You should use the suitable annotation database according to the reference you used for previous step. For example, if your transcriptome is based on GRCh37, you should use "TxDb.Hsapiens.UCSC.hg19.knownGene" but if you use annotation of GRCh38, you should use "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb  ## This DB is based on GENCODE v32. If I would repeat this analysis, I will change my Gencode transcriptome to v32 instead of v36 
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2geneV1 <- na.omit(tx2gene)
# Alternatively, make or get your own table. For example, I downloaded Gencode gencode.v36.metadata.EntrezGene metadata file
tx2geneV2 <- read.table("gencode.v36.metadata.EntrezGene", col.names=c("TXNAME","GENEID"), sep = "\t")
# 3. import and summarize transcript-level estimates by gene. tximport knows the quantification format of many counting programs including salmon
# important note: if the order in the files vector is not the same order in coldata table OR the names of the files vectors are matching rownames of the coldata, the final columns names will be wrong.  
library("tximport")
library("readr")
txi <- tximport(coldata$files, type="salmon", tx2gene=tx2geneV2) ## txi is list of several objects: “counts” - the estimated fragment counts for each gene and sample, “abundance” - the estimated transcript abundances in TPM, and “length” - the effective gene lengths which include changes in length due to biases as well as due to transcript usage
# 4. import into DESeq2
# See more about "DESeqDataSet" later in the import step of tximeta package
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ dex)




## tximeta package: 
# http://bioconductor.org/packages/release/bioc/html/tximeta.html
# tximeta is another Bioconductor package, tximeta (Love et al. 2020), extends tximport, offering the same functionality, plus the additional benefit of automatic addition of annotation metadata for commonly used transcriptomes (GENCODE, Ensembl, RefSeq for human and mouse).
# tximeta produces a "SummarizedExperiment" with additional metadata that can be imported into DESeq2 using DESeqDataSet().
# 1. install and load the tximeta package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("tximeta", quietly = TRUE))
    BiocManager::install("tximeta")
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("readr")
# 2. No need to generate tx2gene bc it has annotation metadata for GENCODE transcriptome of Human
# 3. import and summarize transcript-level estimates by gene. tximeta knows the quantification format of many counting programs including salmon
library("tximeta")
library("readr")
# 3a) import
se <- tximeta(coldata) 
# explore the se object
dim(se)
head(rownames(se))
# 3b) summarize by gene 
gse <- summarizeToGene(se)
# explore the gse object
gse
# several important metadata about gse (e.g. mapping rates, no of filtered reads, counting parameters, etc) can be reterived by
metadata(gse)
# gse (as well as se) are SummarizedExperiment objects. It is composed of:
# 1) assays: list of matrices. “counts” - the estimated fragment counts for each gene and sample, “abundance” - the estimated transcript abundances in TPM, and “length” - the effective gene lengths which include changes in length due to biases as well as due to transcript usage. The first matrix in the list can be pulled out via "assay".
assayNames(gse)
head(assays(gse)[[1]], 3) # equal head(assay(gse), 3)  ## “counts” - the estimated fragment counts for each gene and sample
head(assays(gse)[[2]], 3)                              ## “abundance” - the estimated transcript abundances in TPM
head(assays(gse)[[3]], 3)                              ## “length” - the effective gene lengths which include changes in length due to biases as well as due to transcript usage
colSums(assays(gse)[[1]]) # equal colSums(assay(gse))
colSums(assays(gse)[[2]])
# all the matrices in assays have the same dim, rownames and colnames
dim(gse)
head(rownames(gse))
colnames(gse)
# 2) rowRanges: contains information about the genomic ranges. for our object is the GRanges of the genes (from the left-most position of all the transcripts to the right-most position of all the transcripts).
rowRanges(gse)
seqinfo(rowRanges(gse))
# 3) colData: contains information about the samples. 
colData(gse)

# 4. import into DESeq2
library("DESeq2")
## DESeq2 store the data throughout the analysis in a class named "DESeqDataSet". It is built on top of the SummarizedExperiment class as tximeta output
## One of the two main differences is that DESeqDataSet class enforces that the values in this matrix are non-negative integers.
## A second difference is that the DESeqDataSet has an associated design formula. This formula is specified at the beginning of the analysis, as it will inform many of the DESeq2 functions how to treat the samples in the analysis (one exception is the size factor estimation, i.e., the adjustment for differing library sizes, which does not depend on the design formula). The design formula tells which columns in the sample information table (colData) specify the experimental design and how these factors should be used in the analysis. The simplest design formula for differential expression would be "~ condition", while "~ batch + condition" will control for the effect of the batch column before calculating the DE based on the condition column
ddsTximeta <- DESeqDataSet(gse, design = ~ cell + dex)


## Note on factor levels
# check current levels of comparisons 
levels(ddsTximeta$dex)
# change the names of the levels. It is critical when one renames levels to not change the order
levels(ddsTximeta$dex) <- c("Dexamethasone", "untreated")
# it is prefered in R that the first level of a factor be the reference level (e.g. control, or untreated samples)
ddsTximeta$dex <- relevel(ddsTximeta$dex, ref = "untreated")
# If you subset the columns of a DESeqDataSet removing all samples of one or more levels, you can use droplevels function to remove levels with no samples
ddsTximeta$dex <- droplevels(ddsTximeta$dex)

## Note on Collapsing technical replicates 
## To combine the counts from technical replicates into single columns of the count matrix, use collapseReplicates function 
ddsTximeta$sample_rep = factor(c("S1","S2","S1","S2","S3","S4","S5","S6"))
ddsColl <- collapseReplicates(ddsTximeta, ddsTximeta$sample_rep, ddsTximeta$Run)
colData(ddsColl)
colnames(ddsColl)


#######################################################
## Exploratory analysis and visualization

## 1. Pre-filtering: This step is to reduce the memory size of the DESeq data object. This is different from the critical "independent filtering" required to increase power the analysis before multiple testing
# How many rows do we already have?
nrow(ddsTximeta)
# Filter based on the total number of reads in all samples 
keep <- rowSums(counts(ddsTximeta)) > 10
ddsTximeta <- ddsTximeta[keep,]
nrow(ddsddsTximeta)
# Filter based on having mimimum no of samples with minimal no of reads
keep <- rowSums(counts(ddsTximeta) >= 10) >= 3
ddsTximeta <- ddsTximeta[keep,]
nrow(ddsddsTximeta)

# 1. install and load of visualization package
if (!requireNamespace("vsn", quietly = TRUE))
    BiocManager::install("vsn")
if (!requireNamespace("hexbin", quietly = TRUE))
    BiocManager::install("hexbin")
if (!requireNamespace("dplyr", quietly = TRUE))
    BiocManager::install("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
    BiocManager::install("ggplot2")

# 2. Test for variance stabilizing
# Exploratory analysis of multidimensional data (e.g. clustering and PCA) works best for data that has the same range of variance at different ranges of the mean values (i.e. homoskedastic data). However, for RNA-seq counts, the expected variance grows with the mean. Therefor, the resulting plot typically depends mostly on the genes with highest counts because they show the largest absolute differences between samples
# Let us test the relation between mean and variance in real and simulated RNA-count data
library("hexbin")
library("vsn")

# 1. simulated data
lambda <- 10^seq(from = -1, to = 2, length = 1000) ## jpeg('lambda.jpg'); plot(density(lambda)); dev.off();     
cts_sim <- matrix(rpois(1000*100, lambda), ncol = 100)
jpeg('cts_sim.meanSdPlot.jpg')
meanSdPlot(cts_sim, ranks = FALSE)
dev.off()

# b. normalized real data
# Calc normalization factors
ddsTximeta <- estimateSizeFactors(ddsTximeta) 
# Check for the new matrix of normalization factors 
assayNames(ddsTximeta)
# Obtain a new normalized count matrix
cts_real <- counts(ddsTximeta, normalized=TRUE)
jpeg('cts_real.meanSdPlot.jpg')
meanSdPlot(cts_real, ranks = FALSE)
dev.off()
# We can improve the visualization if we excluded the outliars
cts_real2 <- cts_real[(rowMeans(cts_real) < 50000),]
jpeg('cts_real2.meanSdPlot.jpg')
meanSdPlot(cts_real2, ranks = FALSE)
dev.off()

# A simple and often used strategy to avoid this is to take the logarithm of the normalized count values plus a pseudocount of 1; however, depending on the choice of pseudocount, now the genes with the very lowest counts will contribute a great deal of noise to the resulting plot, because taking the logarithm of small counts actually inflates their variance.
log.cts_sim.one <- log2(cts_sim + 1)
jpeg('log_cts_sim_one.meanSdPlot.jpg')
meanSdPlot(log.cts_sim.one, ranks = FALSE)
dev.off()

log.cts_real.one <- log2(cts_real+1)
jpeg('log_cts_real_one.meanSdPlot.jpg')
meanSdPlot(log.cts_real.one, ranks = FALSE)
dev.off()

# As a solution, DESeq2 offers two transformations for counts to stabilize the variance across the means. They give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value. 
# Variance stabilizing transformation (vst function) (Anders and Huber 2010): faster & less sensitive to high count outliers. Better large datasets (n > 30)
# Regularized-logarithm transformation (rlog function) (Love, Huber, and Anders 2014): work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples
# You can perform both transformations and compare the meanSdPlot or PCA plots generated
# Both vst and rlog return a DESeqTransform object which is based on the SummarizedExperiment class. The transformed values are no longer counts, and are stored in the assay slot.
# Both functions have option called "blind". If TRUE (default), transformation is fully unsupervised. If FALSE, the differences the impact of design's variables (cell line and treatment) will be considered so that it does not affect the expected variance-mean trend of the experiment
vsd <- vst(ddsTximeta, blind = FALSE)
rld <- rlog(ddsTximeta, blind = FALSE)

# let see the difference 
jpeg('vsd_cts_real.meanSdPlot.jpg')
meanSdPlot(assay(vsd), ranks = FALSE)
dev.off()

jpeg('rld_cts_real.meanSdPlot.jpg')
meanSdPlot(assay(rld), ranks = FALSE)
dev.off()

# Remember, we are doing stabilization of the variance across the means to avoid biased impact of high or low genes on sample-sample distance
# So let us see how this would affect such distance between the 1st 2 samples
library("dplyr")
library("ggplot2")

df <- bind_rows(
  as_tibble(log2(cts_real[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

head(df)

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

jpeg('stable_cts_real.jpg')
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)  
dev.off()






