#### Load and prep the metadata object
#     Typically, we have a table with detailed information for each of our samples.
#     For your own project, you might create such a comma-separated value (CSV) file using a text editor or spreadsheet software such as Excel.
#     In this tutorial I copied this file "sample_table.csv" from the airway package (an R package that summarizes RNA-seq & metadata of this study). 

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
#     Package like "tximeta" can automatically do this for commonly used transcriptomes (GENCODE, Ensembl, RefSeq for human and mouse). 
#     However, other packages like "tximport" and other transcriptomes will reuire you to generate this file. 
#     tximport requires a two-column CSV linking transcript id (column 1) to gene id (column 2). the column names are not relevant, but this column order must be used

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
# TxDb.Hsapiens.UCSC.hg38.knownGene exposes an annotation databases generated from UCSC as TxDb objects.
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb  ## This DB is based on GENCODE v32. If I would repeat this analysis, I will change my Gencode transcriptome to v32 instead of v36 
# explore the avlaible key types
columns(txdb)
# generate the two-column table
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2geneV1 <- na.omit(tx2gene)
# Alternatively, make or get your own table. For example, I downloaded Gencode gencode.v36.metadata.EntrezGene metadata file
download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.metadata.EntrezGene.gz","gencode.v36.metadata.EntrezGene.gz")
library(R.utils)
gunzip('gencode.v36.metadata.EntrezGene.gz')
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
# 3. import and summarize transcript-level estimates by gene. tximeta knows the quantification format of many counting programs including salmon which is the default 
library("tximeta")
library("readr")
# 3a) import using tximeta functions. Arguments include:
#  type = c("none", "salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie"),
#  useHub: default is TRUE to first attempt to download a TxDb/EnsDb object from AnnotationHub. If FALSE, it will force downloading the GTF from FTP
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
## A second difference is that the DESeqDataSet has an associated design formula specified at the beginning of the analysis. 
#  The design formula informs many of the DESeq2 functions how to treat the samples in the analysis (one exception is the size factor estimation, i.e., the adjustment for differing library sizes, which does not depend on the design formula). 
#  The design formula tells which columns in the sample information table (colData) specify the experimental design and how these factors should be used in the analysis. 
#  The simplest design formula for differential expression would be "~ condition", while "~ batch + condition" will control for the effect of the batch column before calculating the DE based on the condition column
ddsTximeta <- DESeqDataSet(gse, design = ~ cell + dex)


## Alternatively, the function DESeqDataSetFromMatrix can be used if you already have a matrix of read counts prepared from another source.
## e.g. count matrices from alignment files by the featureCounts function in the Rsubread package. 
## The function needs a matrix of counts as integers, metadata file, and the design formula.
raw_counts <- assay(gse)
raw_counts <- apply(raw_counts, 2, function(x) as.integer(round(x)) )
ddsCountMat <- DESeqDataSetFromMatrix(countData = raw_counts, colData = coldata, design = ~ cell + dex)


## conclusions on the import step:
#  If you are using a transcript based alignment, tximport/tximeta are better to adjust automatically for transcript length while summerizing expression per gene
#     tximeta is better for human and mouse as it automate the annotation with the appropriate annotation version 
#     However for working with self-made annoatation or other species, tximport is the option
#  For the gene-level alignment approaches, importing from a matrix will be the better option


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
nrow(ddsTximeta)
# AND/OR Filter based on having mimimum no of samples with minimal no of reads
keep <- rowSums(counts(ddsTximeta) >= 10) >= 3
ddsTximeta <- ddsTximeta[keep,]
nrow(ddsTximeta)

# My recommendations: initial minimal pre-filtering is ok. Later automatic independent filtering will select for appropriate additional filteration

# 2. Test for variance stabilizing
# Exploratory analysis of multidimensional data (e.g. clustering and PCA) works best for data that has the same range of variance at different ranges of the mean values (i.e. homoskedastic data). However, for RNA-seq counts, the expected variance grows with the mean. Therefor, the resulting plot typically depends mostly on the genes with highest counts because they show the largest absolute differences between samples
# Let us test the relation between mean and variance in real and simulated RNA-count data

if (!requireNamespace("vsn", quietly = TRUE))
    BiocManager::install("vsn")
if (!requireNamespace("hexbin", quietly = TRUE))
    BiocManager::install("hexbin")

library("hexbin")
library("vsn")

# a. simulated data
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

# A simple and often used strategy to avoid variance instability is to take the logarithm of the normalized count values plus a pseudocount of 1; however, depending on the choice of pseudocount, now the genes with the very lowest counts will contribute a great deal of noise to the resulting plot, because taking the logarithm of small counts actually inflates their variance.
log.cts_sim.one <- log2(cts_sim + 1)
jpeg('log_cts_sim_one.meanSdPlot.jpg')
meanSdPlot(log.cts_sim.one, ranks = FALSE)
dev.off()

log.cts_real.one <- log2(cts_real + 1)
jpeg('log_cts_real_one.meanSdPlot.jpg')
meanSdPlot(log.cts_real.one, ranks = FALSE)
dev.off()

# As a solution, DESeq2 offers two transformations for counts to stabilize the variance across the means. They give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value. 
# Variance stabilizing transformation (vst function) (Anders & Huber 2010): faster, less sensitive to high count outliers & Better for large datasets (n > 30)
# Regularized-logarithm transformation (rlog function) (Love, Huber and Anders 2014): work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples
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

if (!requireNamespace("dplyr", quietly = TRUE))
    BiocManager::install("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
    BiocManager::install("ggplot2")

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


# 3. Sample distances (similarity between samples):
# Use the R function dist to calculate the Euclidean distance between samples after VST transformation
# Note: dist Function expects the samples to be in rows, thus we have to transpose the input matrix 
sampleDists <- dist(t(assay(vsd)))
# Let us have a look on the distance matrix
sampleDists

# Visualize as a heatmap
if (!requireNamespace("pheatmap", quietly = TRUE))
    install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
    install.packages("RColorBrewer")

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)  ## Transform the dist object into matrix
rownames(sampleDistMatrix) <- paste(vsd$dex, vsd$cell, sep=" - ") ## change row names to contain treatment & cell line instead of sample ID for better visual
colnames(sampleDistMatrix) <- NULL ## remove column names to avoid redundancy 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) ## specify a blue color palette for painting
jpeg('vsd_cts_real_cluster.jpg')
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

## unsupervised transformation
sampleDists2 <- dist(t(assay(vst(ddsTximeta, blind = TRUE))))
sampleDistMatrix2 <- as.matrix(sampleDists2)  ## Transform the dist object into matrix
rownames(sampleDistMatrix2) <- paste(vsd$dex, vsd$cell, sep=" - ") ## change row names to contain treatment & cell line instead of sample ID for better visual
colnames(sampleDistMatrix2) <- NULL ## remove column names to avoid redundancy 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) ## specify a blue color palette for painting
jpeg('vsd_cts_real_cluster2.jpg')
pheatmap(sampleDistMatrix2,
         clustering_distance_rows = sampleDists2,
         clustering_distance_cols = sampleDists2,
         col = colors)
dev.off()


# Alternatively, we can use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. 
# This measure of dissimilarity takes the inherent variance structure of counts into consideration, thus it takes the original count (not normalized)
# As the dist function, this function expects the samples to be in rows, thus we have to transpose the input matrix 
if (!requireNamespace("PoiClaClu", quietly = TRUE))
    install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsTximeta)))
# the distance matrix is stored in poisd$dd
poisd$dd
# Visualize as a heatmap
samplePoisDistMatrix <- as.matrix(poisd$dd)  ## Transform the dist object into matrix
rownames(samplePoisDistMatrix) <- paste(ddsTximeta$dex, ddsTximeta$cell, sep=" - ") ## change row names to contain treatment & cell line instead of sample ID for better visual
colnames(samplePoisDistMatrix) <- NULL ## remove column names to avoid redundancy 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) ## specify a blue color palette for painting
jpeg('pois_cts_real_cluster.jpg')
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.off()

# But how bad it is going to be if we used dist for "raw counts" or "normalized and log2 transformed"
rawSampleDists <- dist(t(counts(ddsTximeta)))
rawSampleDistMatrix <- as.matrix(rawSampleDists)  ## Transform the dist object into matrix
rownames(rawSampleDistMatrix) <- paste(ddsTximeta$dex, ddsTximeta$cell, sep=" - ") ## change row names to contain treatment & cell line instead of sample ID for better visual
colnames(rawSampleDistMatrix) <- NULL ## remove column names to avoid redundancy 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) ## specify a blue color palette for painting
jpeg('raw_cts_real_cluster.jpg')
pheatmap(rawSampleDistMatrix,
         clustering_distance_rows = rawSampleDists,
         clustering_distance_cols = rawSampleDists,
         col = colors)
dev.off()

logSampleDists <- dist(t(log.cts_real.one))
logSampleDistMatrix <- as.matrix(logSampleDists)  ## Transform the dist object into matrix
rownames(logSampleDistMatrix) <- paste(ddsTximeta$dex, ddsTximeta$cell, sep=" - ") ## change row names to contain treatment & cell line instead of sample ID for better visual
colnames(logSampleDistMatrix) <- NULL ## remove column names to avoid redundancy 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) ## specify a blue color palette for painting
jpeg('log_cts_real_cluster.jpg')
pheatmap(logSampleDistMatrix,
         clustering_distance_rows = logSampleDists,
         clustering_distance_cols = logSampleDists,
         col = colors)
dev.off()


### Conclusion:
### Clustering of raw counts is the most missed up. log2 and vsd unsupervised transformation of normalized counts are almost indifferent and comes next
### The best is the Poisson Distance of raw counts or the vsd supervised transformation of normalized counts. The vsd supervised approach made use of the design info so I think the Poisson Distance is better as - I think - is able to adjust for hidden batch effect


# 4. PCA plot
# Another way to visualize sample-to-sample distances 

# With Supervised VST data
# Using the plotPCA function that comes with DESeq2.
jpeg('vsd_cts_real_pca-DESeq.jpg'); plotPCA(vsd, intgroup = c("dex", "cell")); dev.off();
# Using ggplot2
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar")) ## attr function: get or set specific attributes of an object. 
                                                       ## To see all attributes of you object try: attributes(pcaData)
jpeg('vsd_cts_real_pca-ggplot.jpg'); 
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) + 
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
dev.off();

# Alternatively, with Raw counts
# Using Generalized PCA implemented in the glmpca package (Townes et al. 2019). 
# It performs PCA on raw data that is not Normally distributed (e.g. over-dispersed count data)
if (!requireNamespace("glmpca", quietly = TRUE))
    install.packages("glmpca")
library("glmpca")

gpca <- glmpca(counts(ddsTximeta), L=2)
gpca.dat <- gpca$factors
gpca.dat$dex <- ddsTximeta$dex
gpca.dat$cell <- ddsTximeta$cell

jpeg('glmPCA-ggplot.jpg')
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
dev.off();


# 5. multidimensional scaling (MDS): For more info about difference between MDS and PCA, you can check this page: https://stats.stackexchange.com/questions/14002/whats-the-difference-between-principal-component-analysis-and-multidimensional

# With distances calculated from Supervised VST data
# Using cmdscale function in base R. It needs a matrix of distances (like heatmap). Bind to more metadata as well for better visualization 
mdsData <- as.data.frame(cbind(colData(vsd), cmdscale(sampleDistMatrix))) 
jpeg('vsd_cts_real_MDS-ggplot.jpg')
ggplot(mdsData, aes(x = V1, y = V2, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
dev.off();

# With the Poisson Distance. 
mdsPoisData <- as.data.frame(cbind(colData(ddsTximeta), cmdscale(samplePoisDistMatrix)))
jpeg('pois_cts_real_MDS-ggplot.jpg')
ggplot(mdsPoisData, aes(x = V1, y = V2, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
dev.off();

# 6. Gene clustering
# Typically for the most highly variable genes (let us select the top 20)
# We need the normalized expression values to be stabilizied for variance so we will use the VST data.
# Absolute expression itself is not the point, we care about deviation from the gene’s average across all samples (i.e. mean centered)
if (!requireNamespace("genefilter", quietly = TRUE))
    BiocManager::install("genefilter")

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
jpeg('varGeneCluster.jpg')
pheatmap(mat, annotation_col = anno)
dev.off()

#######################################################
# Differential expression analysis

# 1. perform DE by DESeq funtion
# Steps of analysis:
# 1) estimation of size factors (controlling for differences in the sequencing depth of the samples)
#    this can be achieved by calling estimateSizeFactors function ## ddsTximeta <- estimateSizeFactors(ddsTximeta)   
# 2) estimation of dispersion values for each gene. Dispersion describes the size of the distribution of values a variable and can be measured by different statistics, such as range, variance, and standard deviation
#    this can be achived by calling estimateDispersions function ## ddsTximeta <- estimateDispersions(ddsTximeta) 
# 3) fitting a generalized linear model "Negative Binomial GLM fitting and Wald statistics" 
#    this can be achived by calling nbinomWaldTest function ## ddsTximeta <- nbinomWaldTest(ddsTximeta)
# Input/output: the input is DESeqDataSet & the output is DESeqDataSet that contains all the fitted parameters within it
ddsTximeta <- DESeq(ddsTximeta)

# 2. Building the results table
# Function "results" (if used without any arguments) produces the log2 fold changes and p values for the last variable in the design formula
# It will extract these results for a comparison of the last level over the first level.
res <- results(ddsTximeta)
# You can call the results for any variable in the design formula and define the 2 levels you would like to compare 
res2 <- results(ddsTximeta, contrast=c("cell", "N080611", "N061011"))
# alternatively, you can use the "name" argument for calling your target comparison. This is especially useful for interactions
# The value provided to 'name' must be an element of 'resultsNames(DESeq_DE_object)'
resultsNames(ddsTximeta)
res3 <- results(ddsTximeta, name="cell_N61311_vs_N052611")

# there are two ways to be more strict about which set of genes are considered significant: 
# 1. lower the false discovery rate threshold (the threshold on padj) from default (0.1). This value will be used for optimizing the independent filtering
# 2. raise the log2 fold change threshold from default (0). This value will be used for null-hypothesis testing 
res_restricted <- results(ddsTximeta, alpha = 0.05, lfcThreshold=1)

# results correct for multiple testing using Benjamini-Hochberg (BH) adjustment (also called FDR adjustment)
# this can be modified by assigning a different method to the pAdjustMethod argument e.g. bonferroni (see ?p.adjust)
res_rest_Bonf <- results(ddsTximeta, alpha = 0.05, lfcThreshold=1, pAdjustMethod="bonferroni")

# As result is a DataFrame object, it carries metadata with information on the meaning of the columns:
mcols(res, use.names = TRUE)

# Also. we can get overall summery by functions like:
summary(res)
table(res$padj < 0.05) ## or sum(res$padj < 0.05, na.rm=TRUE)

# Notes about p values
# p values could be NA if: 
#    -  all counts for this gene were zero. 
#    -  the gene was excluded from analysis because it contained an extreme count outlier. 
#
# The DESeq2 software (the results function) automatically performs independent filtering that maximizes the number of genes with adjusted p value less than a critical value (by default, alpha is set to 0.1). 
# In any statistical test, we should filter records that do not have enough observations or if the mean value of all observations is very low because these records will not have enough power. 
# This will alleviate the burden on multiple testing. This is exteremly important before Bonfferni correction. But if we are using FDR approach for correction of multiple testing, the non-significant p-values does not affect the correction process. Filteration is important only if removed significant hits that will be more likely to be false positives
# filtering is permissible only if the statistic that we filter on (here the mean of normalized counts across all samples) is independent of the actual test statistic (the p value) under the null hypothesis. 
# Therefore, filteration based on fold change is a baised strategy. DESeq do not actually filter based on fold change, instead it change the hull hypothesis to consider the sig if its expression has more than 2 fold change but eventually the filteration is based on the mean value of all observations


# Typically we select for records with specific threathold of padj. Also, we also can then order to see most down- or up-regulated genes 
resSig <- subset(res, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
# or use the absolute LFC to get the most significant genes 
head(resSig[ order(abs(resSig$log2FoldChange), decreasing = TRUE), ])


# 3. Building the results table after shrinkage of noisy log2 fold changes
# DESeq2 has a function to shrink noisy log2 fold changes. It is called "lfcShrink" & has 3 types of shrinkage estimators
# A commonly used one is apeglm which require a coefficient from the model to shrink, either by name or number as appears in resultsNames(dds).
# Note that this will only change the log2 fold changes not the pValues
# However, it as argument svalue (default false). If true, p-values and adjusted p-values will be replaced with s-values when using 'apeglm' or 'ashr'. 
# s-values provide the probability of false signs among the tests with equal or smaller s-value than a given given's s-value. See Stephens (2016) reference on s-values.
# I need to do more search to figure out if using lfcShrink (with or without s-value) is useful
if (!requireNamespace("apeglm", quietly = TRUE))
    BiocManager::install("apeglm")

library("apeglm")
resultsNames(ddsTximeta)
res.shr <- lfcShrink(ddsTximeta, coef="dex_Dexamethasone_vs_untreated", type="apeglm")

# you can see the change of LFC after shrinkage
jpeg('plot.jpg')
plot(res.shr$log2FoldChange,res$log2FoldChange)
dev.off()

# Note that this will only change the log2 fold changes not the pValues
# Therefore any filtration or visualization involving log2 fold change will show some changes 
resSig.shr <- subset(res.shr, padj < 0.05) ## this will be indifferent in no of genes compared to resSig
head(resSig.shr[ order(resSig.shr$log2FoldChange), ])                      ## This will be different from resSig
head(resSig.shr[ order(resSig.shr$log2FoldChange, decreasing = TRUE), ])


# 4. Building the results table after Independent Hypothesis Weighting
# I need to read the paper of IHW to see if this is useful (https://bioconductor.org/packages/3.12/IHW)
library("IHW")
res.ihw <- results(dds, filterFun=ihw)

#######################################################
# Annotating and exporting results
# In addition to Ensembl gene IDs, we need to add annotations by alternative gene names 
# Bioconductor provides extensive annotation resources
# GenomicFeatures: A set of tools to download and manipulat transcriptomic annotations from UCSC Genome Browser or BioMart DB
# http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
# AnnotationDbi: Implements a user-friendly interface for querying SQLite-based annotation data packages
# http://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html
# To understand more about Bioconductor annotation resources:
# http://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

# We used TxDb.Hsapiens.UCSC.hg38.knownGene which is a Genome centric GenomicFeatures package that map transcripts to genes
# Here we will use org.Hs.eg.db which is a Gene centric AnnotationDbi package that map between gene IDs
# This is the organism annotation package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”), using Entrez Gene IDs (“eg”) as primary key. 
# To get a list of all available key types, use: columns(org.Hs.eg.db)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

# explore the avlaible key types
columns(org.Hs.eg.db)
# get the gene IDs (without the version no) from the results object
ens.str <- substr(rownames(res), 1, 15) 
# use the mapIds function to add individual columns to the results table.
# We provide the row names of our results table as a key, and specify that keytype=ENSEMBL.
# The column argument tells the mapIds function which information we want, and the multiVals argument tells the function what to do if there are multiple possible values for a single input value
res$symbol <- mapIds(org.Hs.eg.db, keys=ens.str, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db, keys=ens.str, column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Sort the results to prepare for export 
resOrdered <- res[order(res$pvalue),]
# convert the results (or subset of it) into dataframe 
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
# Export as a CSV file
write.csv(resOrderedDF, file = "results.csv")


# A more sophisticated way for exporting results the Bioconductor package ReportingTools (Huntley et al. 2013). 
# ReportingTools will automatically generate dynamic HTML documents, including links to external databases using gene identifiers and boxplots summarizing the normalized counts across groups.
if (!requireNamespace("ReportingTools", quietly = TRUE))
    BiocManager::install("ReportingTools")
library("ReportingTools")
# prepare the html folder 
htmlRep <- HTMLReport(shortName="report", title="My report", reportDirectory="./report")
# add the rownames (Ensembl IDs) as a new column if you need to see the Ensembl IDs in the final output html
resOrderedDF$ensembl <- rownames(resOrderedDF)
# read the results table
publish(resOrderedDF, htmlRep)
# write into the HTML file & save the path in a variable 
url <- finish(htmlRep)
# view in the browser if you are not working on a remote or headless machine
browseURL(url)
 


#######################################################
# Plotting results

# 1. Counts plot of specific gene
if (!requireNamespace("ggbeeswarm", quietly = TRUE))
    install.packages("ggbeeswarm")
library("ggbeeswarm")

topGene <- rownames(res)[which.min(res$padj)]
geneCounts <- plotCounts(ddsTximeta, gene = topGene, intgroup = c("dex","cell"), returnData = TRUE)
geneCounts

jpeg('plotCounts-ggplot.jpg')
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off();

jpeg('plotCounts_wLine-ggplot.jpg')
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()
dev.off();


# 2. MA-plot (mean-difference plot or a Bland-Altman plot): 
# On the y-axis, the “M” stands for “minus” – subtraction of log values (this is equivalent to the log of the ratio)
# On the x-axis, the “A” stands for “average”.
jpeg('plotMA.jpg')
plotMA(res, ylim = c(-5, 5))
# We can label individual points on the MA-plot as well.
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

# Alternatively, We can use the results with the shrinked log2 fold changes
jpeg('plotMA_shr.jpg')
plotMA(res.shr, ylim = c(-5, 5))
dev.off()


# 3. Histogram of the p values
# This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
jpeg('histo_pval.jpg')
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
dev.off()

# We can further have a closer look on genes of small p values (say, less than 0.05), to explore the relation to gene expression
qs <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:20/20))
bins <- cut(res$baseMean, qs)
levels(bins) <- paste0("~", round((qs[-1] + qs[-length(qs)])/2, 0))
fractionSig <- tapply(res$pvalue, bins, function(p) mean(p < .05, na.rm = TRUE))
jpeg('pvalFreq_expression.jpg')
barplot(fractionSig, xlab = "mean normalized count", ylab = "fraction of small p values")
dev.off()

# 4. Gene clustering
# Like we did in the Exploratory analysis section but we will select the top 20 DE genes (I will use the output of lfcShrink function)
# Likewise before, We need the normalized expression values to be stabilizied for variance so we will use the VST data.
# Also, we will show the deviation from the gene’s average across all samples (i.e. mean centered)
# We need the normalized expression values to be stabilizied for variance. Here we can use the output of lfcShrink function for specific comparison.
# Absolute expression itself is not the point, we care about deviation from the gene’s average across all samples (i.e. mean centered)
topDEGenes <- head(rownames(resSig.shr[ order(abs(resSig.shr$log2FoldChange), decreasing = TRUE), ]), 20)
mat_de_vsd  <- assay(vsd)[ topDEGenes, ]
mat_de_vsd  <- mat_de_vsd - rowMeans(mat_de_vsd)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
jpeg('deGeneCluster_vsd.jpg')
pheatmap(mat_de_vsd, annotation_col = anno)
dev.off()

#########################################
## Removing hidden batch effects
#  Suppose we did not know that there were different cell lines. This would represent some hidden and unwanted variation
#  We have several packages to identify and correct for such variation

#  sva package: uses the term surrogate variables for this variation (a surrogate variable A variable that can be measured (or easy to measure) that is used in place of one that cannot be measured (or difficult to measure). For example, whereas it may be difficult to assess the wealth of a household,but it is easier to assess the value of a house. 
if (!requireNamespace("sva", quietly = TRUE))
    BiocManager::install("sva")
library("sva")
# we use a full model matrix with the dex variable, and a reduced, or null, model matrix with only an intercept term. 
# Then we specify that we want to estimate 2 surrogate variables.
mod  <- model.matrix(~ dex, colData(ddsTximeta))
mod0 <- model.matrix(~   1, colData(ddsTximeta))
svseq <- svaseq(cts_real, mod, mod0, n.sv = 2)  ## cts_real is the pre-filtered normalized expression data "cts_real <- counts(ddsTximeta, normalized=TRUE)"
# use SVA to remove any effect on the counts from our surrogate variables
# add these two surrogate variables as columns to the DESeqDataSet and then add them to the design
ddssva <- ddsTximeta
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex


#  RUVSeq package: uses the terms factors of unwanted variation
if (!requireNamespace("RUVSeq", quietly = TRUE))
    BiocManager::install("RUVSeq")
library("RUVSeq")
# run DESeq and results to obtain the p-values for the analysis without knowing about the batches i.e.  using design = ~ dex
ddstemp <- ddsTximeta
design(ddstemp) <- ~ dex
ddstemp <- DESeq(ddstemp)
restemp <- results(ddstemp)
# pull out a set of empirical control genes by looking at the genes that do not have a small p-value.
set <- newSeqExpressionSet(counts(ddsTximeta))  ## counts(ddsTximeta) return the pre-filtered raw value expression data
set <- betweenLaneNormalization(set, which="upper") ## This is pre-filtered normalized data (like cts_real) but I am not sure if the upper quantile norm is specifically required 
not.sig <- rownames(restemp)[which(restemp$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ] ## this should be equal to "not.sig" unless we did any additional filteration in the set object
# estimate factors of unwanted variation
set <- RUVg(set, empirical, k=2) ## the 2 estimated factors will be added to the set object as phenoData (run "set" to see the change)
pData(set)
# As before, if we wanted to control for these factors, we simply add them to the DESeqDataSet and to the design:
ddsruv <- ddsTximeta
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + dex

#############################################
## Time course experiments
## Check the online DESeq2 tutorial


