#
# This program can be run as
#
# cat counts.txt | Rscript deseq2.r
#
# Produces a table with differentially expressed genes. on the standard output.
#
# To install the requirements run the program with the 'install` parameter.


# Read the command line arguments.
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("Experimental design must be specified as: NxM at the command line", call.=FALSE)
}

first = args[1]

if (first == 'install') {
    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    stop("Installation completed", call.=FALSE)
}

# Load the library.
library(DESeq2)

# Extract the experimental design from the command line.
design = unlist(strsplit(first, 'x'))

# Find the desing counts.
cond1_num = as.integer(design[1])
cond2_num = as.integer(design[2])

# Set up the conditions based on the experimental setup.
cond_1 = rep("cond1", cond1_num)
cond_2 = rep("cond2", cond2_num)

# Read the data from the standard input.
countData = read.table("stdin", header=TRUE, sep="\t", row.names=1 )

# Build the dataframe from the conditions
samples = names(countData)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)

#coldata = read.table(coldata_file, header=TRUE, sep="\t", row.names=1 )

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

#Set the reference to be compared
dds$condition = relevel(dds$condition,"cond1")

# Run deseq
dds = DESeq(dds)

# Format the results.
res = results(dds)

# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]

# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)

# Write the table out.
write.table(sorted.df, file="", sep="\t", col.names=NA, quote=FALSE)

# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

# Save the normalize data matrix.
write.table(dt, file="norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)


