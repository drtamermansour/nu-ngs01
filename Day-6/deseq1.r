# This program can be run as
# cat counts.txt | Rscript deseq1.r NxM (where N is the number of groups and M is the number of replicates per group)

# and produces a table with differentially expressed genes.
# To install the requirements run the program with the 'install` parameter.


# Read the command line arguments.
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("Experimental design must be specified as: NxM at the command line", call.=FALSE)
}

first = args[1]

# Extract the experimental design from the command line.
design = unlist(strsplit(first, 'x'))

# Find the design counts.
cond1_num = as.integer(design[1])
cond2_num = as.integer(design[2])

# Set up the conditions based on the experimental setup.
cond_1 = rep("cond1", cond1_num)
cond_2 = rep("cond2", cond2_num)

# Load the library.
library(DESeq)

# Read the data from the file.
counts = read.table("stdin", header=TRUE, row.names=1, sep="\t")

# Since Kallisto generates the estimated counts as real numbers
# and DESeq 1 allows only integers we need to convert real numbers to integers here.
int_counts = as.matrix(counts)
int_counts = apply(int_counts, 2, function(x) as.integer(round(x)) + 1)
row.names(int_counts) <- row.names(counts)

# Replace with integer counts
counts <- int_counts

# Set up the conditions.
conditions = factor(c(cond_1, cond_2))

# Create a count table
cds = newCountDataSet(counts, conditions)

# Estimate size factors.
cds = estimateSizeFactors(cds)

# Estimate dispersions
cds = estimateDispersions(cds)

# Compute a standard comparison
results = nbinomTest(cds, "cond1", "cond2")

# Sort the results data frame by the padj and foldChange columns.
sorted = results[with(results, order(padj, -foldChange)), ]

# Write the results to the standard output
write.table(sorted, file="", sep="\t", row.name=FALSE, quote=FALSE)

# Keep only the differentiall expressed values.
diffs <- subset(sorted, padj < 0.05, select=c(id))

# Get normalized counts and write to a file
nc = counts(cds, normalized=TRUE )

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc), nc)

# Keep only the rows that were differentially expressed before.
keep = subset(dt, id %in% diffs$id)

# Save into the normalize data matrix.
write.table(keep, file="norm-matrix-deseq1.txt", sep="\t", row.name=FALSE, col.names=TRUE, quote=FALSE)

