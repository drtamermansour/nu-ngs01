
# Read the command line arguments.
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("Experimental design must be specified as: NxM at the command line", call.=FALSE)
}

first = args[1]

if (first == 'install') {
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    stop("Installation completed", call.=FALSE)
}

library(edgeR)

# Extract the experimental design from the command line.
design = unlist(strsplit(first, 'x'))

# Find the desing counts.
cond1_num = as.integer(design[1])
cond2_num = as.integer(design[2])

# Set up the conditions based on the experimental setup.
cond_1 = rep("cond1", cond1_num)
cond_2 = rep("cond2", cond2_num)

# Read the input table.
counts = read.table("stdin", header=TRUE, sep="\t", row.names=1 )

# Create the groups.
group=c(cond_1, cond_2)

dge <- DGEList(counts=counts, group=group)
dis <- estimateCommonDisp(dge)
tag <- estimateTagwiseDisp(dis)

# This performs a pairwise comparison.
etx <- exactTest(tag)
etp <- topTags(etx, n=100000, sort.by="p.value")

# Generate the output.
write.table(etp$table, file="", sep="\t", row.name=TRUE, quote=FALSE)

# Get normalized counts and write to a file
scale = dge$samples$lib.size*dge$samples$norm.factors
nc = round(t(t(counts)/scale)*mean(scale))

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

# Save into the normalize data matrix.
write.table(dt, file="norm-matrix-edgeR.txt", sep="\t", row.name=FALSE, col.names=TRUE, quote=FALSE)

