library(DESeq2)



# Upload count matrix
counts <- read.table("results/counts/gene_counts.txt", header=TRUE, row.names=1, check.names=FALSE)

# Сheck first 5 rows
head(counts)
