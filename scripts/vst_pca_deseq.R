# -------------- Performing VST, PCA and DESeq on R -----------------#


# VST (Variance Stabilizing Transformation):

# - RNA-seq counts has mean-dependent variance
# - PCA requires comparable variance across features
# - VST transforms counts so that:
#  Var ≈ independent of mean

# Used for:
# - PCA
# - clustering
# - QC

# NOT used for:
# - differential expression


# ========================================
# Packages installation
# ========================================

if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}


library(DESeq2)
library(pheatmap)
library(ggplot2)


count_data <- read.table(
  "results/counts/featurecounts_counts.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  comment.char = "#"
)

# =========================================
# Check raw counts
# =========================================

count_data <- count_data[, 6:ncol(count_data)]
colnames(count_data) <- gsub("results.hisat2.|_sorted.bam", "", colnames(count_data))

head(count_data)

# =========================================
# OUTPUT: First 5 rows
#==========================================

#                 SRR11856162 SRR11856164 SRR11856165 SRR11856166

# ENSG00000160072         171         128         127          49
# ENSG00000279928           8          46           1           0
# ENSG00000228037          18           5           1           2
# ENSG00000142611          10           2           0           0
# ENSG00000284616           0           0           0           0
# ENSG00000157911          79          32          23           4

# =========================================

# sum(is.na(count_data))
# sum(count_data < 0)
# str(count_data)


# metadata(colData) creation
coldata <- data.frame(
  row.names = colnames(count_data),
  condition = c("Down", "Healthy", "Healthy", "Down"),
  sex = c("male", "male", "female", "female")
)

coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)

coldata

# DESeqDataSet creation
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = coldata,
  design = ~ sex + condition
)

# Low expressed genes filtration (our dataset is small, to avoid additional noise)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# ==========================================
# Estimate size factors and dispersions
# ==========================================


dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds <- estimateDispersions(dds)

# ===========================================
# VST Transformation
# ===========================================

# blind = False (we have already set the design)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Extract the matrix of variance-stabilized counts, transformed count matrix
vsd_mat <- assay(vsd)

# ===========================================
# Mean/SD check
# ===========================================

# Before VST
raw_stats <- data.frame(
  mean = rowMeans(counts(dds, normalized=TRUE)),
  sd = apply(counts(dds, normalized=TRUE), 1, sd),
  type = "Normalized"
)


# After VST
vst_stats <- data.frame(
  mean = rowMeans(vsd_mat),
  sd = apply(vsd_mat, 1, sd),
  type = "VST"
)

combined_stats <- rbind(raw_stats, vst_stats)


p <- ggplot(combined_stats, aes(mean, sd)) +
  geom_point(alpha = 0.3, size = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  labs(
    x = "Mean expression",
    y = "Standard deviation",
    title = "Mean/SD before vs. after VST"
  )


ggsave(
  filename = "results/plots/mean_sd_before_vs_after_VST.png",
  plot = p,
  width = 10,
  height = 5,
  dpi = 300
)

# Correlation check, must be near 0
cor(vst_stats$mean, vst_stats$sd)

p_loess <- ggplot(combined_stats, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "loess", color = "red") +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  labs(
    x = "Mean expression",
    y = "Standard deviation",
    title = "Mean/SD before vs. after VST with LOESS"
  )

ggsave(
  filename = "results/plots/mean_sd_before_vs_after_VST_loess.png",
  plot = p_loess,
  width = 10,
  height = 5,
  dpi = 300
)

# ==========================================
# PCA
# ==========================================

# PCA expects samples in rows, genes in columns
# t() - > matrix transponition, now rows=samples, columns=genes

# Transpose: rows = samples, columns = genes
vsd_t <- t(vsd_mat)

# Calculate variance
gene_var <- apply(vsd_t, 2, var)

# Remove genes with zero variance, otherwise will throw an error.
vsd_t <- vsd_t[, gene_var > 0]
gene_var <- gene_var[gene_var > 0]  # IMPORTANT

# Keep top 500 most variable genes
top_genes <- order(gene_var, decreasing = TRUE)[1:500]
vsd_t <- vsd_t[, top_genes]

# Perform PCA, center + scale, PCA always must be centered
pca_res <- prcomp(vsd_t, center = TRUE, scale. = TRUE)

# Extract % variance explained, what gene's variability explains more?
explained_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# Create a condition factor (Disease vs Healthy)
sample_info <- data.frame(
  sample = colnames(vsd_mat),
  condition = c("Down", "Healthy", "Healthy", "Down")
)

pca_df <- data.frame(pca_res$x[,1:2], condition = sample_info$condition)

# PCA score plot
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5) +
  xlab(paste0("PC1 (", round(explained_var[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(explained_var[2], 1), "%)")) +
  ggtitle("PCA on VST-transformed counts") +
  theme_bw()

ggsave(
  filename = "results/plots/pca_vst.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# Scree plot
pca_var <- pca_res$sdev^2
explained_var <- pca_var / sum(pca_var) * 100
scree_df <- data.frame(
  PC = paste0("PC", 1:length(explained_var)),
  Variance = explained_var
)

p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = Variance, group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = Variance), color = "red", linewidth = 2) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), 
            vjust = -0.5, size = 3) +
  xlab("Principal Component") +
  ylab("% Variance Explained") +
  ggtitle("Scree Plot: PCA on VST-transformed counts") +
  theme_minimal()

ggsave(
  filename = "results/plots/PCA_scree_VST.png",
  plot = p_scree,
  width = 8,
  height = 5,
  dpi = 300
)

# ==========================================
# DESeq 
# ==========================================


# dds (object, already created)  ~ sex + condition
dds <- DESeq(dds)

# Retrieve results (Down vs Healthy)
res <- results(dds, contrast=c("condition", "Down", "Healthy"))


# Genes annotation (ENSG -> Symbol)

#if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# Add genes symbols to the results
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Sort ans save a table
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), file = "results/DESeq2_results_annotated.csv")


# MA-Plot
png("results/plots/ma_plot.png", width = 800, height = 600)
plotMA(res, ylim=c(-5,5), main="MA-plot: Down vs Healthy")
dev.off()

# Volcano Plot
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

p_volcano <- EnhancedVolcano(res,
                             lab = res$symbol,      
                             x = 'log2FoldChange',
                             y = 'padj',
                             title = 'Down Syndrome vs Healthy',
                             subtitle = 'Differential expression (adjusted p-value < 0.05)',
                             pCutoff = 0.05,
                             FCcutoff = 1.5,
                             pointSize = 3.0,
                             labSize = 4.0,
                             col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                             legendLabels = c('NS', 'Log2 FC', 'Adjusted p-value', 'Both'),
                             drawConnectors = TRUE,
                             widthConnectors = 0.5)

ggsave("results/plots/volcano_plot.png", plot = p_volcano, width = 10, height = 8)

# Heatmap
library(pheatmap)

top30_genes <- head(order(res$padj), 30)
mat <- vsd_mat[top30_genes, ]
rownames(mat) <- res$symbol[top30_genes]

mat <- mat - rowMeans(mat)

df_anno <- as.data.frame(colData(dds)[, c("condition", "sex")])

p_heatmap <- pheatmap(mat, 
                      annotation_col = df_anno,
                      main = "Top 30 Differentially Expressed Genes",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                      border_color = NA,
                      fontsize = 10)

save_pheatmap_pdf <- function(x, filename, width=7, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

png("results/plots/heatmap_top30.png", width = 800, height = 1000)
grid::grid.newpage()
grid::grid.draw(p_heatmap$gtable)
dev.off()
