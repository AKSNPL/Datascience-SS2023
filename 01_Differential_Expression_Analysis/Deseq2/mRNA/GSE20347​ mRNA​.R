library(GEOquery)
library(tidyverse)
library(magrittr)

mRNA_GSE20347 <- getGEO("GSE20347", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(mRNA_GSE20347) > 1) idx <- grep("GPL4508", attr(mRNA_GSE20347, "names")) else idx <- 1
mRNA_GSE20347 <- mRNA_GSE20347[[idx]]
GSE20347<- mRNA_GSE20347@assayData[["exprs"]]
colx<- mRNA_GSE20347@featureData@data[["Gene symbol"]]
GSE20347<- cbind(colx,GSE20347)

# Convert matrix to a dataframe
df_GSE20347 <- as.data.frame(GSE20347)

library(DESeq2)

# Extract only the count matrix
count_matrix <- df_GSE20347[, 2:ncol(df_GSE20347)]

# Sample information
sample_names <- colnames(count_matrix)

# Create the sample_info data frame
conditions <- c(rep("normal adjacent esophageal tissue", 17), rep("esophageal squamous cell carcinoma", length(sample_names) - 17))
sample_info <- data.frame(Sample = sample_names, Condition = conditions)

#### SKIRMISH FOR DESEQ2 DATASET ###########
count_matrix <- as.matrix(sapply(count_matrix, as.numeric))
count_matrix <- count_matrix[complete.cases(count_matrix), ]
count_matrix <- round(count_matrix)

#get rownames of previous data
row_names_GSE20347 <- rownames(GSE20347)

# Set gene_id to rownames again
rownames(count_matrix) <- row_names_GSE20347

# Create the DESeqDataSet using the updated sample_info
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Condition)

# Perform differential gene expression analysis
dds <- DESeq(dds)

resultsNames(dds)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

res <- results(dds, name="Condition_normal.adjacent.esophageal.tissue_vs_esophageal.squamous.cell.carcinoma")
summary(res)

# Adjustment to overcome NA-error --> exclude rows containing NA's
filtered_indices <- which(res$log2FoldChange < 1 & res$log2FoldChange > -1 & res$padj < 0.05 & !is.na(res$padj))
filtered_downregulated_genes <- res[filtered_indices, ]
gene_names <- rownames(filtered_downregulated_genes)
gene_names


res2 <- as.data.frame(results(dds, name = "Condition_normal.adjacent.esophageal.tissue_vs_esophageal.squamous.cell.carcinoma"))

# Convert p-value to -log10(p-value)
res2$log10_pvalue <- -log10(res2$pvalue)

# Create the volcano plot
volcano_plot <- ggplot(data = res2, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point() +
  labs(x = "log2 Fold Change", y = "-log10 p-value") +
  theme_minimal()

# Display the volcano plot
print(volcano_plot)

# Add vertical lines for fold change thresholds
volcano_plot <- volcano_plot +
  geom_vline(xintercept = c(-1, 1), color = "red") +
  geom_hline(yintercept = -log10(0.05), color = "red")

# Display the customized volcano plot
print(volcano_plot)



rld <- rlog(dds, blind=TRUE) 
PCA1 <- plotPCA(rld, intgroup = "Condition")
PCA1

vst <- assay(vst(dds))
p <- pca(vst, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
