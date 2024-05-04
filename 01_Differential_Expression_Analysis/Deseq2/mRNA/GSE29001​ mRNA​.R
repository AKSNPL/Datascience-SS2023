library(GEOquery)
library(tidyverse)
library(magrittr)

mRNA_GSE29001 <- getGEO("GSE29001", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(mRNA_GSE29001) > 1) idx <- grep("GPL4508", attr(mRNA_GSE29001, "names")) else idx <- 1
mRNA_GSE29001 <- mRNA_GSE29001[[idx]]
GSE29001 <- mRNA_GSE29001@assayData[["exprs"]]
colx <- mRNA_GSE29001@featureData@data[["Gene symbol"]]
GSE29001 <- cbind(colx,GSE29001)

# Convert matrix to a dataframe
df_GSE29001 <- as.data.frame(GSE29001)

library(DESeq2)

# Extract only the count matrix
count_matrix <- df_GSE29001[, 2:ncol(df_GSE29001)]

# Sample information
sample_names <- colnames(count_matrix)
conditions = c("normal", "normal", "tumor", "normal", "normal", "tumor", "tumor", "normal", "normal", "tumor", "tumor",
                            "normal", "normal", "tumor", "tumor", "normal", "normal", "tumor", "tumor", "normal", "normal",  "tumor", 
                            "normal", "normal", "tumor", "tumor", "normal", "normal", "tumor", "tumor", "normal", "normal", "tumor", 
                            "tumor", "normal", "normal", "tumor", "tumor", "normal", "normal", "tumor", "tumor", "normal", "normal", "tumor")

sample_info <- data.frame(Sample = sample_names, condition = conditions)

#### SKIRMISH FOR DESEQ2 DATASET ###########
count_matrix <- as.matrix(sapply(count_matrix, as.numeric))
count_matrix <- count_matrix[complete.cases(count_matrix), ]
count_matrix <- round(count_matrix)

#get rownames of previous data
row_names_GSE29001 <- rownames(GSE29001)

# Set gene_id to rownames again
rownames(count_matrix) <- row_names_GSE29001

# Create the DESeqDataSet using the updated sample_info
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)

# Perform differential gene expression analysis
dds <- DESeq(dds)

resultsNames(dds)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

res <- results(dds, name="condition_tumor_vs_normal")
summary(res)

# Adjustment to overcome NA-error --> exclude rows containing NA's
filtered_indices <- which(res$log2FoldChange > 1 & res$padj < 0.05 & !is.na(res$padj))
filtered_downregulated_genes <- res[filtered_indices, ]
gene_names <- rownames(filtered_downregulated_genes)
gene_names

rld <- rlog(dds, blind=TRUE) 
PCA1 <- plotPCA(rld, intgroup = "condition")
PCA1

vst <- assay(vst(dds))
p <- pca(vst, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)

biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
