library(clusterProfiler) # For functional enrichment analysis
library(org.Hs.eg.db)    # Annotation database for human genes
library(AnnotationDbi)   # Annotation utilities

# Read data from a CSV file named "go.csv" and store it in the "data" variable.
data <- read.csv("C:/Users/Emre/Desktop/go.csv", sep = ";", header = TRUE)

# Extract the ENSEMBL gene IDs from the "data" variable and store them in "genes_to_test".
genes_to_test <- data$ensembl_gene_id

# Perform Gene Ontology (GO) enrichment analysis for Biological Process (BP) terms using "enrichGO" function.
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results)

# Create a barplot visualizing the top 15 enriched GO terms for Biological Process.
fit <- plot(barplot(GO_results, showCategory = 15))
fit

# Perform GO enrichment analysis for Molecular Function (MF) terms using "enrichGO" function.
GO_results2 <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
as.data.frame(GO_results2)

# Create a barplot visualizing the top 15 enriched GO terms for Molecular Function.
fit2 <- plot(barplot(GO_results2, showCategory = 15))
fit2

# Perform GO enrichment analysis for Cellular Component (CC) terms using "enrichGO" function.
GO_results3 <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")
as.data.frame(GO_results3)

# Create a barplot visualizing the top 15 enriched GO terms for Cellular Component.
fit3 <- plot(barplot(GO_results3, showCategory = 15))
fit3

