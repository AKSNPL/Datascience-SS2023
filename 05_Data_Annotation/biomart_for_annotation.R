# Load required library for using biomaRt
library(biomaRt)

# Create a connection to the Ensembl database and specify the human gene dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Read data from a CSV file named "mrna.csv" and store it in the "mrna" variable.
mrna <- read.csv("C:/Users/Emre/Desktop/mrna.csv")

# Extract the gene symbols from the "mrna" data and store them in "gene_list".
gene_list <- mrna$Gene_symbol

# Display the contents of the "gene_list" variable, which contains the gene symbols.
gene_list

# Use biomaRt to retrieve the Ensembl gene IDs and external gene names (gene symbols)
# corresponding to the gene symbols in "gene_list".
ensembl_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "external_gene_name",
                     values = gene_list,
                     mart = ensembl)

# Print the retrieved Ensembl gene IDs and their corresponding gene symbols.
print(ensembl_ids)

# Load required library for writing data to Excel files
library(writexl)

# Write the retrieved Ensembl gene IDs and gene symbols to an Excel file named "go.csv".
write_xlsx(ensembl_ids, path = "C:/Users/Emre/Desktop/go.csv")
