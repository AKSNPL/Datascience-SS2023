library(readr)
library(data.table)
library(GEOquery)
library(dplyr)
#require(hta20transcriptcluster.db)
###############miRNA1##########################################
#read the top significant genes found by using GEO2R tool:
GSE43732_top_table <- read_delim("GSE43732.top.table.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)
#filter with the help of padj. values
resFilt_GSE43732 <- GSE43732_top_table[which(GSE43732_top_table$adj.P.Val < 0.05 & abs(GSE43732_top_table$logFC) > 0.263), ]


# load series and platform data from GEO
gset1 <- getGEO("GSE43732", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
# choose the index
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
#load the expression data.
t0<- gset1@assayData[["exprs"]]

# extract the information of cancer/normal/adenoma from the metadata.
title0<- gset1@phenoData@data[["characteristics_ch1"]]
#extract the sample names from the expression data.
clname<- colnames(t0)
# make dataframe of the sample names and their corresponding information about cancer/non cancer/adema
d_col<- as.data.frame(title0,clname)
# get metadata out of S4 object
gene_data_frame = fData(gset1)

d<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

d_f<- merge(d,resFilt_GSE43732, by.x="Row.names", by.y= "ID")

duplicated_genes <- d_f$Row.names[duplicated(d_f$Row.names)]


d_f1<- d_f %>% distinct(`Row.names`, .keep_all = T)



###############miRNA2##########################################
#read the top significant genes found by using GEO2R tool:
GSE67269_top_table <- read_delim("GSE67269.top.table.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)
#filter with the help of padj. values
resFilt_GSE67269 <- GSE67269_top_table[which(GSE67269_top_table$adj.P.Val < 0.05 & abs(GSE67269_top_table$logFC) > 0.263), ]


# load series and platform data from GEO
gset1 <- getGEO("GSE67269", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
# choose the index
if (length(gset1) > 1) idx <- grep("GPL19823 ", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[1]]
mi_ids<- as.data.frame(gset1@featureData@data[["miRNA_ID"]])
#load the expression data.
t0<- gset1@assayData[["exprs"]]
h_exprs<- merge(mi_ids,t0,by=0)
h_exprs["Row.names"]<- NULL
head (h_exprs,20)
colnames(h_exprs)[1] <- "miRNA_ID"
h_exprs$miRNA_ID <- gsub("\\*", "", h_exprs$miRNA_ID)
# Aggregate the data by miRNA_ID and take maximum values
aggregated_data <- aggregate(. ~  miRNA_ID, h_exprs, max)

# extract the information of cancer/normal/adenoma from the metadata.
title0<- gset1@phenoData@data[["characteristics_ch1"]]
#extract the sample names from the expression data.
clname<- colnames(aggregated_data[,-1])
# make dataframe of the sample names and their corresponding information about cancer/non cancer/adema
h_col<- as.data.frame(title0,clname)
# get metadata out of S4 object
gene_data_frame = fData(gset1)

h<-aggregated_data

h_f<- merge(h,resFilt_GSE67269, by.x="miRNA_ID", by.y= "miRNA_ID")

duplicated_genes <- h_f$Row.names[duplicated(h_f$Row.names)]

h_f1<- h_f %>% distinct(`miRNA_ID`, .keep_all = T)

#merge all the data with the help of genesymbols.
#ab<- merge(d_f1,h_f1,by.x="Row.names",by.y="miRNA_ID")
####################################################

#clean any duplicate gene symbols:
occurrenceClean <- d_f1[!duplicated(d_f1$Row.names),]
cleany<- occurrenceClean[ , colSums(is.na(occurrenceClean))==0]# remove columns with NA
# remove unnecessary columns
cleanyed<- cleany[, grep("GSM|Row.names", colnames(cleany))]
#again remove more columns
last_cleany<- na.omit(cleanyed)# remove any rows with NA.

fwrite(d_col, file = "miRNA_DS_metadata_col_info.csv", sep = ",",row.names = TRUE)
fwrite(last_cleany, file = "miRNA_DS_preprocessed_data.csv",sep= ",")


######################mRNA############################################################
GSE70409_top_table <- read_delim("GSE70409.top.table.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)
#filter with the help of padj. values
resFilt_GSE70409 <- GSE70409_top_table[which(GSE70409_top_table$adj.P.Val < 0.05 & abs(GSE70409_top_table$logFC) > 0.263), ]


# load series and platform data from GEO
gset1 <- getGEO("GSE70409", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
# choose the index
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
#load the expression data.
t0<- gset1@assayData[["exprs"]]

# extract the information of cancer/normal/adenoma from the metadata.
title0<- gset1@phenoData@data[["characteristics_ch1"]]
#extract the sample names from the expression data.
clname<- colnames(t0)
# make dataframe of the sample names and their corresponding information about cancer/non cancer/adema
a_col<- as.data.frame(title0,clname)

gene_data_frame = fData(gset1)

a<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

a_f<- merge(a,resFilt_GSE70409, by.x="Row.names", by.y= "ID")
a_f1<- a_f %>% distinct(`Gene_symbol`, .keep_all = T)



##############################################################################
GSE20347_top_table <- read_delim("GSE20347.top.table.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)
#filter with the help of padj. values
resFilt_GSE20347 <- GSE20347_top_table[which(GSE20347_top_table$adj.P.Val < 0.05 & abs(GSE20347_top_table$logFC) > 0.263), ]


# load series and platform data from GEO
gset1 <- getGEO("GSE20347", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
# choose the index
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
#load the expression data.
t0<- gset1@assayData[["exprs"]]

# extract the information of cancer/normal/adenoma from the metadata.
title0<- gset1@phenoData@data[["characteristics_ch1"]]
#extract the sample names from the expression data.
clname<- colnames(t0)
# make dataframe of the sample names and their corresponding information about cancer/non cancer/adema
c_col<- as.data.frame(title0,clname)

gene_data_frame = fData(gset1)

cc<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

c_f<- merge(cc,resFilt_GSE20347, by.x="Row.names", by.y= "ID")
c_f1<- c_f %>% distinct(`Gene symbol`, .keep_all = T)


#################do the same for another data##############################################
GSE29001_top_table <- read_delim("GSE29001.top.table.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)
resFilt_GSE29001 <- GSE29001_top_table[which(GSE29001_top_table$adj.P.Val < 0.05 & abs(GSE29001_top_table$logFC) > 0.263), ]

gset1 <- getGEO("GSE29001", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
t0<- gset1@assayData[["exprs"]]

gene_data_frame = fData(gset1)


title0<- gset1@phenoData@data[["characteristics_ch1.1"]]
clname<- colnames(t0)

b_col<- as.data.frame(title0,clname)


b<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

b_f<- merge(b,resFilt_GSE29001, by.x="Row.names", by.y= "ID")
b_f1<- b_f %>% distinct(`Gene symbol`, .keep_all = T)

#############################################
#################do the same for another data##############################################
GSE23400_top_table <- read_delim("GSE23400.top.table.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)
resFilt_GSE23400 <- GSE23400_top_table[which(GSE23400_top_table$adj.P.Val < 0.05 & abs(GSE23400_top_table$logFC) > 0.263), ]

gset1 <- getGEO("GSE23400", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[1]]
t0<- gset1@assayData[["exprs"]]

gene_data_frame = fData(gset1)


title0<- gset1@phenoData@data[["source_name_ch1"]]
clname<- colnames(t0)

g_col<- as.data.frame(title0,clname)


g<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

g_f<- merge(g,resFilt_GSE23400, by.x="Row.names", by.y= "ID")
g_f1<- g_f %>% distinct(`Gene symbol`, .keep_all = T)

###############################################################################
#merge all the data with the help of genesymbols.
ab<- merge(a_f1,b_f1,by.x="Gene_symbol",by.y="Gene symbol")
abc<- merge(ab,c_f1,by.x = "Gene_symbol", by.y= "Gene.symbol")
abcd<- merge(abc,g_f1,by.x = "Gene_symbol", by.y= "Gene.symbol")
############################################################
#order the genesymbol with help of a column:
occurrence <- abcd[order(abcd$Gene_symbol, abcd$GSM1727130, decreasing=TRUE),]
#clean any duplicate gene symbols:
occurrenceClean <- occurrence[!duplicated(occurrence$Gene_symbol),]
cleany<- occurrenceClean[ , colSums(is.na(occurrenceClean))==0]# remove columns with NA
# remove unnecessary columns
cleanyed<- cleany[, grep("GSM|Gene_symbol", colnames(cleany))]
#again remove more columns
last_cleany<- na.omit(cleanyed)# remove any rows with NA.
#bind all metadata cols
binded_col<- rbind(a_col,b_col,c_col,g_col)
#write it and save:
fwrite(binded_col, file = "mRNA_DS_metadata_col_info.csv", sep = ",",row.names = TRUE)
library("readxl")
multi<- read_excel("multimir_results_final.xlsx")
mer<- merge(last_cleany,multi,by.x="Gene_symbol",by.y="target_symbol")
mer1<- mer %>% distinct(`Gene_symbol`, .keep_all = T)
mer1<- mer1[, grep("GSM|Gene_symbol", colnames(mer1))]

#############test_data############################
#GSE38129
GSE38129_top_table <- read_delim("GSE38129.top.table.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)
resFilt_GSE38129 <- GSE38129_top_table[which(GSE38129_top_table$adj.P.Val < 0.05 & abs(GSE38129_top_table$logFC) > 0.263), ]

gset1 <- getGEO("GSE38129", GSEMatrix =TRUE, AnnotGPL=TRUE,getGPL= T)
if (length(gset1) > 1) idx <- grep("GPL6480", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[1]]
t0<- gset1@assayData[["exprs"]]

gene_data_frame = fData(gset1)


title0<- gset1@phenoData@data[["source_name_ch1"]]
clname<- colnames(t0)

t_col<- as.data.frame(title0,clname)


t<-merge(t0,gene_data_frame, by.x=0, by.y= "ID")

t_f<- merge(t,resFilt_GSE38129, by.x="Row.names", by.y= "ID")
t_f1<- t_f %>% distinct(`Gene symbol`, .keep_all = T)
mer11<- t_f1[, grep("GSM|Gene symbol", colnames(t_f1))]
mer11 <- mer11 [1: ncol(mer11)-1 ]
mer11<-na.omit(mer11)


fwrite(t_col, file = "mRNA_DS_metadata_col_test_info.csv", sep = ",",row.names = TRUE)
fwrite(mer11, file = "mRNA_DS_test_data.csv",sep= ",")
mer12<- merge(x =mer1 , y = mer11, by.x = "Gene_symbol", by.y="Gene symbol", all.y = TRUE)
merged_data <- mer12[, names(mer1)]
merged_data<-na.omit(merged_data)
fwrite(merged_data, file = "mRNA_DS_preprocessed_training_data.csv",sep= ",")
# see distribution of one of deuplicated genes.
########################################################
# Install required packages
#Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. 
#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(survminer)
library(survival)

# getting clinical data for TCGA-BRCA cohort -------------------
clinical_escc <- GDCquery_clinic("TCGA-ESCA")
any(colnames(clinical_escc) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_escc) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_escc[,c(8,42,47)]
# looking at some variables associated with survival 
table(clinical_escc$vital_status)

# change certain values the way they are encoded
clinical_escc$deceased <- ifelse(clinical_escc$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_escc$overall_survival <- ifelse(clinical_escc$vital_status == "Alive",
                                         clinical_escc$days_to_last_follow_up,
                                         clinical_escc$days_to_death)

# build a query to get gene expression data for entire cohort
query_escc_all = GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  access = "open")

output_escc <- getResults(query_escc_all)
# get 20 primary tissue sample barcodes
tumor <- output_escc$cases#[1:50]

# # get gene expression data from 20 primary tumors 
query_escc <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor)

# download data
GDCdownload(query_escc)
library(SummarizedExperiment)
# get counts
tcga_escc_data <- GDCprepare(query_escc, summarizedExperiment = T)
escc_matrix <- assay(tcga_escc_data)
escc_matrix[1:10,1:10]

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_escc_data))
coldata <- as.data.frame(colData(tcga_escc_data))


#####################################################
lis<- escc_matrix
comm<- merge(lis,gene_metadata,by.x=0,by.y="gene_id")
comm<- comm[, grep("TCGA|gene_name", colnames(comm))]
test_tcga<- merge(comm, mer1,by.x="gene_name",by.y="Gene_symbol",all.y= T)
test_tcga<- test_tcga[, grep("TCGA|gene_name", colnames(test_tcga))]
sample_ids<- colnames(test_tcga)

# Function to assign groups based on TCGA IDs
# Function to assign groups
assign_group <- function(tcga_id) {
  group <- ""
  parts <- unlist(strsplit(tcga_id, "-"))
  
  if (length(parts) >= 4) {
    fourth_part <- parts[4]
    if (grepl("\\d{2}", fourth_part)) {
      num <- as.numeric(substr(fourth_part, 1, 2))
      if (num >= 10 & num <= 29) {
        group <- "Control"
      } else if (num >= 1 & num <= 9) {
        group <- "Cancer"
      }
    }
  }
  
  return(group)
}

# Assign groups to TCGA IDs
group_assignments <- sapply(sample_ids, assign_group)

comb<- as.data.frame(group_assignments)
#comb<- cbind(c_data,sample_labels)

#fwrite(test_tcga, file = "mRNA_TCGA_DS_test_data.csv",sep= ",")
#fwrite(comb, file = "mRNA_TCGA_DS_col_data.csv",sep= ",",row.names = T)


##############################################


# vst transform counts to be used in survival analysis ---------------
library(DESeq2)
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = escc_matrix,
                              colData = coldata,
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
######################################################################
#

################################################################

# vst 
vsd <- vst(dds, blind=FALSE)
escc_matrix_vst <- assay(vsd)
escc_matrix_vst[1:10,1:10]
library(tidyr)
library(dplyr)
library(tibble)
# Get data for TP53 gene and add gene metadata information to it -------------
escc_gene <- escc_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "RUVBL1")
# get median value
median_value <- median(escc_gene$counts)

# denote which cases have higher   lower expression than median count
escc_gene$strata <- ifelse(escc_gene$counts >= median_value, "HIGH", "LOW")

# Add clinical information to escc_gene
escc_gene$case_id <- gsub('-01.*', '', escc_gene$case_id)
escc_gene <- merge(escc_gene, clinical_escc, by.x = 'case_id', by.y = 'submitter_id')
# Convert days to months for overall_survival variable
escc_gene$overall_survival <- escc_gene$overall_survival / 30

# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = escc_gene)
fit
ggsurvplot(fit,
           data = escc_gene,
           surv.median.line = "hv",
           test.for.trend = FALSE,
           risk.table = T)


fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = escc_gene)