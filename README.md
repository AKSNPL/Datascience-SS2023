# Datascience-SS2023



## Introduction
All results from the project "Integrative Analysis of miRNA and mRNA Expression Profiles in Squamous Cell Carcinoma" are stored in this GitLab repository  
with the corresponding code snippets.

Esophageal Squamous Cell Carcinoma (ESCC) is a type of cancer that originates in the squamous cells lining the esophagus, the muscular tube connecting the throat to the stomach. It is one of the major subtypes of esophageal cancer and is known for its relatively limited treatment options. Early detection and identification of key genes involved in ESCC are crucial for improving patient outcomes and developing targeted therapies.\
One promising way of investigation is the use of miRNA-mRNA interaction networks, which are used to explore the interactions between microRNAs (miRNAs) and messenger RNAs (mRNAs).\
MicroRNAs are small non-coding RNA molecules that play a crucial role in gene regulation. They can bind to specific mRNA sequences, leading to either the degradation or translational repression of the target mRNA. By understanding this complicated miRNA-mRNA regulatory network, we hope to identify key genes that contribute to the development and progression of ESCC.\
This study aims to construct a specific miRNA-mRNA interaction network for ESCC by integrating multi-omics data from patient samples and using machine learning techniques such as SVM and XGBoost to identify key genes, which may have an influence to the pathogenesis of ESCC. A comprehensive analysis of these key genes may provide valuable insights into the underlying mechanisms of the disease and ultimately lead to the development of new diagnostic tools and personalized treatment strategies for patients with ESCC.

## Folder explanation: 
We separated our folders like our workflow in 6 parts: <br />
01_Differential_Expression_Analysis: DGE Analysis and the merging results  <br />
02_MultiMIR : miRNA to target mRNA results  <br />
03_Machine_Learning: Machine learning notebooks <br />
04_Cytoscape : miRNA-mRNA interaction using Cytoscape  <br />
05_Data_Annotation : Functional annotation results  <br />
06_Survival_Analysis: survival analysis scripts <br />