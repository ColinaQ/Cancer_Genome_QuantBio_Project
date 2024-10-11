# Quantative Biology Project Proposal

###### Author: Colina Qiao; Yixuan Huang

### Title

Investigating the Relationship Between Histone Variant Expression and Cancer Survival Rates Using PDMR and TCGA Data 

### Description

This project will analyze the association between specific histone variant gene expression levels and cancer patient survival rates. Employing the National Cancer Institute’s Patient-Derived Models Repository (PDMR) and The Cancer Genome Atlas (TCGA), we aim to identify abnormal expression patterns in cancerous tissues and correlate these with patient outcomes. The project will focus on cancers with available histone variant data, using statistical and bioinformatics tools to generate Kaplan-Meier survival curves and conduct correlation analyses. 

### Example Published Figure

######  Exploring the Prognostic Value, Immune Implication and Biological Function of H2AFY Gene in Hepatocellular Carcinoma

doi: 10.3389/fimmu.2021.723293

FIGURE 1 | The elevated H2AFY expression in HCC

![img](https://www.frontiersin.org/files/Articles/723293/fimmu-12-723293-HTML/image_m/fimmu-12-723293-g001.jpg)

FIGURE 2 | H2AFY expression in sub-groups of different clinical characteristics.

![img](https://www.frontiersin.org/files/Articles/723293/fimmu-12-723293-HTML/image_m/fimmu-12-723293-g002.jpg)

### Datasets

- **PDMR (Patient-Derived Models Repository)**: A resource for molecular and clinical data of patient-derived cancer models.

- **TCGA (The Cancer Genome Atlas)**: A large database of genomic data across multiple cancer types, including histone variant expression.

### Software

- **R **: https://www.r-project.org/ Version 2024.04.2 

- **DESeq2 **: for differential expression analysis. https://bioconductor.org/packages/release/bioc/html/DESeq2.html 

- **Python(possible):** https://www.python.org/downloads/ Version: Python 3.12.3 

### Proposed Steps

##### Step 1: Data Downloading and Preprocessing 

Extract and sort RNA-seq data for histone variants from PDMR and TCGA for selected cancer types (e.g., lung cancer, hepatocellular carcinoma, pancreatic cancer). 

Normalize the expression data using DESeq2 to prepare for differential expression analysis. 

##### Step 2: Experssion analysis in different cancer types

Using DESeq2, conduct differential expression analysis to compare histone variant expression levels in tumor vs. normal tissues. Generate volcano plots to visualize differentially expressed genes.

##### Step 3: Kaplan-Meier Survival Analysis 

Perform Kaplan-Meier survival analysis comparing high vs. low expression groups for selected histone variants (e.g., H2AZ1, H2AFY, H1.3). 

Visualize survival curves and compute statistical significance using log-rank tests. 

##### Step 4: Integration and Functional Analysis of Co-Expressed Genes 

Conduct gene set enrichment analysis (GSEA) and pathway analysis on histone variant co-expression networks to explore functional implications. 

Identify potential drug targets using cBioPortal or STRING database for interactions related to survival outcomes. 