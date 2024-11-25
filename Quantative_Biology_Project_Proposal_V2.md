# Quantative Biology Project Proposal

###### Author: Colina Qiao; Yixuan Huang

### Title

Exploring the Relationship Between Histone Variant Expression and Cancer Survival Rates Using TCGA Data

### Description

This project will investigate the association between histone variant gene expression levels and cancer patient survival outcomes. Using The Cancer Genome Atlas (TCGA), the study focuses on abnormal expression patterns of histone variants in different cancers and their correlation with clinical outcomes. Analytical methods include survival analysis and the visualization of expression patterns across cancer types to identify significant trends. The ultimate goal is to uncover potential prognostic markers and therapeutic targets related to histone variants.

### Example Published Figure

######  Exploring the Prognostic Value, Immune Implication and Biological Function of H2AFY Gene in Hepatocellular Carcinoma

doi: 10.3389/fimmu.2021.723293

FIGURE 1 | The elevated H2AFY expression in HCC

![img](https://www.frontiersin.org/files/Articles/723293/fimmu-12-723293-HTML/image_m/fimmu-12-723293-g001.jpg)

FIGURE 2 | H2AFY expression in sub-groups of different clinical characteristics.

![img](https://www.frontiersin.org/files/Articles/723293/fimmu-12-723293-HTML/image_m/fimmu-12-723293-g002.jpg)


### Datasets

- **TCGA (The Cancer Genome Atlas)**: A large database of genomic data across multiple cancer types, including histone variant expression.


### Software

- **R **: https://www.r-project.org/ Version 2024.04.2 

- **Python(possible):** https://www.python.org/downloads/ Version: Python 3.12.3 

### Proposed Steps

#### Step 1: Data Preprocessing and Extraction

Extract histone variant gene expression data from TCGA.
Clean RNA-seq data for statistical analysis.
Merge RNA-seq data with clinical and survival metadata.

#### Step 2: Expression Analysis of Histone Variants
Identify and analyze expression levels of histone variants across tumor and normal samples.
Create visualizations, such as boxplots, to compare expression patterns across cancer types.

#### Step 3: Kaplan-Meier Survival Analysis
Perform Kaplan-Meier survival analysis for selected histone variants (e.g., H2AFY, H1F0).
Compare survival outcomes for patient groups with high versus low gene expression levels.
Generate cancer-specific Kaplan-Meier survival curves and assess statistical significance using log-rank tests.

#### Step 4: Exploratory Gene Network Analysis
Identify co-expressed genes for significant histone variants.
Investigate potential functional implications using pathway enrichment analysis.

