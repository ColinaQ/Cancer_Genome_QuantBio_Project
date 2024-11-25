# How you’ve addressed prior feedback

### Provide IDs/links to specific PDMR and TCGA datasets
We've uploaded the publication that contains the RNAseq data. The rest can be accessed through https://livejohnshopkins-my.sharepoint.com/:f:/g/personal/yhuan279_jh_edu/Egg3oIt1YdNEnxVyetSnKXQBkvdcYdGixwxKkPMXpUFPBA?e=5yaMEl

### Do you have a GSEA package in mind such as http://bioconductor.org/packages/fgsea or https://bioconductor.org/packages/clusterProfiler?
We tried to identify genes co-expressed with H2AFY, annotate their biological pathways, and visualize their interactions using different R packages.

# New progress since last submission 

we stratified patients into groups based on high or low H2AFY expression and generate Kaplan-Meier survival curves for various cancer types. This analysis allows us to evaluate whether H2AFY expression is associated with significant differences in overall survival or progression-free survival.

We attempted the DESeq2 analysis, but remained undecided about whether to continue in this direction or not.

We tried to examine relationship between the expression levels of the gene H2AFY and other genes across multiple cancer samples. Using RNA-seq data from the TCGA dataset, we hope to identify genes that exhibit a strong correlation with H2AFY expression. This analysis serves as a preliminary step for constructing a gene co-expression network, which can provide insights into the functional roles of H2AFY in cancer biology and its interaction with other genes. We are trying different correlation-based methods and integrating databases such as org.Hs.eg.db and STRINGdb. 

# Project Organization

|-- README.md            <-- project title, description, and overview of analysis objectives
|-- requirements.txt     <-- softwares (and versions)
|-- data/                <-- only small (<1 Mb), instructions on retrieving larger
|-- doc/                 <-- reports, slidedecks
|-- results/             <-- only small (<1 Mb), instructions on retrieving larger
|-- src/                 <-- code

# Struggles you are encountering and questions you would like advice on

One of the key struggles we encountered during the project was attempting to preprocess and normalize RNA-seq count data for DESeq2 analysis. Initially, the dataset contained non-integer values, making it incompatible with DESeq2, which requires integer counts for its statistical models. Then, the dataset included missing (NA) and negative values, further complicating preprocessing.

To address these issues, we filterd out missing values and negative entries by setting them as 0. However, since more than 20% of the values are non-integers, we felt that these adjustments risked changing the biological meaning of the data.


