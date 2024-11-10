# How you’ve addressed prior feedback
### Is there one (or more) histone variants that you would target first?
Yes, we will focus on H2AFY first, since they have known roles in chromatin remodeling, gene regulation, and implications in cancer progression. In future analyses, we will look at other histones and their variants more extensively.
### Provide IDs/links to specific PDMR and TCGA datasets
So far we have only RNA-seq expression and clinical metadata from the TCGA Pan-Cancer Atlas. They can be accessed through https://livejohnshopkins-my.sharepoint.com/:f:/g/personal/yhuan279_jh_edu/Egg3oIt1YdNEnxVyetSnKXQBkvdcYdGixwxKkPMXpUFPBA?e=5yaMEl
### Do you have a GSEA package in mind such as http://bioconductor.org/packages/fgsea or https://bioconductor.org/packages/clusterProfiler?
We may use GSEA package to analyze enrichment in key pathways influenced by H2AFY or other histones and identify gene sets linked to cancer progression and survival outcomes.

# New progress since last submission 
## Obtained the datasets from the TCGA Pan-Cancer Atlas
## Processed the datasets and focused on H2AFY
- Filtered and reshaped RNA-seq data for H2AFY expression and merged it with survival metadata.
- Created and categorized expression levels into high and low groups based on the median H2AFY expression.
- Split the data by cancer type to analyze the relationship between H2AFY expression and survival outcomes individually for each cancer type.
- Generated Kaplan-Meier survival curves with p-values to visually and statistically assess survival differences based on H2AFY expression levels for each cancer type.
## Started to analyze H1F0 expression
# Project Organization

|-- README.md            <-- project title, description, and overview of analysis objectives
|-- requirements.txt     <-- softwares (and versions)
|-- data/                <-- only small (<1 Mb), instructions on retrieving larger
|-- doc/                 <-- reports, slidedecks
|-- results/             <-- only small (<1 Mb), instructions on retrieving larger
|-- src/                 <-- code

# Struggles you are encountering and questions you would like advice on
- We've encountered challenges with missing values in death_days_to and last_contact_days_to columns. Although we've implemented a workaround by treating missing death_days_to as censored if vital_status is "Alive," we are still not sure if this approach fully accounts for all potential data gaps.
- Some histone variants are highly expressed across almost all cancer types in our datasets, but the surviving curves may not indicate any difference. We want to find a biological explaination for this.
