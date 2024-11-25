library(data.table)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)

# Use BiocManager to install packages
#BiocManager::install("")

setwd("/Users/cmdb/Desktop/Quant_Bio_Project/TCGA/")

df_rnaseq <- fread("PANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", header = TRUE)
df_rnaseq <- df_rnaseq %>%
  separate(gene_id, into = c("gene_symbol", "gene_number"), sep = "\\|")

target_gene <- "H2AFY"

target_expression <- df_rnaseq %>%
  filter(gene_symbol == target_gene) %>%
  dplyr::select(starts_with("TCGA")) %>%
  unlist()      

df_transposed <- df_rnaseq %>%
  dplyr::select(starts_with("TCGA")) %>%
  t() %>%  # Transpose the data
  as.data.frame()  # Convert to a data frame
print(dim(df_transposed))

df_filtered <- df_transposed[, apply(df_transposed, 2, var, na.rm = TRUE) > 0]

correlations <- apply(df_filtered, 2, function(gene_expression) {
  cor(gene_expression, target_expression, use = "complete.obs")
})

summary(correlations)

significant_genes <- names(correlations[abs(correlations) > 0.5])


print(head(df_rnaseq[[1]]))

entrez_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)








