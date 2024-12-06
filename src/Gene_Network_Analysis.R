library(data.table)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(openxlsx)
library(ggplot2)
library(ggpubr)


setwd("/Users/cmdb/Desktop/Quant_Bio_Project/TCGA/")

df_rnaseq <- fread("PANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", header = TRUE)

# Separate gene_id into gene_symbol and gene_number
df_rnaseq <- df_rnaseq %>%
  separate(gene_id, into = c("gene_symbol", "gene_number"), sep = "\\|")


target_gene <- "H2AFY"
# Get the target gene expression levels
target_expression <- as.numeric(df_rnaseq[df_rnaseq$gene_symbol == target_gene, grep("^TCGA", colnames(df_rnaseq))])

# Save gene_symbol separately to reattach later
gene_symbols <- df_rnaseq$gene_symbol

# Get only the columns starting with "TCGA" and transpose
df_transposed <- df_rnaseq[, grep("^TCGA", colnames(df_rnaseq))] %>%
  t() %>%
  as.data.frame()

# Add gene_symbol back as rownames or a column for later reference
colnames(df_transposed) <- gene_symbols

# Filter out columns (genes) with zero variance
df_filtered <- df_transposed[, apply(df_transposed, 2, function(x) var(x, na.rm = TRUE) > 0)]

# Calculate correlations
correlations <- apply(df_filtered, 2, function(gene_expression) {
  cor(gene_expression, target_expression, use = "pairwise.complete.obs")  # Handle NAs properly
})

# Summarize correlations
summary(correlations)

# Get significant genes (absolute correlation > 0.3)
significant_genes <- names(correlations[abs(correlations) > 0.3])
print(head(significant_genes))
significant_values <- correlations[abs(correlations) > 0.3]

# Create a data frame with the significant genes and their correlation values
correlation_table <- data.frame(
  Gene = significant_genes,
  Correlation = significant_values
)

# Sort by absolute correlation values in descending order
correlation_table <- correlation_table %>%
  arrange(desc(abs(Correlation)))

# Ensure significant_genes contains valid gene symbols
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
valid_significant_genes <- intersect(significant_genes, valid_symbols)

write.xlsx(correlation_table, file = "/Users/cmdb/Cancer_Genome_QuantBio_Project/H2AFY_Significant_Genes_Correlation.xlsx", rowNames = FALSE)

# Map significant gene symbols to ENTREZ IDs
entrez_ids <- bitr(valid_significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

print(entrez_ids)


enrich_results <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
print(enrich_results)


# Wrap long GO terms and show top terms
enrich_results@result$Description <- str_wrap(enrich_results@result$Description, width = 30)


dotplot(enrich_results, showCategory = 8) +  # Show top 8 terms
  scale_size(range = c(3, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "GO Biological Processes Enrichment",
    x = "Gene Ratio",
    y = "Biological Process",
    color = "Adjusted p-value",
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 15),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")
  )


# Save the plot with larger dimensions
ggsave(path = "/Users/cmdb/Cancer_Genome_QuantBio_Project/results", filename = "GO_enrichment_for_sig_genes.png", width = 10, height = 12, dpi = 300)



