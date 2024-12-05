library("data.table")
library("tidyverse")
library("ggthemes")
library("dplyr")
library("readr")
library("survival")
library("survminer")
library("ggplot2")

setwd("/Users/cmdb/Desktop/Quant_Bio_Project/TCGA/")
df_merged_sample <- fread("merged_sample_quality_annotations.tsv")

df_rnaseq <- fread("PANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")

df_survival <- fread("Survival_SupplementalTable_S1.tsv")

df_rnaseq <- df_rnaseq %>%
  separate(gene_id, into = c("gene_symbol", "gene_number"), sep = "\\|")


# Step 1: Prepare Data for Box Plot
# Make sure `df_long_H2AFY` is already created and has `expression_level`, `sample_id`, and `cancer type abbreviation`
# Check column names in each data frame
df_H2AFY <- df_rnaseq %>%
  filter(gene_symbol == "H2AFY")

df_long_H2AFY <- df_H2AFY %>%
  pivot_longer(
    cols = starts_with("TCGA"),
    names_to = "sample_id",
    values_to = "expression_level"
  )

df_long_H2AFY <- df_long_H2AFY %>%
  mutate(sample_id_trimmed = substr(sample_id, 1, 15))  # Extracts the first 15 characters

colnames(df_long_H2AFY)
colnames(df_merged_sample)

head(df_long_H2AFY$sample_id_trimmed)
head(df_merged_sample$patient_barcode)
head(df_merged_sample$aliquot_barcode)

# Trim sample_id_trimmed in df_long_H2AFY to the first 12 characters
df_long_H2AFY <- df_long_H2AFY %>%
  mutate(patient_id = substr(sample_id_trimmed, 1, 12))

# Trim aliquot_barcode in df_merged_sample to the first 12 characters to get a patient-level ID
df_merged_sample <- df_merged_sample %>%
  mutate(patient_id = substr(aliquot_barcode, 1, 12))

df_H2AFY_with_cancer <- df_long_H2AFY %>%
  inner_join(df_merged_sample, by = "patient_id") %>%
  filter(!is.na(`cancer type`))  # Filter out rows with missing cancer type info


# Enhanced Box Plot for H2AFY Expression Across Cancer Types
ggplot(df_H2AFY_with_cancer, aes(x = `cancer type`, y = expression_level, fill = `cancer type`)) +
  geom_boxplot(outlier.shape = NA) +
#  geom_jitter(width = 0.2, alpha = 0.4, color = "black", size = 1) +
  scale_y_continuous(
    trans = 'log10',  # Log scale for better visibility
    limits = c(300, 30000),
    breaks = c(300, 1000, 3000, 10000, 30000)
  ) +
  labs(
    title = "H2AFY Expression Across Different Cancer Types",
    x = "Cancer Type",
    y = "H2AFY Expression Level (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    legend.position = "none",  # Hide legend
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.grid.major.y = element_line(color = "lightgray", size = 0.5)  # Light gray horizontal lines
  ) +
  # Create a color palette based on the number of unique cancer types
  scale_fill_manual(values = colorRampPalette(ggthemes::tableau_color_pal("Tableau 20")(20))(
    length(unique(df_H2AFY_with_cancer$`cancer type`))
  )) +
  # Add a horizontal line at the median expression level
  geom_hline(yintercept = median(df_H2AFY_with_cancer$expression_level, na.rm = TRUE), linetype = "dashed", color = "red", size = 0.7)


# -------------------------------------------------------------
# find all the histone variant in the dataset
# -------------------------------------------------------------


# Define a regular expression to match common histone variant patterns
histone_variant_pattern <- "^(H1|H2A|H2B|H3|H4)[A-Z0-9]*"

# Filter for histone variants in the dataset
df_histone_variants <- df_rnaseq %>%
  filter(grepl(histone_variant_pattern, gene_symbol))

# Display the unique histone variant gene symbols found
unique_histone_variants <- unique(df_histone_variants$gene_symbol)
print(unique_histone_variants)

# -------------------------------------------------------------
# try to analyze more histone related genes 
# -------------------------------------------------------------

# Step 1: Define a List of Histone-Related Genes
histone_genes <- unique_histone_variants
# Step 2: Filter RNAseq Data for These Histone Genes
df_histones <- df_rnaseq %>%
  filter(gene_symbol %in% histone_genes)

# Step 3: Pivot Data to Long Format
df_long_histones <- df_histones %>%
  pivot_longer(
    cols = starts_with("TCGA"),
    names_to = "sample_id",
    values_to = "expression_level"
  ) %>%
  mutate(sample_id_trimmed = substr(sample_id, 1, 15))  # Trim sample_id to match survival data

# Check the first few entries in sample_id_trimmed in df_long_histones
head(df_long_histones$sample_id_trimmed)

# Check the first few entries in patient_id created from aliquot_barcode in df_merged_sample
df_merged_sample <- df_merged_sample %>%
  mutate(patient_id = substr(aliquot_barcode, 1, 12))
head(df_merged_sample$patient_id)

# Create a new column in df_long_histones to match patient_id format
df_long_histones <- df_long_histones %>%
  mutate(patient_id = substr(sample_id_trimmed, 1, 12))

# Join on patient_id
df_histones_with_cancer <- df_long_histones %>%
  inner_join(df_merged_sample, by = "patient_id") %>%
  filter(!is.na(`cancer type`))

# Verify the result
nrow(df_histones_with_cancer)  # Should be greater than 0 if the join is successful
head(df_histones_with_cancer)   # Inspect the first few rows


# Loop through each histone gene and create a box plot
for (gene in histone_genes) {
  # Filter for the current gene and remove non-positive expression values
  df_gene <- df_histones_with_cancer %>%
    filter(gene_symbol == gene, expression_level > 0)  # Keep only positive expression values
  
  # Calculate the median expression level, ensuring NA values are removed
  median_expression <- median(df_gene$expression_level, na.rm = TRUE)
  
  # Generate the box plot without individual points
  p <- ggplot(df_gene, aes(x = `cancer type`, y = expression_level, fill = `cancer type`)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(
      trans = 'log10',
      limits = c(300, 30000),
      breaks = c(300, 1000, 3000, 10000, 30000)
    ) +
    labs(
      title = paste("Expression of", gene, "Across Different Cancer Types"),
      x = "Cancer Type",
      y = paste(gene, "Expression Level (log scale)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "lightgray", size = 0.5)
    ) +
    scale_fill_manual(values = colorRampPalette(ggthemes::tableau_color_pal("Tableau 20")(20))(
      length(unique(df_gene$`cancer type`))
    ))
  
  # Conditionally add geom_hline if median is within the y-axis range
  if (!is.na(median_expression) && median_expression >= 300 && median_expression <= 30000) {
    p <- p + geom_hline(yintercept = median_expression, linetype = "dashed", color = "red", size = 0.7)
  }
  
  # Explicitly print the plot to display it in the loop
  print(p)
}

# H1F0, H3f3B have high expression level across different cancer type



###################################################################
#Analysis on H2AF2


df_high_exp_survival_H2AFY <- df_long_H2AFY %>%
  inner_join(df_survival, by = c("sample_id_trimmed" = "sample"))

median_exp_H2AFY <- median(df_high_exp_survival_H2AFY$expression_level, na.rm = TRUE)

df_high_exp_survival_H2AFY <- df_high_exp_survival_H2AFY %>%
  mutate(expression_group = ifelse(expression_level >= median_exp_H2AFY, "High", "Low"))

df_high_exp_survival_H2AFY <- df_high_exp_survival_H2AFY %>%
  mutate(
    death_days_to = as.numeric(death_days_to),  # Convert to numeric
    last_contact_days_to = as.numeric(last_contact_days_to),  # Ensure numeric conversion for last_contact_days_to
    survival_time = ifelse(is.na(death_days_to) & vital_status == "Alive", last_contact_days_to, death_days_to),  # Use last_contact_days_to if Alive and death_days_to is NA
    event = ifelse(vital_status == "Dead", 1, 0)  # 1 for death, 0 for censored
  ) %>%
  filter(!is.na(survival_time))  # Remove rows with NA survival_time


df_cancer_H2AFY <- split(df_high_exp_survival_H2AFY, df_high_exp_survival_H2AFY$"cancer type abbreviation")

for (cancer in names(df_cancer_H2AFY)) {
  assign(paste0("H2AFY_", cancer), df_cancer_H2AFY[[cancer]], envir = .GlobalEnv)
}


output_dir <- "/Users/cmdb/Cancer_Genome_QuantBio_Project/results/Kaplan_Meier_Survival_Curves/H2AFY"

for (cancer in names(df_cancer_H2AFY)) {
  
  df_cancer <- get(paste0("H2AFY_", cancer))
  
  if (nrow(df_cancer[df_cancer$expression_group == "High", ]) > 0 &&
      nrow(df_cancer[df_cancer$expression_group == "Low", ]) > 0) {
    
    surv_object <- Surv(as.numeric(df_cancer$survival_time), df_cancer$event)
    
    fit <- survfit(surv_object ~ expression_group, data = df_cancer)
    
    # Create the Kaplan-Meier plot
    plot <- ggsurvplot(
      fit,
      data = df_cancer,
      pval = TRUE,
      title = paste("Kaplan-Meier Survival Curve for", cancer, "Cancer"),
      legend.title = "H2AFY Expression",
      legend.labs = c("Low", "High"),
      palette = c("#1b9e77", "#d95f02")  # Color-blind-friendly palette
    )
    
    output_file <- file.path(output_dir, paste0("H2AFY_", cancer, ".png"))
    ggsave(filename = output_file, plot = plot$plot)
  }
}




###################################################################
#Analysis on H1F0

df_H1F0 <- df_rnaseq %>%
  filter(gene_symbol == "H1F0")

df_long_H1F0 <- df_H1F0 %>%
  pivot_longer(
    cols = starts_with("TCGA"),
    names_to = "sample_id",
    values_to = "expression_level"
  )

df_long_H1F0 <- df_long_H1F0 %>%
  mutate(sample_id_trimmed = substr(sample_id, 1, 15))  # Extracts the first 15 characters

df_high_exp_survival_H1F0 <- df_long_H1F0 %>%
  inner_join(df_survival, by = c("sample_id_trimmed" = "sample"))

median_exp_H1F0 <- median(df_high_exp_survival_H1F0$expression_level, na.rm = TRUE)

df_high_exp_survival_H1F0 <- df_high_exp_survival_H1F0 %>%
  mutate(expression_group = ifelse(expression_level >= median_exp_H1F0, "High", "Low"))

df_high_exp_survival_H1F0 <- df_high_exp_survival_H1F0 %>%
  mutate(
    death_days_to = as.numeric(death_days_to),  # Convert to numeric
    last_contact_days_to = as.numeric(last_contact_days_to),  # Ensure numeric conversion for last_contact_days_to
    survival_time = ifelse(is.na(death_days_to) & vital_status == "Alive", last_contact_days_to, death_days_to),  # Use last_contact_days_to if Alive and death_days_to is NA
    event = ifelse(vital_status == "Dead", 1, 0)  # 1 for death, 0 for censored
  ) %>%
  filter(!is.na(survival_time))  # Remove rows with NA survival_time

df_cancer_H1F0 <- split(df_high_exp_survival_H1F0, df_high_exp_survival_H1F0$"cancer type abbreviation")



