# Load necessary libraries
library(tidyverse)   # For data manipulation and plotting
library(ggplot2)     # For advanced plotting
library(dplyr)       # For data manipulation
library(readr)       # For reading TSV files

# Load DESeq2 results dataset
deseq2_data <- read_tsv("C:/Users/HP/Desktop/R Studio/Computational Assignment/A_vs_E.deseq2.results.tsv")


# Inspect the dataset structure
head(deseq2_data)
str(deseq2_data)

# Check for missing values
missing_data_summary <- sum(is.na(deseq2_data))
cat("Number of missing values:", missing_data_summary, "\n")

# Remove rows with missing values (if any)
deseq2_data <- na.omit(deseq2_data)

# Define thresholds
logfc_threshold <- 1    # Log fold change threshold
pvalue_threshold <- 0.05 # P-value threshold

# Calculate upregulated and downregulated genes
summary_stats <- deseq2_data %>%
  mutate(Significant = padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold,
         Regulation = case_when(
           padj < pvalue_threshold & log2FoldChange > logfc_threshold ~ "Upregulated",
           padj < pvalue_threshold & log2FoldChange < -logfc_threshold ~ "Downregulated",
           TRUE ~ "Not Significant"
         )) %>%
  group_by(Regulation) %>%
  summarise(Count = n())

print(summary_stats)

# MA plot
ggplot(deseq2_data, aes(x = baseMean, y = log2FoldChange, color = padj < pvalue_threshold)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  labs(title = "MA Plot", x = "Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

# Save the plot
ggsave("ma_plot.png")



# Histogram of p-values
ggplot(deseq2_data, aes(x = padj)) +
  geom_histogram(bins = 50, fill = "green", color = "black", alpha = 0.8) +
  labs(title = "Distribution of Adjusted P-Values", x = "Adjusted P-Value", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("pvalue_histogram.png")


# Define thresholds for significance
logfc_threshold <- 1
pvalue_threshold <- 0.05

# Add significance classification
deseq2_data <- deseq2_data %>%
  mutate(Significance = case_when(
    padj < pvalue_threshold & log2FoldChange > logfc_threshold ~ "Upregulated",
    padj < pvalue_threshold & log2FoldChange < -logfc_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Volcano plot using ggplot2
ggplot(deseq2_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "green", "Downregulated" = "red", "Not Significant" = "blue")) +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Significance") +
  theme_minimal() +
  theme(legend.position = "top")

# Save the plot
ggsave("volcano_plot_ggplot.png", width = 8, height = 6)


# Prepare data for heatmap
# Select the top 20 genes based on significance
top_genes <- deseq2_data %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  select(gene_id, log2FoldChange)

# Reshape the data for ggplot (if there are multiple comparisons)
heatmap_data <- top_genes %>%
  mutate(row_id = row_number()) %>%
  pivot_wider(names_from = gene_id, values_from = log2FoldChange) %>%
  pivot_longer(-row_id, names_to = "gene", values_to = "log2FoldChange")

# Generate heatmap using ggplot2
ggplot(heatmap_data, aes(x = gene, y = factor(row_id), fill = log2FoldChange)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "yellow", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "Heatmap of Top Differentially Expressed Genes",
       x = "Gene",
       y = "Rank (Top 20 by Significance)",
       fill = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("heatmap_ggplot.png", width = 8, height = 6)

# Create lists of significant genes
significant_genes <- deseq2_data %>%
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>%
  arrange(desc(log2FoldChange))
significant_genes
head(significant_genes)
# Save to CSV
write_csv(significant_genes, "significant_genes.csv")


