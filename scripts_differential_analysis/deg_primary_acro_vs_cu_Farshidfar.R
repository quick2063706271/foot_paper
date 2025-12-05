# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(glue)
library(dplyr)

# Load metadata and survival data
source("/Users/kyu/Documents/rebecca_lab/github/foot_back_scripts/rnaseq/differential_analysis/load_metadata_survival.R")

# Set output directory
out_dir <- "/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/"

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Generate and save PCA plot
generate_pca_plot <- function(dds, output_name, intgroup = "Subtype") {
  vst_data <- vst(dds)
  
  pca_dir <- glue("{out_dir}/pca/")
  dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(glue("{pca_dir}/{output_name}.pdf"))
  print(plotPCA(vst_data, intgroup = intgroup))
  dev.off()
  
  return(vst_data)
}

#' Extract and save DESeq2 results (using result name)
extract_results_by_name <- function(dds, result_name, output_name) {
  res <- results(dds, name = result_name, independentFiltering = TRUE)
  res_df <- as.data.frame(res)
  
  # Extract gene names from rownames (format: ID|GeneName)
  res_df$Gene.name <- sapply(strsplit(rownames(res_df), "[|]"), `[`, 2)
  
  # Save full results
  deg_dir <- glue("{out_dir}/deg/")
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res_df, glue("{deg_dir}/DEG_analysis_result_table_{output_name}.csv"))
  
  # Filter complete cases and sort by adjusted p-value
  res_filtered <- res_df[complete.cases(res_df$padj), ]
  res_sorted <- res_filtered[order(res_filtered$padj), ]
  
  return(res_sorted)
}

#' Generate volcano plot
generate_volcano_plot <- function(dif_result, output_name, plot_title, 
                                   padj_threshold = 0.05, log2fc_threshold = 1) {
  # Calculate -log10(padj)
  dif_result$log10Padj <- -log10(dif_result$padj)
  
  # Identify significant genes
  sig_genes <- dif_result[which(dif_result$log10Padj > -log10(padj_threshold)), ]
  sig_upregulated <- sig_genes[which(sig_genes$log2FoldChange > log2fc_threshold), ]
  sig_downregulated <- sig_genes[which(sig_genes$log2FoldChange < -log2fc_threshold), ]
  
  # Count DE genes
  num_upregulated <- nrow(sig_upregulated)
  num_downregulated <- nrow(sig_downregulated)
  
  # Combine for labeling
  sig_for_labeling <- rbind(sig_upregulated, sig_downregulated)
  
  # Create plot
  volcano_dir <- glue("{out_dir}/volcano_plot/")
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(glue("{volcano_dir}/{output_name}.pdf"))
  
  p <- ggplot(dif_result, aes(x = log2FoldChange, y = log10Padj)) + 
    # All genes (gray)
    geom_point(colour = "grey", alpha = 0.6) + 
    # Upregulated genes
    geom_point(data = sig_upregulated, 
               aes(colour = "up-regulated")) +
    # Downregulated genes
    geom_point(data = sig_downregulated, 
               aes(colour = "down-regulated")) +
    # Threshold lines
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
               col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_threshold), 
               col = "black", linetype = "dashed") + 
    # Gene labels
    geom_text_repel(data = sig_for_labeling, 
                    aes(label = Gene.name), 
                    size = 2, 
                    max.overlaps = 20) +
    # Labels and theme
    ggtitle(plot_title) + 
    xlab("log2 Fold Change") +
    ylab("-log10(Adjusted P-value)") +
    theme_bw() + 
    theme(
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 14), 
      text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(values = c("up-regulated" = "#bb0c00", 
                                   "down-regulated" = "#00AFBB")) +
    # Gene count annotations
    annotate("text", x = Inf, y = Inf,
             label = paste("Upregulated genes:", num_upregulated),
             hjust = 1.1, vjust = 2, size = 3.5, color = "#bb0c00") +
    annotate("text", x = -Inf, y = Inf,
             label = paste("Downregulated genes:", num_downregulated),
             hjust = -0.1, vjust = 2, size = 3.5, color = "#00AFBB") + 
    ylim(c(0, max(dif_result$log10Padj, na.rm = TRUE) + 2))
  
  print(p)
  dev.off()
}

# ==============================================================================
# DESeq2 Analysis: Primary Skin Samples - Acral vs SunExposed
# ==============================================================================

# Filter to primary skin samples only
coldata_primary_skin <- coldata %>% 
  dplyr::filter(Site.of.tumor.sampling == "Primary") %>% 
  dplyr::filter(Body_Site == "skin")

cts_primary_skin <- cts[, rownames(coldata_primary_skin)]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = cts_primary_skin,
  colData = coldata_primary_skin,
  design = ~ Batch + Sex + Subtype + AJCC.Stage
)

# Set reference level for Subtype
dds$Subtype <- relevel(dds$Subtype, ref = "SunExposed")

# Run DESeq2 analysis
dds <- DESeq(dds)

# ==============================================================================
# PCA Plot
# ==============================================================================

generate_pca_plot(dds, "primary_all_Acral_vs_SunExposed", intgroup = "Subtype")

# ==============================================================================
# Differential Expression Analysis
# ==============================================================================

res_sorted <- extract_results_by_name(
  dds = dds,
  result_name = "Subtype_Acral_vs_SunExposed",
  output_name = "primary_all_Acral_vs_SunExposed"
)

# ==============================================================================
# Volcano Plot
# ==============================================================================

generate_volcano_plot(
  dif_result = res_sorted,
  output_name = "primary_acro_vs_cutaneous_volcano_plot",
  plot_title = "Primary Acral vs Cutaneous RNAseq\nDifferential Analysis"
)