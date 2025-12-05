# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(glue)
library(dplyr)

# Set output directory
out_dir <- "~/Downloads/"

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Extract and save DESeq2 results (using contrast)
extract_results_by_contrast <- function(dds, contrast_vec, output_name) {
  res <- results(dds, contrast = contrast_vec, independentFiltering = TRUE)
  res_sorted <- res[order(res$padj), ]
  res_df <- as.data.frame(res_sorted)
  
  # Extract gene names and IDs from rownames (format: ID|GeneName)
  res_df$Gene.name <- sapply(strsplit(rownames(res_df), "[|]"), `[`, 2)
  res_df$ID <- sapply(strsplit(rownames(res_df), "[|]"), `[`, 1)
  
  # Save full results (including NA padj values)
  deg_dir <- glue("{out_dir}/deg/")
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res_df, glue("{deg_dir}/DEG_analysis_result_table_{output_name}.csv"))
  
  # Filter complete cases for downstream analysis
  res_filtered <- res_df[complete.cases(res_df$padj), ]
  
  return(res_filtered)
}

#' Generate volcano plot with custom ggrepel settings
generate_volcano_plot_custom <- function(dif_result, output_name, plot_title, 
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
                    size = 4, 
                    max.overlaps = 20,
                    force = 1,
                    min.segment.length = 0,
                    nudge_y = 1) +
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
# Load and Preprocess Data
# ==============================================================================

# Load count data
count_data_file <- "~/Documents/rebecca_lab/acral_melanoma/data/raw_count/YUMM4-1_raw_counts.csv"
count_data <- read.csv(count_data_file)

# Create rownames with gene IDs and names
rownames(count_data) <- paste(count_data$ID, count_data$Gene.name, sep = "|")
count_data <- count_data %>% dplyr::select(-ID, -Gene.name)

# Create sample metadata
coldata <- data.frame(
  cond = c("foot", "foot", "foot", "back", "back", "back"),
  row.names = colnames(count_data)
)

# ==============================================================================
# DESeq2 Analysis: YUMM4-1 Foot vs Back
# ==============================================================================

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = coldata,
  design = ~ cond
)

# Set reference level
dds$cond <- relevel(dds$cond, ref = "back")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Optional: Save/load DESeq2 object
# saveRDS(dds, glue("{out_dir}/deseq_objects/YUMM4.1.rds"))
# dds <- readRDS(glue("{out_dir}/deseq_objects/YUMM4.1.rds"))

# ==============================================================================
# Differential Expression Analysis
# ==============================================================================

res_sorted <- extract_results_by_contrast(
  dds = dds,
  contrast_vec = c("cond", "foot", "back"),
  output_name = "YUMM4-1"
)

# ==============================================================================
# Volcano Plot
# ==============================================================================

generate_volcano_plot_custom(
  dif_result = res_sorted,
  output_name = "YUMM4_1",
  plot_title = "YUMM4.1 Foot vs Flank",
  log2fc_threshold = 0.5
)
