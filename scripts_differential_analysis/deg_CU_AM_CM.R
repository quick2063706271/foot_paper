# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(glue)

# Set output directory
out_dir <- "/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/"

# Load and preprocess data
source("~/Documents/rebecca_lab/acral_melanoma/preprocess_CU_AM_CM.R")
rownames(count_data_CU_AM_CM) <- count_data_CU_AM_CM$Gene.name
count_data_CU_AM_CM$Gene.name <- NULL

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

#' Extract and save DESeq2 results
extract_results <- function(dds, contrast_vec, output_name) {
  res <- results(dds, contrast = contrast_vec, independentFiltering = TRUE)
  res_sorted <- res[order(res$padj), ]
  res_df <- as.data.frame(res_sorted)
  res_df$Gene.name <- rownames(res_df)
  
  # Save results
  deg_dir <- glue("{out_dir}/deg/")
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res_df, glue("{deg_dir}/DEG_analysis_result_table_{output_name}.csv"))
  
  return(res_df)
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
# Analysis 1: All samples (AM vs CM)
# ==============================================================================

dds_all <- DESeqDataSetFromMatrix(
  countData = count_data_CU_AM_CM,
  colData = coldata_CU_AM_CM,
  design = ~ Subtype + Tissue.Stage
)
dds_all$Subtype <- relevel(dds_all$Subtype, ref = "CM")
dds_all <- DESeq(dds_all)
# saveRDS(dds_all, glue("{out_dir}/deseq_objects/CU_AM_CM.rds"))

generate_pca_plot(dds_all, "CU_AM_CM", intgroup = "Subtype")

res_all <- extract_results(
  dds = dds_all,
  contrast_vec = c("Subtype", "AM", "CM"),
  output_name = "CU_AM_CM"
)

generate_volcano_plot(
  dif_result = res_all,
  output_name = "CU_AM_CM",
  plot_title = "CU AM CM RNAseq\nDifferential Analysis"
)

# ==============================================================================
# Analysis 2: Primary samples only (AM vs CM)
# ==============================================================================

# Subset to primary samples
coldata_primary <- coldata_CU_AM_CM[coldata_CU_AM_CM$Tissue.Stage == "Primary", ]
count_data_primary <- count_data_CU_AM_CM[, rownames(coldata_primary)]

dds_primary <- DESeqDataSetFromMatrix(
  countData = count_data_primary,
  colData = coldata_primary,
  design = ~ Subtype
)
dds_primary$Subtype <- relevel(dds_primary$Subtype, ref = "CM")
dds_primary <- DESeq(dds_primary)
# saveRDS(dds_primary, glue("{out_dir}/deseq_objects/primary_CU_AM_CM.rds"))

generate_pca_plot(dds_primary, "primary_CU_AM_CM", intgroup = "Subtype")

res_primary <- extract_results(
  dds = dds_primary,
  contrast_vec = c("Subtype", "AM", "CM"),
  output_name = "primary_CU_AM_CM"
)

generate_volcano_plot(
  dif_result = res_primary,
  output_name = "primary_CU_AM_CM",
  plot_title = "Primary CU AM CM RNAseq\nDifferential Analysis"
)

# ==============================================================================
# Analysis 3: Metastatic samples only (AM vs CM)
# ==============================================================================

# Subset to metastatic samples
coldata_meta <- coldata_CU_AM_CM[coldata_CU_AM_CM$Tissue.Stage != "Primary", ]
count_data_meta <- count_data_CU_AM_CM[, rownames(coldata_meta)]

dds_meta <- DESeqDataSetFromMatrix(
  countData = count_data_meta,
  colData = coldata_meta,
  design = ~ Subtype
)
dds_meta$Subtype <- relevel(dds_meta$Subtype, ref = "CM")
dds_meta <- DESeq(dds_meta)
# saveRDS(dds_meta, glue("{out_dir}/deseq_objects/metastatic_CU_AM_CM.rds"))

generate_pca_plot(dds_meta, "metastatic_CU_AM_CM", intgroup = "Subtype")

res_meta <- extract_results(
  dds = dds_meta,
  contrast_vec = c("Subtype", "AM", "CM"),
  output_name = "metastatic_CU_AM_CM"
)

generate_volcano_plot(
  dif_result = res_meta,
  output_name = "metastatic_CU_AM_CM",
  plot_title = "Metastatic CU AM CM RNAseq\nDifferential Analysis"
)
