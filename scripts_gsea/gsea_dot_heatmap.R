library(glue)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)

mouse_human_genes <- read.csv("/Users/kyu/Documents/rebecca_lab/utils/HOM_MouseHumanSequence.rpt",sep="\t")
mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
# remove some columns
mouse <- mouse[,c(1,4,5)]
human <- human[,c(1,4,5)]
# merge the 2 dataset  (note that the human list is longer than the mouse one)
mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 

file_dir <- "/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/gsea/result/"
cell_lines <-
  c(
    "WM4324",
    "YUMM4-1",
    "B16_foot_vs_skin",
    "new_foot_ear",
    "all_Acral_vs_SunExposed",
    "primary_all_Acral_vs_SunExposed",
    "metastatic_all_Acral_vs_SunExposed",
    "acral_cutaneous_weiss",
    "primary_acral_cutaneous_weiss",
    "metastatic_acral_cutaneous_weiss",
    "CU_AM_CM",
    "primary_CU_AM_CM",
    "belote_sole_heel_vs_other",
    "Fibroblasts_acral_vs_arm",
    "KC_acral_vs_arm",
    "Melanocytes_acral_vs_arm",
    "primary_all_Acral_vs_SunExposed_foot_trunk"
  )

cell_lines <-
  c(
    "WM4324",
    "YUMM4-1",
    "B16_foot_vs_skin",
    "new_foot_ear",
    "all_Acral_vs_SunExposed",
    "primary_all_Acral_vs_SunExposed",
    "metastatic_all_Acral_vs_SunExposed",
    "acral_cutaneous_weiss",
    "primary_acral_cutaneous_weiss",
    "metastatic_acral_cutaneous_weiss",
    "CU_AM_CM",
    "primary_CU_AM_CM",
    "primary_all_Acral_vs_SunExposed_foot_trunk"
  )

# file_cell_lines <- c("WM4324", "YUMM1-7", "YUMM4-1", "B16_foot_vs_skin", "foot_ear", "new_foot_ear","all_Acral_vs_SunExposed", "primary_all_Acral_vs_SunExposed", "metastatic_all_Acral_vs_SunExposed", "primary_stage12_Acral_vs_SunExposed", "acral_cutaneous_weiss", "primary_acral_cutaneous_weiss", "CU_AM_CM", "primary_CU_AM_CM", 'belote_sole_heel_vs_other')

gmt_databases <- c("c2", "c5gobp", "c8", "h")
gmt_databases <- c("c2", "c5gobp",  "h")
# gmt_databases <- c("c2", "c5gobp", "h")
fdr_q_val_thresh <- 0.25
dif_results <- list()
gsea_dfs <- list()
for (cell_line in cell_lines) {
  dif_results[[glue("up_{cell_line}")]] <- c()
  dif_results[[glue("down_{cell_line}")]] <- c()
  for (gmt_database in gmt_databases) {
    gsea_pos_dir <- glue("{file_dir}/{cell_line}_{gmt_database}.GseaPreranked.*/gsea_report_for_na_pos_*.tsv")
    gsea_neg_dir <- glue("{file_dir}/{cell_line}_{gmt_database}.GseaPreranked.*/gsea_report_for_na_neg_*.tsv")
    gsea_pos_dir <- Sys.glob(gsea_pos_dir)
    gsea_neg_dir <- Sys.glob(gsea_neg_dir)
    gsea_pos_res <- read.table(gsea_pos_dir, sep = "\t", header = TRUE)
    gsea_neg_res <- read.table(gsea_neg_dir, sep = "\t", header = TRUE)
    gsea_pos_res <- gsea_pos_res[complete.cases(gsea_pos_res$NAME), ]
    gsea_neg_res <- gsea_neg_res[complete.cases(gsea_neg_res$NAME), ]
    gsea_pos_res$NES <- as.numeric(gsea_pos_res$NES)
    gsea_neg_res$NES <- as.numeric(gsea_neg_res$NES)
    gsea_pos_res$FDR.q.val <- as.numeric(gsea_pos_res$FDR.q.val)
    gsea_neg_res$FDR.q.val <- as.numeric(gsea_neg_res$FDR.q.val)
    gsea_pos_res <- gsea_pos_res[!is.na(gsea_pos_res$NAME),]
    gsea_neg_res <- gsea_neg_res[!is.na(gsea_neg_res$NAME),]
    gsea_pos_names <- gsea_pos_res[gsea_pos_res$FDR.q.val <= fdr_q_val_thresh, "NAME"]
    gsea_pos_names <- gsea_pos_names[!is.na(gsea_pos_names)]
    gsea_neg_names <- gsea_neg_res[gsea_neg_res$FDR.q.val <= fdr_q_val_thresh, "NAME"]
    dif_results[[glue("up_{cell_line}")]] <- append(dif_results[[glue("up_{cell_line}")]], gsea_pos_names)
    if ("" %in% gsea_neg_names) {
      print(cell_line)
      print(gsea_neg_names)
    }
    dif_results[[glue("down_{cell_line}")]] <- append(dif_results[[glue("down_{cell_line}")]], gsea_neg_names)
    gsea_res <- rbind(gsea_pos_res, gsea_neg_res)
    gsea_res <- gsea_res[gsea_res$NES != "---",]
    gsea_res <- gsea_res[!is.na(gsea_res$NAME),]
    if (is.null(gsea_dfs[[glue("{cell_line}")]])) {
      gsea_dfs[[glue("{cell_line}")]] <- gsea_res
    } else {
      gsea_dfs[[glue("{cell_line}")]] <- rbind(gsea_dfs[[glue("{cell_line}")]], gsea_res)
    }
  }
}
# cell_lines <- c("WM4324", "YUMM4-1", "B16_foot_vs_skin", "foot_ear", "primary_all_Acral_vs_SunExposed", "primary_acral_cutaneous_weiss", "primary_CU_AM_CM")

cell_lines <- c("WM4324", "YUMM4-1", "B16_foot_vs_skin","new_foot_ear", "primary_all_Acral_vs_SunExposed", "primary_acral_cutaneous_weiss")

genesets_1 <- unique(unname(unlist(dif_results[glue("up_{cell_lines}")])))
genesets_2 <- unique(unname(unlist(dif_results[glue("down_{cell_lines}")])))
genesets <- unique(c(genesets_1, genesets_2))
gsea_results_df <- data.frame(matrix(0, nrow = length(genesets), ncol=length(cell_lines)))
rownames(gsea_results_df) <- genesets
colnames(gsea_results_df) <- cell_lines
for (cell_line in cell_lines) {
  gsea_df <- gsea_dfs[[cell_line]]
  rownames(gsea_df) <- gsea_df$NAME
  shared_genesets <- intersect(rownames(gsea_results_df), rownames(gsea_df))
  gsea_results_df[shared_genesets, cell_line] <- gsea_df[shared_genesets, "NES"]
}

gsea_sig_df <- data.frame(matrix(1, nrow = length(genesets), ncol=length(cell_lines)))
rownames(gsea_sig_df) <- genesets
colnames(gsea_sig_df) <- cell_lines
for (cell_line in cell_lines) {
  gsea_df <- gsea_dfs[[cell_line]]
  rownames(gsea_df) <- gsea_df$NAME
  shared_genesets <- intersect(rownames(gsea_sig_df), rownames(gsea_df))
  gsea_sig_df[shared_genesets, cell_line] <- gsea_df[shared_genesets, "FDR.q.val"]
}


keywords <- c(
  "GOBP_EMBRYONIC_APPENDAGE_MORPHOGENESIS|GOBP_APPENDAGE_DEVELOPMENT|GOBP_APPENDAGE_MORPHOGENESIS|GOBP_EMBRYONIC_FORELIMB_MORPHOGENESIS|GOBP_FORELIMB_MORPHOGENESIS|GOBP_EMBRYONIC_HINDLIMB_MORPHOGENESIS|GOBP_HINDLIMB_MORPHOGENESIS",
  "GOBP_ANTERIOR_POSTERIOR_PATTERN_SPECIFICATION|GOBP_EMBRYONIC_DIGIT_MORPHOGENESIS|GOBP_PATTERN_SPECIFICATION_PROCESS|GOBP_REGIONALIZATION"
)
keywords_group <- paste(keywords, collapse ="|")

matches <- unique(grep(keywords_group, 
                        rownames(gsea_results_df), value=TRUE, ignore.case = TRUE))
gsea_results_filtered_df_selected <- gsea_results_df[matches, cell_lines]
gsea_sig_filtered_df_selected <- gsea_sig_df[matches, cell_lines]



prepare_data_for_dot_plot <- function(gsea_results_filtered_df_selected, gsea_sig_filtered_df_selected) {
  gsea_results_filtered_df_selected$geneset <- rownames(gsea_results_filtered_df_selected)
  gsea_results_filtered_df_selected_long <- reshape2::melt(gsea_results_filtered_df_selected, id.vars = "geneset")
  gsea_results_filtered_df_selected_long$variable <- factor(gsea_results_filtered_df_selected_long$variable, levels = cell_lines)
  gsea_results_filtered_df_selected_long$geneset <- factor(gsea_results_filtered_df_selected_long$geneset, levels=gsea_results_filtered_df_selected$geneset)
  gsea_results_filtered_df_selected_long[order(gsea_results_filtered_df_selected_long$variable, gsea_results_filtered_df_selected_long$geneset),]
  
  gsea_sig_filtered_df_selected$geneset <- rownames(gsea_sig_filtered_df_selected)
  gsea_sig_filtered_df_selected_long <- reshape2::melt(gsea_sig_filtered_df_selected, id.vars = "geneset")
  gsea_sig_filtered_df_selected_long$variable <- factor(gsea_sig_filtered_df_selected_long$variable, levels = cell_lines)
  gsea_sig_filtered_df_selected_long$geneset <- factor(gsea_sig_filtered_df_selected_long$geneset, levels=gsea_results_filtered_df_selected$geneset)
  gsea_sig_filtered_df_selected_long[order(gsea_sig_filtered_df_selected_long$variable, gsea_sig_filtered_df_selected_long$geneset),]
  gsea_results_filtered_df_selected_long$fdr <- gsea_sig_filtered_df_selected_long$value
  gsea_results_filtered_df_selected_long$log10_fdr <- -log10(gsea_sig_filtered_df_selected_long$value)
  return(gsea_results_filtered_df_selected_long)
}



gsea_results_filtered_df_selected_long <- prepare_data_for_dot_plot(gsea_results_filtered_df_selected, gsea_sig_filtered_df_selected)

gsea_results_filtered_df_selected_long$geneset_short <- gsea_results_filtered_df_selected_long$geneset


gsea_results_filtered_df_selected_long$geneset_short<- gsub("GOBP_", "", gsea_results_filtered_df_selected_long$geneset_short)
gsea_results_filtered_df_selected_long$geneset_short<- gsub("_", " ", gsea_results_filtered_df_selected_long$geneset_short)
gsea_results_filtered_df_selected_long$geneset_short <- factor(gsea_results_filtered_df_selected_long$geneset_short, 
                                                               levels = rev(c("EMBRYONIC HINDLIMB MORPHOGENESIS",
                                                                          "EMBRYONIC FORELIMB MORPHOGENESIS",
                                                                          "APPENDAGE DEVELOPMENT",
                                                                          "EMBRYONIC APPENDAGE MORPHOGENESIS",
                                                                          "FORELIMB MORPHOGENESIS",
                                                                          "EMBRYONIC DIGIT MORPHOGENESIS",
                                                                          "PATTERN SPECIFICATION PROCESS",
                                                                          "HINDLIMB MORPHOGENESIS",
                                                                          "APPENDAGE MORPHOGENESIS")))
p <- ggplot(
  gsea_results_filtered_df_selected_long,
  aes(
    x = variable,
    y = geneset_short,
    fill = value,
    size = fdr
  )
) +
  geom_point(shape = 21, color="transparent") + scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    x = "",
    y = "",
    fill = "NES",
    # Adjust based on the meaning of 'value'
    size = "FDR",
    # Adjust legend label for clarity
    title = "GSEA Dot Plot"
  ) + 
  scale_size_continuous(
    range = c(2, 13),
    trans = "reverse",
    limits = c(1, 0),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) + 
  guides(size = guide_legend(override.aes = list(color = "black", fill = "black", shape = 21))) + 
  theme_minimal(base_size = 14) +  # Minimal theme with a clean look
  theme(
    axis.text.y = element_text(size = 16),
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA),   # White background outside panel
    panel.grid.major = element_line(color = "gray80"),  # Light gray major grid lines
    panel.grid.minor = element_line(color = "gray90"),  # Lighter gray minor grid lines
    legend.key = element_rect(fill = "white", color = NA),  # White background for legend
    axis.text.x = element_text(angle = 45,  hjust = 1,size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
ggsave("/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/gsea/Heatmap/shared_dotplot_limb_morpho_wo_meta_all.pdf", plot = p, width = 10, height = 7.5, dpi = 300, bg = "transparent")


p <- ggplot(
  gsea_results_filtered_df_selected_long,
  aes(
    x = variable,
    y = geneset_short,
    fill = value,
    size = fdr
  )
) +
  geom_point(shape = 21, color="transparent") + scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    x = "",
    y = "",
    fill = "NES",
    # Adjust based on the meaning of 'value'
    size = "FDR",
    # Adjust legend label for clarity
    title = "GSEA Dot Plot"
  ) + 
  scale_size_continuous(
    range = c(2, 13),
    trans = "reverse",
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(1, 0)
  ) + 
  guides(size = guide_legend(override.aes = list(color = "black", fill = "black", shape = 21))) + 
  theme_minimal(base_size = 14) +  # Minimal theme with a clean look
  theme(
    axis.text.y = element_text(size = 16),
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA),   # White background outside panel
    panel.grid.major = element_line(color = "gray80"),  # Light gray major grid lines
    panel.grid.minor = element_line(color = "gray90"),  # Lighter gray minor grid lines
    legend.key = element_rect(fill = "white", color = NA),  # White background for legend
    axis.text.x = element_text(angle = 45,  hjust = 1,size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
ggsave("/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/gsea/Heatmap/shared_dotplot_limb_morpho.pdf", plot = p, width = 7, height = 5, dpi = 300, bg = "transparent")



cell_lines = c("YUMM4-1", "B16_foot_vs_skin","new_foot_ear","WM4324")
collagen_selected_pathways <-
  c(
    "REACTOME_COLLAGEN_DEGRADATION",
    "REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES",
    "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",
    "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES"
  )

gsea_results_filtered_df_selected <-
  gsea_results_df[collagen_selected_pathways, cell_lines]
gsea_sig_filtered_df_selected <-
  gsea_sig_df[collagen_selected_pathways, cell_lines]

gsea_results_filtered_df_selected_long <- prepare_data_for_dot_plot(gsea_results_filtered_df_selected, gsea_sig_filtered_df_selected)
gsea_results_filtered_df_selected_long <- gsea_results_filtered_df_selected_long[gsea_results_filtered_df_selected_long$variable %in% cell_lines,]
gsea_results_filtered_df_selected_long$geneset <- factor(gsea_results_filtered_df_selected_long$geneset, levels = (collagen_selected_pathways))
p <- ggplot(
  gsea_results_filtered_df_selected_long,
  aes(
    x = variable,
    y = geneset,
    fill = value,
    size = fdr
  )
) +
  geom_point(shape = 21, color="transparent") + scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    x = "",
    y = "",
    fill = "NES",
    # Adjust based on the meaning of 'value'
    size = "FDR",
    # Adjust legend label for clarity
    title = "GSEA Dot Plot"
  ) + 
  scale_size_continuous(
    range = c(2, 13),
    trans = "reverse",
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(1, 0)
  ) + 
  guides(size = guide_legend(override.aes = list(color = "black", fill = "black", shape = 21))) + 
  theme_minimal(base_size = 14) +  # Minimal theme with a clean look
  theme(
    axis.text.y = element_text(size = 16),
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA),   # White background outside panel
    panel.grid.major = element_line(color = "gray80"),  # Light gray major grid lines
    panel.grid.minor = element_line(color = "gray90"),  # Lighter gray minor grid lines
    legend.key = element_rect(fill = "white", color = NA),  # White background for legend
    axis.text.x = element_text(angle = 45,  hjust = 1,size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
ggsave("/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/gsea/Heatmap/shared_dotplot_collagen_selected_4.pdf", plot = p, width = 14, height = 4.5, dpi = 300, bg = "transparent")

