library(glue)
library(biomaRt)

setwd("~/Documents/rebecca_lab/acral_melanoma")

cell_lines <- c("new_foot_ear", "WM4324",  "YUMM4-1", "B16_foot_vs_skin",  "all_Acral_vs_SunExposed", "primary_all_Acral_vs_SunExposed", "metastatic_all_Acral_vs_SunExposed",  "acral_cutaneous_weiss", "primary_acral_cutaneous_weiss", "CU_AM_CM", "primary_CU_AM_CM")


for (i in seq_along(cell_lines)) {
  cell_line = cell_lines[i]
  dif_result <- read.csv(glue("./data/deg/DEG_analysis_result_table_{cell_line}.csv"), 
                         row.names = 1)
  dif_result <- dif_result[rowSums(is.na(dif_result)) == 0,]
  # dif_result$fcsign <- sign(dif_result$log2FoldChange)
  # dif_result$logP=-log10(dif_result$padj)
  # dif_result$metric= dif_result$logP*dif_result$log2FoldChange
  dif_result$metric <- dif_result$stat
  if (!"ID" %in% colnames(dif_result)) {
    dif_result$ID <- do.call(rbind, strsplit(rownames(dif_result), "[|]"))[,1]
  }
  y <-dif_result[,c("Gene.name", "metric")]
  # head(y)
  save_dir <- "./data/gsea/"
  dir.create(save_dir, recursive = TRUE)
  write.table(y, file=glue("{save_dir}/{cell_line}_DE_gene_names.rnk"), quote=F, sep="\t", row.names=F, col.names = F)
}




