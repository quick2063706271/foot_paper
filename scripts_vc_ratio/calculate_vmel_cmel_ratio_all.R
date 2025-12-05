source("~/Documents/rebecca_lab/acral_melanoma/scripts/calculate_vmel_cmel_ratio.R")
library(glue)
library(ggpubr)
library(ggplot2)
library(ComplexHeatmap)
library(DESeq2)
library(biomaRt)

deseq_obj_dir = "/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/deseq_objects/"
cell_lines <- c(
  "new_foot_ear",
  "B16",
  "YUMM4.1",
  "CU_AM_CM",
  "primary_CU_AM_CM",
  "metastatic_CU_AM_CM",
  "WM4324",
  "acral_cutaneous_weiss",
  "all_Acral_vs_SunExposed",
  "metastatic_acral_cutaneous_weiss",
  "metastatic_all_Acral_vs_SunExposed",
  "primary_acral_cutaneous_weiss",
  "primary_all_Acral_vs_SunExposed"
)

species <- c(
  "mouse",
  "mouse",
  "mouse",
  "human",
  "human",
  "human",
  "human",
  "human",
  "human",
  "human",
  "human",
  "human",
  "human"
)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # For human genes

get_cnts_ready <- function(normalized_counts) {
  if (sum(grepl("|", rownames(normalized_counts), fixed = T))>1) {
    rownames(normalized_counts) <- make.unique(do.call(rbind, strsplit(rownames(normalized_counts), "[|]"))[,2])
  }
  if (sum(grepl("ENSG", rownames(normalized_counts), fixed = T))==nrow(normalized_counts)) {
    gene_ids <- rownames(normalized_counts)
    gene_mapping <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = ensembl
    )
    
    # Merge the gene names with your gene count table
    # Ensure that the order matches
    gene_mapping <- gene_mapping[match(gene_ids, gene_mapping$ensembl_gene_id), ]
    
    # Replace rownames with gene names, if available
    new_row_names <- ifelse(
      is.na(gene_mapping$hgnc_symbol) | gene_mapping$hgnc_symbol == "",
      gene_ids,  # Keep gene ID if gene name is missing
      gene_mapping$hgnc_symbol
    )
    
    # Make the new row names unique before assigning
    new_row_names <- make.unique(new_row_names)
    
    # Assign the unique row names to the gene count table
    rownames(normalized_counts) <- new_row_names
  }
  return(normalized_counts)
}
subtype_colnames <- c("Subtype", "cond", "subtype", "Cond", "condition", "sample_type")

foot_subtypes <- c("Acral", "Foot", "foot", "AM", "acral")
plots <- list()

out_dir <- "/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/boxplots/8_vst_clean"
for (i in seq_along(cell_lines)) {
  
  cell_line <- cell_lines[i]
  plots[[cell_line]] <- list()
  specie <- species[i]
  deseq_object_file <- glue("{deseq_obj_dir}/{cell_line}.rds")
  dds <- readRDS(deseq_object_file)
  
  normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
  normalized_counts <- get_cnts_ready(normalized_counts)
  
  print(normalized_counts)
  vds<-vst(dds, blind = FALSE)
  normalized_vds_counts<-as.data.frame(assay(vds))
  normalized_vds_counts <- get_cnts_ready(normalized_vds_counts)
  
  nt <- normTransform(dds)
  nt_counts <- as.data.frame(assay(nt))
  nt_counts <- get_cnts_ready(nt_counts)
  
  for (gene_set in c("top_8")) {
    vc_ratio <- calculate_vmel_cmel_ratio(normalized_vds_counts, specie = specie, gene_sets = gene_set)
    subtype_colname <- subtype_colnames[subtype_colnames %in% colnames(colData(dds))]
    if (length(subtype_colname) > 1) {
      subtype_colname=subtype_colname[1]
    }
    dat <- data.frame(reads=as.numeric(vc_ratio), Subtype=as.data.frame(colData(dds))[[subtype_colname]])
    foot_subtype = foot_subtypes[foot_subtypes %in% dat$Subtype]
    back_subtype = as.character(unique(dat$Subtype)[unique(dat$Subtype)!=foot_subtype])
    dat$Subtype <- relevel(dat$Subtype, ref = back_subtype)
    p <- ggboxplot(dat, x = "Subtype", y = "reads", fill = "Subtype", palette = c("#00AFBB", "#E7B800"), add = "jitter", legend = "right") +
      ggtitle(glue("v-mel:c-mel score {gene_set}\n {cell_line}"))  + ylab("v-mel:c-mel score")+ stat_compare_means(method = "t.test", paired = FALSE, label.x = 1.5, label.y = max(dat$reads))
    plots[[cell_line]][[gene_set]][["vst"]] <- p
    pdf(glue("{out_dir}/vmel_cmel_ratio_{cell_line}_vds.pdf"), width = 5, height = 5)
    print(p)
    dev.off()
    write.csv(dat, glue("{out_dir}/vmel_cmel_ratio_{cell_line}_vds.csv"))
  }
}
