
library(Seurat)
library(data.table)
library(arrow)
library(dplyr)
library(glue)



meta_file = "/dcs05/hongkai/data/kyu/rebecca_lab/public_data/Wiedemann_2023_Cell_Reports/metadata.csv"
count_file = "~/raw_counts.csv"
# count_file = "/Users/kyu/Documents/rebecca_lab/public_data/Wiedemann_2023_Cell_Reports/raw_counts.feather"
counts <- as.data.frame(fread(count_file, sep=" "))
rownames(counts) <- counts$V1

counts <- counts %>% select(-V1)  # Remove by column name

colnames(counts) <- gsub(".", "-", colnames(counts), fixed = TRUE)

# counts <- as.data.frame(counts)

meta <- read.csv(meta_file, row.names = 1, sep = " ")  # If metadata is a CSV file
meta <- meta[colnames(counts),]
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta)

saveRDS(seurat_obj, "/dcs05/hongkai/data/kyu/rebecca_lab/public_data/Wiedemann_2023_Cell_Reports/seurat.rds")

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


seurat_subset <- subset(seurat_obj, subset = celltype == "Melanocytes")

de_results <- FindMarkers(seurat_subset, ident.1 = "acral", ident.2 = "arm", group.by = "bodysite", test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
head(de_results)

