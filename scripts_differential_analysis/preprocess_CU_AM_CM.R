library(readxl)

input_file <- "/Users/kyu/Documents/rebecca_lab/acral_melanoma/external_dataset/CU_AM_CM/KASEY_PT_RNAseq.xlsx"

excel_df <- as.data.frame(read_excel(input_file, sheet = "counts_table_combined-KC format",  skip = 10))

count_data_CU_AM_CM <- excel_df[, (colnames(excel_df)[16:ncol(excel_df)])]
count_data_CU_AM_CM$Gene.name <- excel_df$GENE


excel_df <- as.data.frame(read_excel(input_file, sheet = "counts_table_combined-KC format",  col_names = FALSE))
coldata_CU_AM_CM <- excel_df[c(1:9,11), c(15:ncol(excel_df))]
coldata_CU_AM_CM <- t(coldata_CU_AM_CM)
sample_names <- as.character(coldata_CU_AM_CM[, 10])
sample_names <- sample_names[2: length(sample_names)]
colnames_CU_AM_CM <- as.character(coldata_CU_AM_CM[1,])
colnames_CU_AM_CM <- colnames_CU_AM_CM[1: length(colnames_CU_AM_CM)-1]
coldata_CU_AM_CM <- coldata_CU_AM_CM[2: nrow(coldata_CU_AM_CM), 1: ncol(coldata_CU_AM_CM)-1]
colnames(coldata_CU_AM_CM) <- colnames_CU_AM_CM
rownames(coldata_CU_AM_CM) <- sample_names
coldata_CU_AM_CM[is.na(coldata_CU_AM_CM)] <- "<na>"
colnames(coldata_CU_AM_CM) <- make.names(colnames(coldata_CU_AM_CM))
coldata_CU_AM_CM <- as.data.frame(coldata_CU_AM_CM)
