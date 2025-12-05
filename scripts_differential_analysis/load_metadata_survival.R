library(DESeq2)
library(readxl)
library(dplyr)
library(tidyverse)
library(glue)
library(limma)
library(pheatmap)
library(survival)
library(survminer)

meta_data_dir <- 
  "/Users/kyu/Documents/rebecca_lab/github/foot_back_scripts/rnaseq/count_data_meta/farshidfar/RNA-seq_batches.xlsx"
meta_data_dir2 <- "/Users/kyu/Documents/rebecca_lab/github/foot_back_scripts/rnaseq/count_data_meta/farshidfar/sra_metadata_34ca34ab-9c4a-4620-bed9-74600ecb15cc.csv"
count_file_dir <- "/Users/kyu/Documents/rebecca_lab/github/foot_back_scripts/rnaseq/count_data_meta/farshidfar/gene_count_matrix.csv"

df_1a <- read_excel(meta_data_dir, sheet = "1a", skip=1)
df_1a <- as.data.frame(df_1a)
colnames(df_1a) <- make.names(colnames(df_1a))
df_1a <- df_1a %>% dplyr::select(Patient.ID, Sex, Age.at.diagnosis..years., Mutation.category, Subtype, Site.of.tumor.sampling, Primary.tumor.location, Ethnicity, AJCC.Stage)
colnames(df_1a)[3] <- "Age.at.diagnosis"
df_1c <- read_excel(meta_data_dir, sheet = "1c", skip=1)
df_1c <- as.data.frame(df_1c)
colnames(df_1c) <- make.names(colnames(df_1c))
df_1c <- df_1c %>% dplyr::select(Patient.ID, Batch)
df_1c$Patient.ID[!grepl("-", df_1c$Patient.ID)] <- paste0("YU", df_1c$Patient.ID[!grepl("-", df_1c$Patient.ID)])
df_1ac <- inner_join(df_1a, df_1c, by = "Patient.ID")


df_meta_data <- read.csv(meta_data_dir2)
df_meta_data <- df_meta_data %>% dplyr::select(SampleName, Run, Body_Site, LibraryLayout)
df_meta_data <- df_meta_data %>% 
  separate(SampleName, c("Patient.ID", "rnaseq_type"), remove = FALSE) %>%
  dplyr::filter(rnaseq_type == "rnaseq") %>%
  dplyr::select(Patient.ID,  Run, Body_Site, LibraryLayout)
df_meta_data$Patient.ID[1] <- "YUWW165"
df_meta_data <- inner_join(df_1ac, df_meta_data, by = "Patient.ID")
df_meta_data$Batch <- paste0(df_meta_data$Batch, df_meta_data$LibraryLayout)
df_meta_data <- df_meta_data %>% dplyr::select(-LibraryLayout)
patient_id <- df_meta_data$Patient.ID
df_meta_data <- df_meta_data[, -1]
rownames(df_meta_data) <- patient_id


df_meta_data$Body_Site <- gsub(" ", "_", df_meta_data$Body_Site)

run_to_name <- rownames(df_meta_data)
names(run_to_name) <- df_meta_data$Run
coldata <- df_meta_data %>% dplyr::select(Sex, 
                                   Mutation.category,    
                                   Subtype, 
                                   Site.of.tumor.sampling,    
                                   Batch, 
                                   Body_Site,
                                   AJCC.Stage)
coldata$AJCC.Stage <- as.character(coldata$AJCC.Stage)
# write.table(coldata, "/Users/kyu/Documents/rebecca_lab/acral_melanoma/external_dataset/processed_raw_AM_CM_patients_RNAseq/coldata.tsv", sep = "\t", quote = FALSE)
rownames(coldata) <- df_meta_data$Run
# write.table(coldata, "/Users/kyu/Documents/rebecca_lab/acral_melanoma/external_dataset/processed_raw_AM_CM_patients_RNAseq/coldata.tsv", sep = "\t", quote = FALSE)
cts <- read.csv(count_file_dir, row.names =  1)
shared_samples <- base::intersect(rownames(coldata), colnames(cts))
coldata <- coldata[shared_samples, ]
coldata <- mutate_if(coldata, is.character, as.factor)
AM_CM_patients_coldata <- coldata
cts <- cts[, shared_samples]


# meta_data_dir <- 
#   "/Users/kyu/Documents/rebecca_lab/github/foot_back_scripts/rnaseq/count_data_meta/farshidfar/RNA-seq_batches.xlsx"


df_1a <- read_excel(meta_data_dir, sheet = "1a", skip=1)
df_1a <- as.data.frame(df_1a)
colnames(df_1a) <- make.names(colnames(df_1a))

df_1a <- df_1a %>% dplyr::select(Patient.ID, "Overall.survival.event",
                          "Overall.survival.period.from.diagnosis.date..year.",
                          "Overall.survival.period.from.resection.date..year.",
                          "Time.between.diagnosis.and.resection..year.",
                          "Mutation.category",
                          "Subtype",
                          "Site.of.tumor.sampling",
                          "Ethnicity", "AJCC.Stage") 

colnames(df_1a) <- c("Patient.ID", "survival.event",
                     "survival.period.from.diagnosis",
                     "survival.period.from.resection",
                     "diagnosis.and.resection",
                     "Mutation.category",
                     "Subtype",
                     "Site.of.tumor.sampling",
                     "Ethnicity", "AJCC.Stage")


df_1a$survival.event <- as.numeric(df_1a$survival.event)
df_1a$survival.period.from.diagnosis <- as.numeric(df_1a$survival.period.from.diagnosis)
df_1a$survival.period.from.resection <- as.numeric(df_1a$survival.period.from.resection)
df_1a <- df_1a[rowSums(is.na(df_1a)) == 0,]
rownames(df_1a) <- df_1a$Patient.ID

df_1b <- as.data.frame(read_excel(meta_data_dir, sheet = "1b", skip=1))
