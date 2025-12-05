library(DESeq2)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggfortify)
library(RColorBrewer)
library(viridis)
library(glue)
library(biomaRt)
library(GenomicFeatures)
library(DESeq2)
library(tidyverse)
library("readxl")

setwd("~/Documents/rebecca_lab/public_data/weiss_2022_nat_dataset/rnaseq/")
sampleinfo <- read.csv("AM_vs_CM_samples.csv", header= TRUE, stringsAsFactors=F)
sampleinfo <- sampleinfo[2:115,1:6] #simplify description and remove patient 4550 who may have either acral or cutaneous
sampleinfo$Subtype[sampleinfo$Subtype == "acral"] <- "Acral"
sampleinfo$Subtype[sampleinfo$Subtype == "cutaneous"] <- "Cutaneous"
sampleinfo$Specimen.Type[sampleinfo$Specimen.Type == "Local recurrence"] <- "Local Recurrence"
sampleinfo$group <- paste(sampleinfo$Subtype, sampleinfo$Specimen.Type, sep = " ") 



AM_vs_CM_counts <- read.csv("Melanoma_AM_vs_CM.genes.ExpectedCounts.csv",
                            header=TRUE,
                            comment = "#")
dim(AM_vs_CM_counts)
AM_vs_CM_counts$MELA_4550 <- NULL # remove this sample from the dataset
dim(AM_vs_CM_counts)

countdata <- AM_vs_CM_counts %>%
  column_to_rownames("Ensembl.Gene.ID") # turn the geneid column into rownames
countdata <- data.matrix(countdata[,8:121,], rownames.force = NA)
head(countdata)
dim(countdata)

#filter out genes with few to no counts
keep <- rowSums(countdata, na.rm = TRUE) > 20
countdata <- countdata[keep,]
dim(countdata)

#plot library size by sample name
librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=20e6, lty=2)

#need to round the countdata to integer values to make it compatible with vst transformation
countdata <- round(countdata, digits = 0)
dim(countdata)
head(countdata)

#vst normalization instead of rlog due to having >30 samples
vst_count <- vst(countdata, blind = FALSE)
head(vst_count)

pcDat <- prcomp(t(vst_count))

#specify order of groups on legend
sampleinfo$group <- factor(sampleinfo$group, 
                           levels = c("Acral Primary", 
                                      "Acral Local Recurrence", 
                                      "Acral ITM", 
                                      "Acral Regional Lymph Nodes", 
                                      "Acral Distant Lymph Nodes",
                                      "Cutaneous Primary", 
                                      "Cutaneous Local Recurrence", 
                                      "Cutaneous ITM", 
                                      "Cutaneous Regional Lymph Nodes", 
                                      "Cutaneous Distant Lymph Nodes"))
autoplot(pcDat,
         data = sampleinfo,
         colour="group",
         size=3,
         x = 1,
         y = 2,
         scale = TRUE) +
  theme_bw() +
  scale_color_brewer(palette = "RdYlBu", direction = -1) +
  ggtitle("PCA of all samples")


# first lets check that our rows and columns match
all(sampleinfo$donorLabel == colnames(countdata))
# create the design formula
# compare gene expression between all acral and cutaneous melanoma samples
design <- as.formula(~ Subtype) 
# compare gene expression between all acral and cutaneous melanoma samples and taking into account differences between sample type (ie. primary tumor, lymph node metastasis, distant metastsis, etc)
design_complex <- as.formula(~ Specimen.Type + Subtype) 


# Primary samples
sampleinfo_primary <- sampleinfo[sampleinfo$Specimen.Type == "Primary",]
countdata_primary <- countdata[, sampleinfo_primary$donorLabel]
ddsObj <- DESeqDataSetFromMatrix(countData = countdata_primary,
                                 colData = sampleinfo_primary,
                                 design = ~Subtype)

# Run DESeq
ddsObj <- DESeq(ddsObj, betaPrior = TRUE)
ddsObj$Subtype <- relevel(ddsObj$Subtype, ref = "Cutaneous")
resAvC <- results(ddsObj, contrast = c("Subtype", "Acral", "Cutaneous"))





mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))


resAvC_list = list(resAvC = as.data.frame(resAvC),
                   resAvC_complex = as.data.frame(resAvC_complex)
)

resAvC_list_named = lapply(resAvC_list, function(j) {
  j %>% 
    rownames_to_column(., "Ensembl") %>%
    merge(x = .,
          y = getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","description", "wikigene_name"),
                    values = .$Ensembl, 
                    mart = mart),
          by.x = "Ensembl",
          by.y = "ensembl_gene_id") %>%
    .[.$wikigene_name != "",]
})

resAvC_data <- resAvC_list_named[["resAvC"]]
resAvC_data_new <- resAvC_data
colnames(resAvC_data_new)[1] <- "ID"
colnames(resAvC_data_new)[ncol(resAvC_data_new)] <- "Gene.name"
write.csv(resAvC_data_new, "~/Documents/rebecca_lab/acral_melanoma/data/deg/DEG_analysis_result_table_primary_acral_cutaneous_weiss.csv")

# saveRDS(ddsObj, "/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/deseq_objects/primary_acral_cutaneous_weiss.rds")
ddsObj <- readRDS("/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/deseq_objects/primary_acral_cutaneous_weiss.rds")
vst <- vst(ddsObj)
# assay(vst) <- limma::removeBatchEffect(assay(vst))
pdf("/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/pca/primary_acral_cutaneous_weiss.pdf")
plotPCA(vst, intgroup = "Subtype")
dev.off()




# Metastatic samples
sampleinfo_metastatic <- sampleinfo[sampleinfo$Specimen.Type != "Primary",]
countdata_metastatic <- countdata[, sampleinfo_metastatic$donorLabel]
ddsObj <- DESeqDataSetFromMatrix(countData = countdata_metastatic,
                                 colData = sampleinfo_metastatic,
                                 design = design)

# Run DESeq
ddsObj <- DESeq(ddsObj, betaPrior = TRUE)
ddsObj$Subtype <- relevel(ddsObj$Subtype, ref = "Cutaneous")
resAvC <- results(ddsObj, alpha=0.05, contrast = c("Subtype", "Acral", "Cutaneous"))
library(biomaRt)
library(GenomicFeatures)
library(DESeq2)
library(tidyverse)
# library(xlsx)
library("readxl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))


resAvC_list = list(resAvC = as.data.frame(resAvC),
                   resAvC_complex = as.data.frame(resAvC_complex)
)

resAvC_list_named = lapply(resAvC_list, function(j) {
  j %>% 
    rownames_to_column(., "Ensembl") %>%
    merge(x = .,
          y = getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","description", "wikigene_name"),
                    values = .$Ensembl, 
                    mart = mart),
          by.x = "Ensembl",
          by.y = "ensembl_gene_id") %>%
    .[.$wikigene_name != "",]
})

resAvC_data <- resAvC_list_named[["resAvC"]]
resAvC_data_new <- resAvC_data
colnames(resAvC_data_new)[1] <- "ID"
colnames(resAvC_data_new)[ncol(resAvC_data_new)] <- "Gene.name"
write.csv(resAvC_data_new, "~/Documents/rebecca_lab/acral_melanoma/data/deg/DEG_analysis_result_table_metastatic_acral_cutaneous_weiss.csv")




vsd <- vst(ddsObj, blind = F)
pdf(glue("/Users/kyu/Documents/rebecca_lab/acral_melanoma/data/pca/primary_acral_cutaneous_weiss.pdf"), width = 5, height = 5)
print(plotPCA(vsd, intgroup=c("Subtype")))
dev.off()



##For all samples
# rebuild a clean DDS object
ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleinfo,
                                 design = design)

# Run DESeq
ddsObj <- DESeq(ddsObj, betaPrior = TRUE)
ddsObj$Subtype <- relevel(ddsObj$Subtype, ref = "Cutaneous")
##########################################################################################

###Do the same but with the complex design
ddsObj_complex <- DESeqDataSetFromMatrix(countData = countdata,
                                         colData = sampleinfo,
                                         design = design_complex)

# Run DESeq
ddsObj_complex <- DESeq(ddsObj_complex, betaPrior = TRUE)
ddsObj_complex$Subtype <- relevel(ddsObj_complex$Subtype, ref = "Cutaneous")

