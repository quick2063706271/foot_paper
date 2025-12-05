# Load required libraries
library(readxl)
library(matrixStats)

# ==============================================================================
# Load and Process Mouse-Human Gene Mapping
# ==============================================================================

# Load mouse-human gene mapping data
mouse_human_genes <- read.csv("~/Documents/rebecca_lab/github/foot_paper/scripts_vc_ratio/HOM_MouseHumanSequence.rpt", sep = "\t")

# Split by organism
mouse <- split.data.frame(mouse_human_genes, mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes, mouse_human_genes$Common.Organism.Name)[[1]]

# Keep only relevant columns
mouse <- mouse[, c(1, 4, 5)]
human <- human[, c(1, 4, 5)]

# Merge datasets (note: human list is longer than mouse)
mh_data <- merge.data.frame(mouse, human, by = "DB.Class.Key", all.y = TRUE)

# ==============================================================================
# Define Gene Signatures
# ==============================================================================

# VMEL (Vascular Melanoma) signature - top 8 genes
vmel_genes_human_top_8 <- c(
  "ID3",
  "NTRK2",
  "ID2",
  "MEG3",
  "RAB3B",
  "IGDCC4",
  "MIA",
  "PDLIM4"
)

# CMEL (Cutaneous Melanoma) signature - top 9 genes
cmel_genes_human_top_8 <- c(
  "AKAP12",
  "SLC45A2",
  "HPGD",
  "MCOLN3",
  "RGL1",
  "SEMA5A",
  "ACP5",
  "APCDD1",
  "GALNT18"
)

# Convert to mouse gene symbols
vmel_genes_mouse_top_8 <- mh_data[mh_data$Symbol.y %in% vmel_genes_human_top_8, "Symbol.x"]
cmel_genes_mouse_top_8 <- mh_data[mh_data$Symbol.y %in% cmel_genes_human_top_8, "Symbol.x"]

# ==============================================================================
# VMEL/CMEL Ratio Calculation Function
# ==============================================================================

#' Calculate VMEL/CMEL ratio from normalized gene expression data
#'
#' @param normalized_counts Matrix of normalized expression values (e.g., VST from DESeq2)
#' @param specie Character. Species of the data: "human" or "mouse"
#' @param gene_sets Character. Gene set to use: "top_8" (default)
#' @param method Character. Aggregation method: "mean" (default) or "prod"
#'
#' @return Numeric vector of VMEL/CMEL ratios for each sample
calculate_vmel_cmel_ratio <- function(normalized_counts, 
                                       specie = "human", 
                                       gene_sets = "top_8", 
                                       method = "mean") {
  
  # Select gene sets based on species
  if (gene_sets == "top_8") {
    if (specie == "human") {
      Vmel_genes <- vmel_genes_human_top_8
      Cmel_genes <- cmel_genes_human_top_8
      
      # Handle missing GALNT18 gene
      if (!"GALNT18" %in% rownames(normalized_counts)) {
        Cmel_genes[Cmel_genes == "GALNT18"] <- "GALNT14"
        message("Note: GALNT18 not found in data, using GALNT14 as substitute")
      }
    } else if (specie == "mouse") {
      Vmel_genes <- vmel_genes_mouse_top_8
      Cmel_genes <- cmel_genes_mouse_top_8
    } else {
      stop("Species must be 'human' or 'mouse'")
    }
  }
  
  # Print genes being used
  cat("VMEL genes:\n")
  cat(paste(Vmel_genes, collapse = ", "), "\n\n")
  cat("CMEL genes:\n")
  cat(paste(Cmel_genes, collapse = ", "), "\n\n")
  
  # Ensure matrix format
  normalized_counts <- as.matrix(normalized_counts)
  
  # Filter to genes present in data
  vmel_present <- rownames(normalized_counts) %in% Vmel_genes
  cmel_present <- rownames(normalized_counts) %in% Cmel_genes
  
  # Calculate scores based on method
  if (method == "mean") {
    v_score <- colMeans(normalized_counts[vmel_present, , drop = FALSE])
    c_score <- colMeans(normalized_counts[cmel_present, , drop = FALSE])
  } else if (method == "prod") {
    v_score <- matrixStats::colProds(normalized_counts[vmel_present, , drop = FALSE])
    c_score <- matrixStats::colProds(normalized_counts[cmel_present, , drop = FALSE])
  } else {
    stop("Method must be 'mean' or 'prod'")
  }
  
  # Calculate ratio
  vc_ratio <- v_score / c_score
  
  return(vc_ratio)
}

