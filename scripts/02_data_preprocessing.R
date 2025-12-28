# Script 02: Data Preprocessing
# Quality control and preparation of RNA-seq count data for DESeq2 analysis

# Load required libraries
library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

#' Load count matrix from file
#' 
#' @param count_file Path to count matrix file (CSV format)
#' @return Count matrix with genes as rows, samples as columns
load_count_matrix <- function(count_file) {
    cat("Loading count matrix from:", count_file, "\n")
    counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
    cat("Loaded", nrow(counts), "genes and", ncol(counts), "samples.\n")
    return(as.matrix(counts))
}

#' Load metadata from file
#' 
#' @param metadata_file Path to metadata file (CSV format)
#' @return Data frame with sample metadata
load_metadata <- function(metadata_file) {
    cat("Loading metadata from:", metadata_file, "\n")
    metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
    cat("Loaded metadata for", nrow(metadata), "samples.\n")
    return(metadata)
}

#' Perform quality control checks
#' 
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @param output_dir Directory to save QC plots
perform_qc <- function(counts, metadata, output_dir = "results/plots") {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    cat("Performing quality control checks...\n")
    
    # Calculate library sizes
    library_sizes <- colSums(counts)
    
    # Calculate number of genes detected per sample
    genes_detected <- colSums(counts > 0)
    
    # Create QC data frame
    qc_data <- data.frame(
        Sample = names(library_sizes),
        LibrarySize = library_sizes,
        GenesDetected = genes_detected
    )
    
    # Plot library sizes
    p1 <- ggplot(qc_data, aes(x = Sample, y = LibrarySize)) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "Library Sizes", y = "Total Counts")
    ggsave(paste0(output_dir, "/qc_library_sizes.png"), p1, width = 10, height = 6)
    
    # Plot genes detected
    p2 <- ggplot(qc_data, aes(x = Sample, y = GenesDetected)) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "Genes Detected per Sample", y = "Number of Genes")
    ggsave(paste0(output_dir, "/qc_genes_detected.png"), p2, width = 10, height = 6)
    
    # Summary statistics
    cat("\nQC Summary:\n")
    cat("Library sizes - Mean:", mean(library_sizes), "Median:", median(library_sizes), "\n")
    cat("Genes detected - Mean:", mean(genes_detected), "Median:", median(genes_detected), "\n")
    
    return(qc_data)
}

#' Prepare data for DESeq2
#' 
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @param condition_col Column name in metadata indicating condition (e.g., "condition", "group")
#' @param condition_levels Vector of condition levels (e.g., c("Control", "Disease"))
#' @return DESeqDataSet object
prepare_deseq_data <- function(counts, metadata, condition_col, condition_levels = NULL) {
    cat("Preparing DESeqDataSet...\n")
    
    # Ensure sample names match
    if (!all(colnames(counts) %in% rownames(metadata))) {
        cat("Warning: Sample names don't match. Attempting to match...\n")
        common_samples <- intersect(colnames(counts), rownames(metadata))
        counts <- counts[, common_samples]
        metadata <- metadata[common_samples, ]
    }
    
    # Ensure metadata is in same order as counts
    metadata <- metadata[colnames(counts), ]
    
    # Set condition factor
    if (!is.null(condition_levels)) {
        metadata[[condition_col]] <- factor(metadata[[condition_col]], levels = condition_levels)
    } else {
        metadata[[condition_col]] <- factor(metadata[[condition_col]])
    }
    
    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = metadata,
        design = as.formula(paste("~", condition_col))
    )
    
    # Pre-filter: remove genes with very low counts
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat("After filtering:", nrow(dds), "genes retained.\n")
    
    return(dds)
}

#' Filter samples based on QC metrics
#' 
#' @param dds DESeqDataSet object
#' @param min_library_size Minimum library size threshold
#' @param min_genes_detected Minimum number of genes detected
#' @return Filtered DESeqDataSet
filter_samples <- function(dds, min_library_size = 100000, min_genes_detected = 5000) {
    cat("Filtering samples based on QC metrics...\n")
    
    library_sizes <- colSums(counts(dds))
    genes_detected <- colSums(counts(dds) > 0)
    
    keep_samples <- library_sizes >= min_library_size & genes_detected >= min_genes_detected
    
    if (sum(!keep_samples) > 0) {
        cat("Removing", sum(!keep_samples), "samples that don't meet QC criteria.\n")
        dds <- dds[, keep_samples]
    } else {
        cat("All samples pass QC criteria.\n")
    }
    
    return(dds)
}

# Example usage:
# counts <- load_count_matrix("data/counts/GSE12345_counts.csv")
# metadata <- load_metadata("data/metadata/GSE12345_metadata.csv")
# qc_data <- perform_qc(counts, metadata)
# dds <- prepare_deseq_data(counts, metadata, condition_col = "condition", 
#                          condition_levels = c("Control", "Disease"))
# dds <- filter_samples(dds)
# save(dds, file = "data/processed/deseq_dataset.RData")

cat("\n=== Data Preprocessing Script Complete ===\n")
cat("Next steps:\n")
cat("1. Load your count matrix and metadata\n")
cat("2. Run QC checks\n")
cat("3. Prepare DESeqDataSet\n")
cat("4. Proceed to 03_deg_analysis.R\n")
