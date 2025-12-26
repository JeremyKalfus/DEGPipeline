library(GEOquery)
library(ArrayExpress)
library(dplyr)

library(readr)
#' @param gse_id GEO Series ID (e.g., "GSE12345")
#' @param dest_dir Destination directory for downloaded data
#' @return GEO object containing expression data and metadata
download_geo_dataset <- function(gse_id, dest_dir = "data/raw") {
    cat("Downloading GEO dataset:", gse_id, "\n")
    
    gse <- getGEO(GEO = gse_id, destdir = dest_dir, getGPL = TRUE)
    
    if (is.list(gse)) {
        gse <- gse[[1]]
    }
    
    cat("Download complete. Dataset contains ", ncol(exprs(gse)), " samples.\n")
    return(gse)
}

#' @param gse GEO object
#' @param output_file Output file path for count matrix
extract_count_matrix <- function(gse, output_file) {
    expr_matrix <- exprs(gse)
    
    write.csv(expr_matrix, file = output_file, row.names = TRUE)
    cat("Count matrix saved to:", output_file, "\n")
}

#' 
#' @param gse GEO object
#' @param output_file Output file path for metadata
extract_metadata <- function(gse, output_file) {
    metadata <- pData(gse)
    
    write.csv(metadata, file = output_file, row.names = TRUE)
    cat("Metadata saved to:", output_file, "\n")
    return(metadata)
}

#' Search for RNA-seq datasets in GEO
#' 
#' @param search_terms Search terms (default: "RNA-seq")
#' @return Information about searching GEO
search_datasets <- function(search_terms = "RNA-seq") {
    cat("Searching GEO for:", search_terms, "\n")
    cat("Note: Manual search recommended at https://www.ncbi.nlm.nih.gov/geo/\n")
    cat("Look for datasets with:\n")
    cat("  - RNA-seq data type\n")
    cat("  - Appropriate sample types for your analysis\n")
    cat("  - Disease vs control comparisons\n")
    cat("\nExample search queries:\n")
    cat("  - 'Parkinson[Title] AND blood[Title] AND RNA-seq[Study Type]'\n")
    cat("  - 'disease[Title] AND RNA-seq[Study Type]'\n")
    cat("  - Use GEO advanced search for more options\n")
}

cat("\n=== Data Download Script Complete ===\n")