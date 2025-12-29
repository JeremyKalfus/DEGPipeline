# Main Workflow Script for Parkinsons DEG Analysis Pipeline
# This script orchestrates the entire biomarker discovery pipeline

library(DESeq2)
library(dplyr)
library(ggplot2)

set.seed(42)

cat("========================================\n")
cat("DEG Analysis Pipeline\n")
cat("========================================\n\n")

cat("Loading analysis scripts...\n")
source("scripts/01_data_download.R")
source("scripts/02_data_preprocessing.R")
source("scripts/03_deg_analysis.R")
source("scripts/04_visualization.R")
source("scripts/05_pathway_enrichment.R")
source("scripts/06_biomarker_identification.R")
cat("Scripts loaded.\n\n")

# ============================================================================
# STEP 1: Data Download
# ============================================================================
cat("STEP 1: Data Download\n")
cat("---------------------\n")
cat("Available commands:\n\n")
cat("1. Download a specific GEO dataset:\n")
cat("   gse_id <- \"GSE12345\"  # Replace with your GSE ID\n")
cat("   gse <- download_geo_dataset(gse_id)\n")
cat("   gse <- download_geo_dataset(gse_id, dest_dir = \"data/raw\")\n\n")
cat("3. Extract count matrix from downloaded dataset:\n")
cat("   extract_count_matrix(gse, \"data/counts/GSE12345_counts.csv\")\n")
cat("   extract_count_matrix(gse, paste0(\"data/counts/\", gse_id, \"_counts.csv\"))\n\n")
cat("4. Extract metadata from downloaded dataset:\n")
cat("   metadata <- extract_metadata(gse, \"data/metadata/GSE12345_metadata.csv\")\n")
cat("   metadata <- extract_metadata(gse, paste0(\"data/metadata/\", gse_id, \"_metadata.csv\"))\n\n")
cat("Complete example workflow:\n")
cat("   gse_id <- \"GSE12345\"\n")
cat("   gse <- download_geo_dataset(gse_id)\n")
cat("   extract_count_matrix(gse, paste0(\"data/counts/\", gse_id, \"_counts.csv\"))\n")
cat("   metadata <- extract_metadata(gse, paste0(\"data/metadata/\", gse_id, \"_metadata.csv\"))\n\n")
