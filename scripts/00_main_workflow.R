# Main Workflow Script for Parkinsons DEG Analysis Pipeline
# This script orchestrates the entire biomarker discovery pipeline

# Load required libraries
library(DESeq2)
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

cat("========================================\n")
cat("DEG Analysis Pipeline\n")
cat("========================================\n\n")

# Source all analysis scripts
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

# ============================================================================
# STEP 2: Data Preprocessing
# ============================================================================
cat("STEP 2: Data Preprocessing\n")
cat("--------------------------\n")
cat("Load your count matrix and metadata, then run preprocessing:\n\n")

# Example (uncomment and modify):
# counts <- load_count_matrix("data/counts/GSE12345_counts.csv")
# metadata <- load_metadata("data/metadata/GSE12345_metadata.csv")
# qc_data <- perform_qc(counts, metadata)
# dds <- prepare_deseq_data(counts, metadata, condition_col = "condition", 
#                           condition_levels = c("Control", "PD"))
# dds <- filter_samples(dds)
# save(dds, file = "data/processed/deseq_dataset.RData")

# ============================================================================
# STEP 3: DEG Analysis
# ============================================================================
cat("STEP 3: DEG Analysis\n")
cat("--------------------\n")
cat("Run DESeq2 analysis and extract biomarkers:\n\n")

# Example (uncomment and modify):
# load("data/processed/deseq_dataset.RData")
# analysis <- run_deseq_analysis(dds, contrast = c("condition", "PD", "Control"))
# res_shrunk <- apply_lfc_shrinkage(analysis$dds, analysis$results)
# biomarkers <- extract_biomarkers(res_shrunk, padj_threshold = 0.05, 
#                                 lfc_threshold = 1, top_n = 100)
# save_deg_results(res_shrunk, biomarkers, dataset_name = "GSE12345")
# generate_summary(res_shrunk, biomarkers)

# ============================================================================
# STEP 4: Visualization
# ============================================================================
cat("STEP 4: Visualization\n")
cat("---------------------\n")
cat("Generate visualizations:\n\n")

# Example (uncomment and modify):
# res_df <- read.csv("results/deg_tables/GSE12345_all_results.csv")
# biomarkers <- read.csv("results/deg_tables/GSE12345_top_biomarkers.csv")
# create_volcano_plot(res_df, output_file = "results/plots/volcano.png")
# create_ma_plot(res_df, output_file = "results/plots/ma.png")
# create_heatmap(dds, biomarkers, n_genes = 50)
# create_pca_plot(dds, condition_col = "condition")
# create_distance_heatmap(dds)

# ============================================================================
# STEP 5: Pathway Enrichment
# ============================================================================
cat("STEP 5: Pathway Enrichment\n")
cat("--------------------------\n")
cat("Run pathway enrichment analysis:\n\n")

# Example (uncomment and modify):
# biomarkers <- read.csv("results/deg_tables/GSE12345_top_biomarkers.csv")
# biomarker_genes <- biomarkers$gene
# run_pathway_enrichment(biomarker_genes, output_prefix = "PD_biomarkers")

# ============================================================================
# STEP 6: Cross-Dataset Biomarker Identification
# ============================================================================
cat("STEP 6: Cross-Dataset Biomarker Identification\n")
cat("----------------------------------------------\n")
cat("Identify consensus biomarkers across multiple datasets:\n\n")

# Example (uncomment and modify):
# biomarker_files <- c(
#     "results/deg_tables/GSE12345_top_biomarkers.csv",
#     "results/deg_tables/GSE67890_top_biomarkers.csv"
# )
# dataset_names <- c("GSE12345", "GSE67890")
# consensus_result <- run_biomarker_identification(biomarker_files, dataset_names)

cat("\n========================================\n")
cat("Pipeline Setup Complete!\n")
cat("========================================\n")
cat("\nNext steps:\n")
cat("1. Download your RNA-seq datasets\n")
cat("2. Run each step sequentially\n")
cat("3. Review results in the results/ directory\n")
cat("4. Check biomarker rankings in results/biomarkers/\n\n")

