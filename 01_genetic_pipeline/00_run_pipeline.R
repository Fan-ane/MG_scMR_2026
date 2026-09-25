#!/usr/bin/env Rscript
# =============================================================================
# MG-scMR Pipeline: Single-Cell Immune eQTL-Driven Causal Gene Identification
#                   and Drug Target Discovery for Myasthenia Gravis
# =============================================================================
# Master pipeline script — sources all analysis steps in order
# =============================================================================

cat("========================================\n")
cat("MG-scMR Analysis Pipeline\n")
cat("========================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Set working directory
setwd("/ragene/workspace")

# Step 1: Data Preprocessing
cat(">>> Step 1: Data Preprocessing\n")
source("results/scripts/01_data_preprocessing.R")
cat("<<< Step 1 complete\n\n")

# Step 2: MR Analysis
cat(">>> Step 2: MR Analysis\n")
source("results/scripts/02_mr_analysis.R")
cat("<<< Step 2 complete\n\n")

# Step 3: Colocalization
cat(">>> Step 3: Colocalization\n")
source("results/scripts/03_colocalization.R")
cat("<<< Step 3 complete\n\n")

# Step 4: DICE Validation
cat(">>> Step 4: DICE Independent Validation\n")
source("results/scripts/05_dice_validation.R")
cat("<<< Step 4 complete\n\n")

# Step 5: PheWAS (Gene-level)
cat(">>> Step 5: PheWAS Gene-level\n")
source("results/scripts/06_phewas.R")
cat("<<< Step 5 complete\n\n")

# Step 5b: PheWAS (SNP-level)
cat(">>> Step 5b: PheWAS SNP-level\n")
source("results/scripts/06b_phewas_snp.R")
cat("<<< Step 5b complete\n\n")

# Step 6: Drug Repurposing
cat(">>> Step 6: Drug Repurposing\n")
source("results/scripts/07_drug_repurposing.R")
cat("<<< Step 6 complete\n\n")

# Step 7: Visualization (after PheWAS data is available)
cat(">>> Step 7: Visualization\n")
source("results/scripts/04_visualization.R")
cat("<<< Step 7 complete\n\n")

cat("========================================\n")
cat("Pipeline complete!\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n")
