#!/usr/bin/env Rscript
# =============================================================================
# Step 2: Mendelian Randomization Analysis
# - Wald ratio (single SNP instruments)
# - IVW (multiple SNP instruments)
# - FDR correction (Benjamini-Hochberg)
# - Odds ratio calculation
# =============================================================================

library(data.table)
library(TwoSampleMR)

DATA_DIR <- "data"
OUTPUT_DIR <- "data/processed"
RESULTS_DIR <- "results"

# Load harmonized data (instruments + GWAS outcome merged during harmonization)
instruments <- readRDS(file.path(OUTPUT_DIR, "instruments.rds"))
cell_types_meta <- readRDS(file.path(OUTPUT_DIR, "cell_types_meta.rds"))

# Load GWAS for outcome data
gwas <- fread(file.path(DATA_DIR, "mg_gwas/GCST90432156.h.tsv.gz"))

# Merge GWAS outcome columns into instruments via RSID
gwas_subset <- gwas[, .(hm_rsid, beta_gwas = hm_beta, se_gwas = standard_error,
                         pval_gwas = p_value)]
instruments <- merge(instruments, gwas_subset, by.x = "RSID", by.y = "hm_rsid",
                     all.x = TRUE)
instruments <- instruments[!is.na(beta_gwas)]  # Keep only matched SNPs

# --- MR Analysis ---
cat("Running MR analysis...\n")

# Group instruments by gene-cell_type pair
gene_ct_pairs <- unique(instruments[, .(GENE, cell_type)])
cat(sprintf("Testing %d gene-cell type pairs\n", nrow(gene_ct_pairs)))

mr_results <- list()

for (i in seq_len(nrow(gene_ct_pairs))) {
  gene <- gene_ct_pairs$GENE[i]
  ct <- gene_ct_pairs$cell_type[i]

  # Get instruments for this gene-cell pair
  inst <- instruments[GENE == gene & cell_type == ct]
  nsnps <- nrow(inst)

  if (nsnps == 1) {
    # Wald ratio: beta_MR = beta_outcome / beta_exposure
    b_mr <- inst$beta_gwas / inst$beta_eqtl
    se_mr <- abs(inst$se_gwas / inst$beta_eqtl)
    pval_mr <- 2 * pnorm(-abs(b_mr / se_mr))
    method <- "Wald ratio"
  } else {
    # IVW (Burgess 2013): beta_IVW = sum(bx*by/se_y^2) / sum(bx^2/se_y^2)
    bx <- inst$beta_eqtl; by <- inst$beta_gwas; sy2 <- inst$se_gwas^2
    denom <- sum(bx^2 / sy2)
    b_mr <- sum(bx * by / sy2) / denom
    se_mr <- sqrt(1 / denom)
    pval_mr <- 2 * pnorm(-abs(b_mr / se_mr))
    method <- "Inverse variance weighted"
  }

  mr_results[[i]] <- data.table(
    id.exposure = gene,
    gene = gene,
    cell_type = ct,
    method = method,
    nsnp = nsnps,
    b = b_mr,
    se = se_mr,
    pval = pval_mr
  )
}

mr_all <- rbindlist(mr_results)
cat(sprintf("Total MR tests: %d\n", nrow(mr_all)))

# --- FDR Correction ---
mr_all[, fdr := p.adjust(pval, method = "BH")]

# --- Odds Ratio Calculation ---
mr_all[, or := exp(b)]
mr_all[, or_lci := exp(b - 1.96 * se)]
mr_all[, or_uci := exp(b + 1.96 * se)]

# --- Significant results ---
mr_sig <- mr_all[fdr < 0.05]
cat(sprintf("Significant pairs (FDR < 0.05): %d\n", nrow(mr_sig)))
cat(sprintf("Unique significant eGenes: %d\n", uniqueN(mr_sig$gene)))
cat(sprintf("Cell types with hits: %d/14\n", uniqueN(mr_sig$cell_type)))

# --- Save ---
saveRDS(mr_all, file.path(OUTPUT_DIR, "mr_results_all.rds"))
saveRDS(mr_sig, file.path(OUTPUT_DIR, "mr_results_sig_fdr.rds"))

# CSV exports
fwrite(mr_all, file.path(RESULTS_DIR, "Supplementary_Table_S2_MR_all.csv"))
fwrite(mr_sig, file.path(RESULTS_DIR, "Supplementary_Table_S2_MR_significant.csv"))

cat("Step 2 complete.\n")
