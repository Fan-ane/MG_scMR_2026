#!/usr/bin/env Rscript
# =============================================================================
# Step 1: Data Preprocessing
# - Load OneK1K sc-eQTL eSNP data (14 immune cell types)
# - Convert Spearman's rho to beta/SE
# - Filter: P < 5e-8, F-statistic > 10
# - LD clumping: r2 < 0.001, 10Mb window
# - Harmonize with MG GWAS (Braun et al. 2024)
# =============================================================================

library(data.table)
# Note: TwoSampleMR is loaded only for ld_clump() via genetics.binaRies::get_plink_binary()
suppressPackageStartupMessages(library(TwoSampleMR))

# --- Configuration ---
DATA_DIR <- "data"
ONEK1K_ESNP_DIR <- file.path(DATA_DIR, "onek1k_esnp")
GWAS_FILE <- file.path(DATA_DIR, "mg_gwas/GCST90432156.h.tsv.gz")
OUTPUT_DIR <- "data/processed"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Cell type metadata (sample sizes from OneK1K)
cell_types_meta <- data.table(
  cell_type = c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 NC", "CD8 ET",
                "CD8 S100B", "NK", "NK R", "B IN", "B Mem",
                "Mono C", "Mono NC", "DC", "Plasma"),
  n_samples = c(982L, 982L, 857L, 982L, 982L,
                981L, 928L, 967L, 982L, 982L,
                969L, 932L, 968L, 643L),
  esnp_file = paste0(c("CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
                        "CD8_S100B", "NK", "NK_R", "B_IN", "B_Mem",
                        "Mono_C", "Mono_NC", "DC", "Plasma"), "_eSNPs.tsv")
)
saveRDS(cell_types_meta, file.path(OUTPUT_DIR, "cell_types_meta.rds"))

# --- Load and process each cell type ---
cat("Loading OneK1K eSNP data...\n")

all_instruments <- list()
for (i in seq_len(nrow(cell_types_meta))) {
  ct <- cell_types_meta$cell_type[i]
  n <- cell_types_meta$n_samples[i]
  f <- file.path(ONEK1K_ESNP_DIR, cell_types_meta$esnp_file[i])

  cat(sprintf("  Processing %s (n=%d)...\n", ct, n))

  esnp <- fread(f)

  # Spearman's rho to beta/SE conversion:
  # beta ≈ rho, SE ≈ sqrt((1 - rho^2) / (n - 2))
  esnp[, beta_eqtl := CORRELATION]
  esnp[, se_eqtl := sqrt((1 - CORRELATION^2) / (n - 2))]
  esnp[, pval_eqtl := P_VALUE]

  # F-statistic = (beta / se)^2
  esnp[, fstat := (beta_eqtl / se_eqtl)^2]

  # Filter: P < 5e-8 and F > 10
  esnp_sig <- esnp[pval_eqtl < 5e-8 & fstat > 10]

  if (nrow(esnp_sig) == 0) {
    cat(sprintf("    No significant instruments for %s\n", ct))
    next
  }

  esnp_sig[, cell_type := ct]
  esnp_sig[, n_samples := n]

  all_instruments[[ct]] <- esnp_sig
  cat(sprintf("    %d significant eSNP-gene pairs\n", nrow(esnp_sig)))
}

instruments <- rbindlist(all_instruments, fill = TRUE)
cat(sprintf("Total instruments before clumping: %d\n", nrow(instruments)))

# --- LD Clumping ---
cat("Performing LD clumping (r2 < 0.001, 10Mb window)...\n")

# Format for TwoSampleMR clumping
clump_input <- data.frame(
  SNP = instruments$RSID,
  pval = instruments$pval_eqtl,
  id = paste(instruments$GENE, instruments$cell_type, sep = "_")
)

clumped <- ld_clump(
  clump_input,
  clump_kb = 10000,
  clump_r2 = 0.001,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "data/ld_ref/EUR"
)

instruments <- instruments[RSID %in% clumped$SNP]
cat(sprintf("Instruments after clumping: %d\n", nrow(instruments)))

# Note: GWAS outcome data is loaded and merged on-demand in script 02 (MR)
# and script 03 (coloc). We persist only the eQTL instruments here.
saveRDS(instruments, file.path(OUTPUT_DIR, "instruments.rds"))
cat("Step 1 complete.\n")
