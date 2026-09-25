#!/usr/bin/env Rscript
# =============================================================================
# Step 0: Data Download Helper
# - Downloads public data automatically:
#     * MG GWAS (Braun 2024, GWAS Catalog GCST90432156)
#     * 1000 Genomes EUR LD reference (MRC IEU)
# - Prints registration instructions for restricted data:
#     * OneK1K sc-eQTL (eSNP + eQTL files)
#     * DICE eQTL VCFs
# - Verifies file existence and basic integrity (size > 0)
# - Safe to re-run: skips already-downloaded files
# =============================================================================

cat("========================================\n")
cat("MG-scMR Data Download Helper\n")
cat("========================================\n\n")

# --- Configuration ---
DATA_DIR <- "data"
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DATA_DIR, "mg_gwas"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DATA_DIR, "ld_ref"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DATA_DIR, "onek1k_esnp"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DATA_DIR, "onek1k_eqtl"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DATA_DIR, "dice_eqtl"), recursive = TRUE, showWarnings = FALSE)

# Increase download timeout for large files (default 60s is too short)
options(timeout = max(3600, getOption("timeout")))

# --- Helper: download with skip-if-exists ---
download_if_missing <- function(url, dest, label, expected_min_mb = 0) {
  if (file.exists(dest)) {
    sz_mb <- file.info(dest)$size / 1024^2
    if (sz_mb >= expected_min_mb) {
      cat(sprintf("  [SKIP] %s already exists (%.1f MB)\n", label, sz_mb))
      return(invisible(TRUE))
    } else {
      cat(sprintf("  [WARN] %s exists but size %.1f MB < expected %d MB; re-downloading\n",
                  label, sz_mb, expected_min_mb))
      file.remove(dest)
    }
  }
  cat(sprintf("  [DL]   %s\n", label))
  cat(sprintf("         from: %s\n", url))
  cat(sprintf("         to:   %s\n", dest))
  result <- tryCatch(
    {
      download.file(url, dest, mode = "wb", quiet = FALSE)
      sz_mb <- file.info(dest)$size / 1024^2
      cat(sprintf("         done: %.1f MB\n", sz_mb))
      TRUE
    },
    error = function(e) {
      cat(sprintf("         FAIL: %s\n", e$message))
      if (file.exists(dest)) file.remove(dest)
      FALSE
    }
  )
  invisible(result)
}

# =============================================================================
# Part 1: MG GWAS (public, no registration required)
# =============================================================================
cat("--- Part 1: MG GWAS (Braun et al. 2024, GCST90432156) ---\n")

# Harmonised summary statistics from GWAS Catalog FTP
mg_gwas_url <- paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
  "GCST90432001-GCST90433000/GCST90432156/harmonised/",
  "GCST90432156.h.tsv.gz"
)
mg_gwas_dest <- file.path(DATA_DIR, "mg_gwas", "GCST90432156.h.tsv.gz")
download_if_missing(mg_gwas_url, mg_gwas_dest, "MG GWAS harmonised", expected_min_mb = 100)

cat("\n")

# =============================================================================
# Part 2: 1000 Genomes EUR LD reference (public, MRC IEU)
# =============================================================================
cat("--- Part 2: 1000 Genomes EUR LD reference ---\n")

ld_url <- "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
ld_archive <- file.path(DATA_DIR, "ld_ref", "1kg.v3.tgz")
ld_marker <- file.path(DATA_DIR, "ld_ref", "EUR.bim")

if (file.exists(ld_marker)) {
  cat("  [SKIP] EUR.bim already exists (LD reference present)\n")
} else {
  ok <- download_if_missing(ld_url, ld_archive, "1000G LD reference archive",
                            expected_min_mb = 500)
  if (ok && file.exists(ld_archive)) {
    cat("  [EXT]  Extracting EUR plink files...\n")
    untar(ld_archive, exdir = file.path(DATA_DIR, "ld_ref"))
    # Move EUR.* up if extracted into a subdirectory
    eur_files <- list.files(file.path(DATA_DIR, "ld_ref"),
                             pattern = "^EUR\\.(bed|bim|fam)$",
                             recursive = TRUE, full.names = TRUE)
    for (f in eur_files) {
      target <- file.path(DATA_DIR, "ld_ref", basename(f))
      if (f != target) file.rename(f, target)
    }
    if (file.exists(ld_marker)) {
      cat("  [OK]   EUR.{bed,bim,fam} ready\n")
      file.remove(ld_archive)  # cleanup archive
    } else {
      cat("  [WARN] Extraction complete but EUR.bim not found; check archive layout\n")
    }
  }
}

cat("\n")

# =============================================================================
# Part 3: OneK1K sc-eQTL (registration required)
# =============================================================================
cat("--- Part 3: OneK1K sc-eQTL ---\n")
cat("\n")
cat("OneK1K data requires registration at the OneK1K portal.\n")
cat("Please follow these steps manually:\n\n")
cat("  1. Register at:  https://onek1k.org/\n")
cat("  2. Navigate to:  Downloads -> sc-eQTL summary statistics\n")
cat("  3. For each of the 14 immune cell types listed below, download both:\n")
cat("       (a) The eSNP file (lead SNPs per gene)\n")
cat("       (b) The full eQTL file (all SNPs per gene)\n")
cat("  4. Place files into the local directories:\n")
cat(sprintf("       %s/<CELL>_eSNPs.tsv\n", file.path(DATA_DIR, "onek1k_esnp")))
cat(sprintf("       %s/<CELL>_eQTLs.tsv.gz\n", file.path(DATA_DIR, "onek1k_eqtl")))
cat("\n")
cat("  Required cell types (use these exact filename prefixes):\n")
ct_files <- c("CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "NK", "NK_R", "B_IN", "B_Mem",
              "Mono_C", "Mono_NC", "DC", "Plasma")
for (ct in ct_files) cat(sprintf("    - %s\n", ct))

# Check what is present
esnp_present <- list.files(file.path(DATA_DIR, "onek1k_esnp"), pattern = "_eSNPs\\.tsv$")
eqtl_present <- list.files(file.path(DATA_DIR, "onek1k_eqtl"), pattern = "_eQTLs\\.tsv\\.gz$")
cat(sprintf("\n  Currently present: %d/14 eSNP files, %d/14 eQTL files\n",
            length(esnp_present), length(eqtl_present)))

cat("\n")

# =============================================================================
# Part 4: DICE eQTL VCF (registration required)
# =============================================================================
cat("--- Part 4: DICE eQTL VCF ---\n")
cat("\n")
cat("DICE data requires registration at the DICE database.\n")
cat("Please follow these steps manually:\n\n")
cat("  1. Register at:  https://dice-database.org/\n")
cat("  2. Navigate to:  Downloads -> eQTL summary statistics (VCF format)\n")
cat("  3. Download VCF files for the 10 cell types listed below.\n")
cat("  4. Place files into:\n")
cat(sprintf("       %s/<CELL>_eQTL.vcf.gz\n", file.path(DATA_DIR, "dice_eqtl")))
cat("\n")
cat("  Required cell types (use these exact filename prefixes):\n")
dice_files <- c("CD4_Naive", "CD4_Stim", "CD8_Naive", "CD8_Stim", "NK",
                "B_Naive", "B_Mem", "Mono_Classical", "Mono_NonClassical", "DC")
for (dc in dice_files) cat(sprintf("    - %s\n", dc))

dice_present <- list.files(file.path(DATA_DIR, "dice_eqtl"), pattern = "_eQTL\\.vcf\\.gz$")
cat(sprintf("\n  Currently present: %d/10 DICE VCF files\n", length(dice_present)))

cat("\n")

# =============================================================================
# Final readiness check
# =============================================================================
cat("========================================\n")
cat("Data Readiness Summary\n")
cat("========================================\n\n")

mg_ok <- file.exists(mg_gwas_dest) && file.info(mg_gwas_dest)$size > 100 * 1024^2
ld_ok <- file.exists(ld_marker)
onek1k_ok <- length(esnp_present) == 14 && length(eqtl_present) == 14
dice_ok <- length(dice_present) == 10

cat(sprintf("  [%s] MG GWAS\n",       ifelse(mg_ok,     "OK", "MISSING")))
cat(sprintf("  [%s] 1000G EUR LD ref\n", ifelse(ld_ok,  "OK", "MISSING")))
cat(sprintf("  [%s] OneK1K eSNP+eQTL (need 14+14)\n",
            ifelse(onek1k_ok, "OK", sprintf("INCOMPLETE: %d+%d/14",
                                            length(esnp_present), length(eqtl_present)))))
cat(sprintf("  [%s] DICE eQTL VCF (need 10)\n",
            ifelse(dice_ok, "OK", sprintf("INCOMPLETE: %d/10", length(dice_present)))))

cat("\n")
if (mg_ok && ld_ok && onek1k_ok && dice_ok) {
  cat("All data ready. You can now run:\n")
  cat("  Rscript results/scripts/00_run_pipeline.R\n")
} else {
  cat("Some data is still missing. Re-run this script after manual downloads,\n")
  cat("or proceed with the parts of the pipeline that have data available.\n")
}
cat("\n")
