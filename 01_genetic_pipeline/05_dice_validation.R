#!/usr/bin/env Rscript
# =============================================================================
# Step 5: DICE Independent Validation
# - Load DICE eQTL VCF files (10 immune cell types)
# - Match colocalized genes to DICE cell types
# - Test replication: P < 0.05 in DICE
# - Assess effect direction concordance where instruments overlap
# =============================================================================

library(data.table)

OUTPUT_DIR <- "data/processed"
RESULTS_DIR <- "results"
DICE_DIR <- "data/dice_eqtl"

# DICE to OneK1K cell type mapping
dice_to_onek1k <- c(
  "CD4_Naive" = "CD4 NC",
  "CD4_Stim" = "CD4 ET",
  "CD8_Naive" = "CD8 NC",
  "CD8_Stim" = "CD8 ET",
  "NK" = "NK",
  "B_Naive" = "B IN",
  "B_Mem" = "B Mem",
  "Mono_Classical" = "Mono C",
  "Mono_NonClassical" = "Mono NC",
  "DC" = "DC"
)

# Parse DICE VCF file
parse_dice_vcf <- function(vcf_path) {
  dt <- fread(vcf_path, skip = "#CHROM", header = TRUE, sep = "\t")
  names(dt)[1] <- "CHROM"
  dt$GeneSymbol <- gsub(".*GeneSymbol=([^;]+).*", "\\1", dt$INFO)
  dt$Pvalue <- as.numeric(gsub(".*Pvalue=([^;]+).*", "\\1", dt$INFO))
  dt$Beta <- as.numeric(gsub(".*Beta=([^;]+).*", "\\1", dt$INFO))
  dt[, .(rsid = ID, gene = GeneSymbol, pvalue = Pvalue, beta = Beta)]
}

# Load colocalized results (PP.H4 > 50%)
coloc_merged <- readRDS(file.path(OUTPUT_DIR, "coloc_merged.rds"))
colocalized <- coloc_merged[PP.H4 > 0.5]
instruments <- readRDS(file.path(OUTPUT_DIR, "instruments.rds"))

cat(sprintf("Testing %d colocalized gene-cell pairs in DICE\n", nrow(colocalized)))

# --- Validation ---
validation_results <- list()
val_idx <- 0

for (dice_ct in names(dice_to_onek1k)) {
  onek1k_ct <- dice_to_onek1k[dice_ct]

  # Genes to test in this cell type
  genes_to_test <- colocalized[cell_type == onek1k_ct]$gene
  if (length(genes_to_test) == 0) next

  # Load DICE VCF once per cell type
  vcf_file <- file.path(DICE_DIR, paste0(dice_ct, "_eQTL.vcf.gz"))
  if (!file.exists(vcf_file)) {
    cat(sprintf("  DICE file not found: %s\n", vcf_file))
    next
  }

  cat(sprintf("  Loading DICE %s (%d genes to test)...\n", dice_ct, length(genes_to_test)))
  dice_data <- parse_dice_vcf(vcf_file)

  for (gene_val in genes_to_test) {
    gene_hits <- dice_data[gene == gene_val]
    if (nrow(gene_hits) == 0) next

    # Best DICE association for this gene
    best_hit <- gene_hits[which.min(pvalue)]

    # Check instrument overlap
    inst_snps <- instruments[GENE == gene_val & cell_type == onek1k_ct,
                             .(rsid = RSID, onek1k_beta = beta_eqtl)]
    snp_overlap <- merge(inst_snps, gene_hits[, .(rsid, beta)],
                         by = "rsid")

    effect_concordant <- NA
    n_overlap <- nrow(snp_overlap)
    if (n_overlap > 0) {
      effect_concordant <- all(sign(snp_overlap$onek1k_beta) == sign(snp_overlap$beta))
    }

    val_idx <- val_idx + 1
    validation_results[[val_idx]] <- data.table(
      gene = gene_val,
      cell_type_onek1k = onek1k_ct,
      cell_type_dice = dice_ct,
      dice_n_snps = nrow(gene_hits),
      dice_lead_pval = best_hit$pvalue,
      dice_lead_beta = best_hit$beta,
      dice_lead_rsid = best_hit$rsid,
      instrument_overlap = n_overlap,
      effect_concordant = effect_concordant,
      validated = ifelse(best_hit$pvalue < 0.05, "Replicated", "Not replicated")
    )
  }
}

dice_validation <- rbindlist(validation_results)

cat(sprintf("\nDICE Validation Results:\n"))
cat(sprintf("  Pairs tested: %d\n", nrow(dice_validation)))
cat(sprintf("  Replicated (P < 0.05): %d (%.0f%%)\n",
            sum(dice_validation$validated == "Replicated"),
            100 * mean(dice_validation$validated == "Replicated")))
n_concordant <- sum(dice_validation$effect_concordant == TRUE, na.rm = TRUE)
n_testable <- sum(!is.na(dice_validation$effect_concordant))
cat(sprintf("  Effect direction concordance: %d/%d (%.0f%%)\n",
            n_concordant, n_testable, 100 * n_concordant / n_testable))

# --- Save ---
saveRDS(dice_validation, file.path(OUTPUT_DIR, "dice_validation.rds"))
fwrite(dice_validation, file.path(RESULTS_DIR, "Supplementary_Table_DICE_validation.csv"))

cat("Step 5 complete.\n")
