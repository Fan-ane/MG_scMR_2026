#!/usr/bin/env Rscript
# =============================================================================
# Step 3: Bayesian Colocalization Analysis
# - coloc.abf() for each significant gene-cell type pair
# - Batch extraction approach: one grep pass per cell type
# - RSID-based matching (OneK1K GRCh37 vs GWAS GRCh38)
# - Tier classification: PP.H4 > 80% = Tier 1, 50-80% = Tier 2
# =============================================================================

library(data.table)
library(coloc)

OUTPUT_DIR <- "data/processed"
RESULTS_DIR <- "results"
EQTL_DIR <- "data/onek1k_eqtl"
GWAS_FILE <- "data/mg_gwas/GCST90432156.h.tsv.gz"

# GWAS parameters
gwas_N <- 5708 + 432028
gwas_s <- 5708 / gwas_N

# Load significant MR results
sig_fdr <- readRDS(file.path(OUTPUT_DIR, "mr_results_sig_fdr.rds"))
cell_types_meta <- readRDS(file.path(OUTPUT_DIR, "cell_types_meta.rds"))

cat(sprintf("Testing colocalization for %d significant gene-cell pairs\n",
            nrow(sig_fdr)))

# --- Load GWAS summary statistics ---
cat("Loading GWAS data...\n")
gwas <- fread(GWAS_FILE)
gwas_rsids <- gwas$hm_rsid

# --- Batch extraction per cell type ---
cell_type_to_file <- c(
  "CD4 NC" = "CD4_NC", "CD4 ET" = "CD4_ET", "CD4 SOX4" = "CD4_SOX4",
  "CD8 NC" = "CD8_NC", "CD8 ET" = "CD8_ET", "CD8 S100B" = "CD8_S100B",
  "NK" = "NK", "NK R" = "NK_R", "B IN" = "B_IN", "B Mem" = "B_Mem",
  "Mono C" = "Mono_C", "Mono NC" = "Mono_NC", "DC" = "DC", "Plasma" = "Plasma"
)

coloc_results <- list()
result_idx <- 0

for (ct in unique(sig_fdr$cell_type)) {
  genes_in_ct <- sig_fdr[cell_type == ct]$gene
  ct_n <- cell_types_meta[cell_type == ct]$n_samples
  ct_file_prefix <- cell_type_to_file[ct]
  eqtl_file <- file.path(EQTL_DIR, paste0(ct_file_prefix, "_eQTLs.tsv.gz"))

  cat(sprintf("\n--- %s: %d genes to test ---\n", ct, length(genes_in_ct)))

  # Batch extract: one grep pass for all genes in this cell type.
  # Use -wF (fixed string) with patterns from a temp file to avoid regex
  # metacharacter issues (e.g. dots in gene symbols like "KB-1980E6.3").
  gene_file <- tempfile(fileext = ".txt")
  writeLines(genes_in_ct, gene_file)
  cmd <- sprintf("zcat '%s' | grep -wFf '%s'", eqtl_file, gene_file)

  cat("  Extracting eQTL data (batch)...\n")
  eqtl_lines <- tryCatch(
    fread(cmd = cmd, header = FALSE),
    error = function(e) { cat("  Error:", e$message, "\n"); NULL }
  )
  unlink(gene_file)

  if (is.null(eqtl_lines) || nrow(eqtl_lines) == 0) {
    cat("  No eQTL data extracted, skipping cell type\n")
    next
  }

  # Read header from file and apply as column names
  header_cmd <- sprintf("zcat '%s' | head -1", eqtl_file)
  header <- fread(cmd = header_cmd, header = FALSE, nrows = 1)
  setnames(eqtl_lines, as.character(unlist(header[1, ])))

  # Process each gene
  for (gene in genes_in_ct) {
    eqtl_gene <- eqtl_lines[GENE == gene]
    if (nrow(eqtl_gene) == 0) next

    # RSID-based matching with GWAS
    common_rsids <- intersect(eqtl_gene$RSID, gwas_rsids)
    if (length(common_rsids) < 10) {
      cat(sprintf("  %s: only %d overlapping SNPs, skipping\n",
                  gene, length(common_rsids)))
      next
    }

    eqtl_sub <- eqtl_gene[RSID %in% common_rsids]
    gwas_sub <- gwas[hm_rsid %in% common_rsids]

    # Align by RSID
    setkey(eqtl_sub, RSID)
    setkey(gwas_sub, hm_rsid)
    snp_ids <- intersect(eqtl_sub$RSID, gwas_sub$hm_rsid)

    eqtl_sub <- eqtl_sub[snp_ids]
    gwas_sub <- gwas_sub[snp_ids]

    # Compute beta/varbeta for eQTL (Spearman's rho conversion)
    n_eqtl <- ct_n
    eqtl_beta <- eqtl_sub$CORRELATION
    eqtl_se <- sqrt((1 - eqtl_sub$CORRELATION^2) / (n_eqtl - 2))
    eqtl_varbeta <- eqtl_se^2
    eqtl_maf <- pmin(eqtl_sub$A2_FREQ_ONEK1K, 1 - eqtl_sub$A2_FREQ_ONEK1K)

    # GWAS beta/varbeta
    gwas_beta <- gwas_sub$hm_beta
    gwas_varbeta <- gwas_sub$standard_error^2

    # Run coloc
    result <- tryCatch(
      coloc.abf(
        dataset1 = list(
          beta = eqtl_beta,
          varbeta = eqtl_varbeta,
          N = n_eqtl,
          MAF = eqtl_maf,
          type = "quant",
          snp = snp_ids
        ),
        dataset2 = list(
          beta = gwas_beta,
          varbeta = gwas_varbeta,
          N = gwas_N,
          s = gwas_s,
          type = "cc",
          snp = snp_ids
        )
      ),
      error = function(e) { cat(sprintf("  Error for %s: %s\n", gene, e$message)); NULL }
    )

    if (is.null(result)) next

    pp <- result$summary
    result_idx <- result_idx + 1
    coloc_results[[result_idx]] <- data.table(
      gene = gene,
      cell_type = ct,
      nsnps = length(snp_ids),
      PP.H0 = pp["PP.H0.abf"],
      PP.H1 = pp["PP.H1.abf"],
      PP.H2 = pp["PP.H2.abf"],
      PP.H3 = pp["PP.H3.abf"],
      PP.H4 = pp["PP.H4.abf"]
    )

    cat(sprintf("  %s: %d SNPs, PP.H4 = %.3f\n",
                gene, length(snp_ids), pp["PP.H4.abf"]))
  }
}

coloc_dt <- rbindlist(coloc_results)

# --- Tier classification ---
coloc_dt[, tier := ifelse(PP.H4 > 0.8, "Tier 1 (Strong)",
                   ifelse(PP.H4 > 0.5, "Tier 2 (Moderate)",
                          "Not colocalized"))]

cat(sprintf("\nColocalization results:\n"))
cat(sprintf("  Tier 1 (PP.H4 > 80%%): %d pairs\n", sum(coloc_dt$tier == "Tier 1 (Strong)")))
cat(sprintf("  Tier 2 (PP.H4 50-80%%): %d pairs\n", sum(coloc_dt$tier == "Tier 2 (Moderate)")))
cat(sprintf("  Not colocalized: %d pairs\n", sum(coloc_dt$tier == "Not colocalized")))

# --- Merge with MR results ---
coloc_merged <- merge(coloc_dt, sig_fdr[, .(gene, cell_type, b, se, pval, fdr, or, or_lci, or_uci)],
                      by = c("gene", "cell_type"), all.x = TRUE)

# --- Save ---
saveRDS(coloc_dt, file.path(OUTPUT_DIR, "coloc_results.rds"))
saveRDS(coloc_merged, file.path(OUTPUT_DIR, "coloc_merged.rds"))

# CSV layouts match delivered supplementary tables
coloc_csv <- coloc_dt[, .(cell_type, gene, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, tier)]
fwrite(coloc_csv, file.path(RESULTS_DIR, "Supplementary_Table_Coloc.csv"))
fwrite(coloc_merged, file.path(RESULTS_DIR, "Supplementary_Table_Coloc_merged.csv"))

cat("Step 3 complete.\n")
