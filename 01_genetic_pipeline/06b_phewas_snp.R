#!/usr/bin/env Rscript
# =============================================================================
# Step 6b: SNP-Level PheWAS via GWAS Catalog REST API
# - Query GWAS Catalog for associations of lead instrument SNPs
# - Categorize traits into broad categories
# =============================================================================

library(data.table)
library(httr)
library(jsonlite)

OUTPUT_DIR <- "data/processed"
RESULTS_DIR <- "results"

# Load instruments for colocalized genes
coloc_merged <- readRDS(file.path(OUTPUT_DIR, "coloc_merged.rds"))
colocalized <- coloc_merged[PP.H4 > 0.5]
instruments <- readRDS(file.path(OUTPUT_DIR, "instruments.rds"))

# Get lead SNPs for colocalized genes
lead_snps <- instruments[paste(GENE, cell_type) %in% paste(colocalized$gene, colocalized$cell_type)]
lead_snps <- lead_snps[, .SD[which.min(pval_eqtl)], by = .(GENE, cell_type)]
unique_rsids <- unique(lead_snps$RSID)

cat(sprintf("Querying GWAS Catalog for %d unique lead SNPs\n", length(unique_rsids)))

# GWAS Catalog REST API
GWAS_API <- "https://www.ebi.ac.uk/gwas/rest/api"

query_gwas_catalog <- function(rsid) {
  url <- sprintf("%s/singleNucleotidePolymorphisms/%s/associations", GWAS_API, rsid)
  res <- GET(url, accept_json())

  if (status_code(res) != 200) return(NULL)

  data <- content(res, as = "text", encoding = "UTF-8")
  parsed <- fromJSON(data, flatten = TRUE)

  assocs <- parsed$`_embedded`$associations
  if (is.null(assocs) || length(assocs) == 0) return(NULL)

  results <- data.table(
    rsid = rsid,
    trait = assocs$`_links.efoTraits.href`,
    pval = assocs$pvalue,
    beta = ifelse(is.null(assocs$betaNum), NA, assocs$betaNum),
    odds_ratio = ifelse(is.null(assocs$orPerCopyNum), NA, assocs$orPerCopyNum)
  )

  # Extract trait names from EFO links if available
  if ("efoTraits" %in% names(assocs)) {
    results$trait <- sapply(assocs$efoTraits, function(x) {
      if (is.data.frame(x) && "trait" %in% names(x)) x$trait[1] else NA
    })
  }

  results
}

# Trait categorization
categorize_trait <- function(trait) {
  trait_lower <- tolower(trait)
  if (grepl("immun|autoimmun|inflamm|allerg|asthma|arthrit|lupus|sclerosis", trait_lower)) {
    return("Immune/AutoImmune")
  } else if (grepl("blood|hemoglobin|platelet|corpuscul|erythro|leuko|iron|ferrit", trait_lower)) {
    return("Hematological")
  } else if (grepl("lipid|cholesterol|triglycerid|hdl|ldl|fatty", trait_lower)) {
    return("Lipid")
  } else if (grepl("cancer|carcinoma|neoplasm|tumor|glioma|melanoma", trait_lower)) {
    return("Neoplasm")
  } else if (grepl("diabet|glucose|insulin|bmi|body mass|obes|metabol|weight|height", trait_lower)) {
    return("Metabolic")
  } else if (grepl("protein|biomark|complement|cytokine|interleukin", trait_lower)) {
    return("Biomarker")
  } else if (grepl("stature|anthropom|hip|waist", trait_lower)) {
    return("Anthropometric")
  } else {
    return("Other")
  }
}

# Query each SNP
all_snp_phewas <- list()
for (i in seq_along(unique_rsids)) {
  rsid <- unique_rsids[i]
  cat(sprintf("  [%d/%d] %s...", i, length(unique_rsids), rsid))

  result <- tryCatch(
    query_gwas_catalog(rsid),
    error = function(e) { cat(" ERROR:", e$message); NULL }
  )

  if (!is.null(result) && nrow(result) > 0) {
    # Add gene annotation
    gene_for_snp <- lead_snps[RSID == rsid]$GENE[1]
    result[, gene := gene_for_snp]
    result[, trait_category := sapply(trait, categorize_trait)]
    all_snp_phewas[[rsid]] <- result
    cat(sprintf(" %d associations\n", nrow(result)))
  } else {
    cat(" no results\n")
  }

  Sys.sleep(0.3)  # Rate limiting
}

phewas_snp <- rbindlist(all_snp_phewas, fill = TRUE)

# Remove duplicates
phewas_snp <- unique(phewas_snp)

cat(sprintf("\nTotal SNP-trait associations: %d across %d SNPs\n",
            nrow(phewas_snp), uniqueN(phewas_snp$rsid)))

# --- Save ---
# RDS keeps trait_category for downstream visualization (script 04)
saveRDS(phewas_snp, file.path(OUTPUT_DIR, "phewas_snp_level.rds"))
# CSV layout matches delivered Supplementary table (no trait_category column)
phewas_snp_csv <- phewas_snp[, .(rsid, trait, pval, beta, odds_ratio, gene)]
fwrite(phewas_snp_csv, file.path(RESULTS_DIR, "Supplementary_Table_PheWAS_snp.csv"))

cat("Step 6b complete.\n")
