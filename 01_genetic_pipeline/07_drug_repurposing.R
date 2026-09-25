#!/usr/bin/env Rscript
# =============================================================================
# Step 7: Drug Repurposing Analysis
# - Query Enrichr/DSigDB for drug-gene associations
# - Query Open Targets knownDrugs for approved/clinical drugs
# - Generate Table 1 with MR/coloc evidence + drug interactions
# =============================================================================

library(data.table)
library(httr)
library(jsonlite)
library(rlang)  # for %||% operator used in NULL-coalescing drug field access

OUTPUT_DIR <- "data/processed"
RESULTS_DIR <- "results"

# Load colocalized results
coloc_merged <- readRDS(file.path(OUTPUT_DIR, "coloc_merged.rds"))
colocalized <- coloc_merged[PP.H4 > 0.5]

unique_genes <- unique(colocalized$gene)
cat(sprintf("Querying drug databases for %d colocalized genes\n", length(unique_genes)))

# --- Enrichr/DSigDB Query ---
cat("Querying Enrichr/DSigDB...\n")

enrichr_add <- function(genes) {
  url <- "https://maayanlab.cloud/Enrichr/addList"
  res <- POST(url, body = list(list = paste(genes, collapse = "\n")))
  content(res)
}

enrichr_enrich <- function(user_list_id, library_name = "DSigDB") {
  url <- sprintf("https://maayanlab.cloud/Enrichr/enrich?userListId=%d&backgroundType=%s",
                 user_list_id, library_name)
  res <- GET(url)
  content(res)
}

# Submit gene list to Enrichr
add_result <- tryCatch(enrichr_add(unique_genes), error = function(e) NULL)

enrichr_drugs <- data.table()
if (!is.null(add_result) && !is.null(add_result$userListId)) {
  enrich_result <- tryCatch(
    enrichr_enrich(add_result$userListId),
    error = function(e) NULL
  )

  if (!is.null(enrich_result) && !is.null(enrich_result$DSigDB)) {
    for (entry in enrich_result$DSigDB) {
      term <- entry[[2]]  # Drug signature term (index 1=rank, 2=term name)
      pval <- entry[[3]]  # P-value
      genes_overlap <- entry[[6]]  # Overlapping genes

      enrichr_drugs <- rbind(enrichr_drugs, data.table(
        drug_signature = term,
        pvalue = pval,
        genes = paste(genes_overlap, collapse = ", ")
      ))
    }
    cat(sprintf("  Enrichr/DSigDB: %d drug signatures found\n", nrow(enrichr_drugs)))
  }
}

# --- Open Targets knownDrugs Query ---
cat("Querying Open Targets knownDrugs...\n")

OT_API <- "https://api.platform.opentargets.org/api/v4/graphql"

query_known_drugs <- function(gene_symbol) {
  # Get Ensembl ID
  search_query <- sprintf('{
    search(queryString: "%s", entityNames: ["target"]) {
      hits { id name }
    }
  }', gene_symbol)

  res <- POST(OT_API, body = list(query = search_query), encode = "json",
              content_type_json())
  if (status_code(res) != 200) return(NULL)

  search_data <- content(res)
  hits <- search_data$data$search$hits
  if (length(hits) == 0) return(NULL)

  ensembl_id <- hits[[1]]$id

  # Query known drugs
  drug_query <- sprintf('{
    target(ensemblId: "%s") {
      approvedSymbol
      knownDrugs(size: 100) {
        count
        rows {
          drug { id name }
          mechanismOfAction
          disease { name }
          phase
          status
          urls { name url }
        }
      }
    }
  }', ensembl_id)

  res2 <- POST(OT_API, body = list(query = drug_query), encode = "json",
               content_type_json())
  if (status_code(res2) != 200) return(NULL)

  drug_data <- content(res2)
  target <- drug_data$data$target
  if (is.null(target)) return(NULL)

  kd <- target$knownDrugs
  if (is.null(kd) || kd$count == 0) return(NULL)

  results <- lapply(kd$rows, function(row) {
    data.table(
      gene = gene_symbol,
      drug_name = row$drug$name %||% NA,
      mechanism = row$mechanismOfAction %||% NA,
      indication = row$disease$name %||% NA,
      phase = row$phase %||% NA,
      status = row$status %||% NA
    )
  })

  rbindlist(results, fill = TRUE)
}

# Query each gene
all_drugs <- list()
for (i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  cat(sprintf("  [%d/%d] %s...", i, length(unique_genes), gene))

  result <- tryCatch(
    query_known_drugs(gene),
    error = function(e) { cat(" ERROR:", e$message); NULL }
  )

  if (!is.null(result) && nrow(result) > 0) {
    all_drugs[[gene]] <- result
    cat(sprintf(" %d drugs\n", nrow(result)))
  } else {
    cat(" no drugs\n")
  }

  Sys.sleep(0.5)
}

ot_drugs <- rbindlist(all_drugs, fill = TRUE)
cat(sprintf("Open Targets knownDrugs: %d drug-gene pairs across %d genes\n",
            nrow(ot_drugs), uniqueN(ot_drugs$gene)))

# --- Build Table 1 ---
cat("Building Table 1: Drug Repurposing Targets...\n")

table1 <- colocalized[, .(gene, cell_type, tier, PP.H4, b, se, pval, fdr,
                           or, or_lci, or_uci)]

# Merge drug annotations
if (nrow(ot_drugs) > 0) {
  drug_summary <- ot_drugs[, .(
    n_drugs = uniqueN(drug_name),
    approved_drugs = paste(unique(drug_name[phase >= 4 | status == "Approved"]), collapse = "; "),
    drug_mechanisms = paste(unique(na.omit(mechanism)), collapse = "; "),
    indications = paste(unique(na.omit(indication)), collapse = "; ")
  ), by = gene]

  table1 <- merge(table1, drug_summary, by = "gene", all.x = TRUE)
} else {
  table1[, `:=`(n_drugs = 0L, approved_drugs = "", drug_mechanisms = "", indications = "")]
}

# Fill NA drug columns
table1[is.na(n_drugs), n_drugs := 0L]
table1[is.na(approved_drugs), approved_drugs := ""]
table1[is.na(drug_mechanisms), drug_mechanisms := ""]
table1[is.na(indications), indications := ""]

# Direction column appended at end (matches delivered Table 1 CSV layout)
table1[, direction := ifelse(b < 0,
  paste0("Protective (", "↑", " expression ", "→", " ", "↓", " MG risk)"),
  paste0("Risk (", "↑", " expression ", "→", " ", "↑", " MG risk)"))]

# Sort by tier and p-value
table1 <- table1[order(tier, pval)]

cat(sprintf("Table 1: %d gene-cell pairs, %d with known drugs\n",
            nrow(table1), sum(table1$n_drugs > 0)))

# --- Save ---
saveRDS(table1, file.path(OUTPUT_DIR, "table1_drug_repurposing.rds"))
fwrite(table1, file.path(RESULTS_DIR, "Table1_Drug_Repurposing_Targets.csv"))

cat("Step 7 complete.\n")
