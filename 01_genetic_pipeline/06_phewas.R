#!/usr/bin/env Rscript
# =============================================================================
# Step 6: Gene-Level PheWAS via Open Targets Platform API
# - Query associated diseases for each colocalized gene
# - Uses Open Targets Platform GraphQL API
# - Note: Originally specified as "AstraZeneca PheWAS Portal" in requirements;
#   Open Targets Platform API used as equivalent alternative because
#   AstraZeneca's PheWAS Portal does not provide a public programmatic API.
#   Open Targets is a comprehensive gene-disease association database that
#   encompasses and extends AstraZeneca PheWAS data.
# =============================================================================

library(data.table)
library(httr)
library(jsonlite)

OUTPUT_DIR <- "data/processed"
RESULTS_DIR <- "results"

# Load colocalized genes
coloc_merged <- readRDS(file.path(OUTPUT_DIR, "coloc_merged.rds"))
colocalized <- coloc_merged[PP.H4 > 0.5]
unique_genes <- unique(colocalized$gene)

cat(sprintf("Querying Open Targets for %d colocalized genes\n", length(unique_genes)))

# Open Targets GraphQL API
OT_API <- "https://api.platform.opentargets.org/api/v4/graphql"

query_open_targets <- function(gene_symbol) {
  # First, get Ensembl ID from gene symbol
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

  # Query associated diseases
  assoc_query <- sprintf('{
    target(ensemblId: "%s") {
      id
      approvedSymbol
      associatedDiseases(page: {index: 0, size: 500}) {
        count
        rows {
          disease { id name therapeuticAreas { id name } }
          score
          datasourceScores { id score }
        }
      }
    }
  }', ensembl_id)

  res2 <- POST(OT_API, body = list(query = assoc_query), encode = "json",
               content_type_json())

  if (status_code(res2) != 200) return(NULL)

  assoc_data <- content(res2)
  target <- assoc_data$data$target
  if (is.null(target)) return(NULL)

  rows <- target$associatedDiseases$rows
  if (length(rows) == 0) return(NULL)

  results <- lapply(rows, function(row) {
    ta <- row$disease$therapeuticAreas
    ta_names <- if (length(ta) > 0) paste(sapply(ta, function(x) x$name), collapse = "; ") else "Other"
    data.table(
      gene = gene_symbol,
      ensembl_id = ensembl_id,
      disease_id = row$disease$id,
      disease = row$disease$name,
      therapeutic_area = ta_names,
      overall_score = row$score
    )
  })

  rbindlist(results)
}

# Query each gene
all_phewas <- list()
for (i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  cat(sprintf("  [%d/%d] %s...", i, length(unique_genes), gene))

  result <- tryCatch(
    query_open_targets(gene),
    error = function(e) { cat(" ERROR:", e$message); NULL }
  )

  if (!is.null(result) && nrow(result) > 0) {
    all_phewas[[gene]] <- result
    cat(sprintf(" %d associations\n", nrow(result)))
  } else {
    cat(" no results\n")
  }

  Sys.sleep(0.5)  # Rate limiting
}

phewas_gene <- rbindlist(all_phewas, fill = TRUE)
cat(sprintf("\nTotal gene-disease associations: %d across %d genes\n",
            nrow(phewas_gene), uniqueN(phewas_gene$gene)))

# --- Save ---
saveRDS(phewas_gene, file.path(OUTPUT_DIR, "phewas_gene_level.rds"))
fwrite(phewas_gene, file.path(RESULTS_DIR, "Supplementary_Table_PheWAS_gene.csv"))

cat("Step 6 complete.\n")
