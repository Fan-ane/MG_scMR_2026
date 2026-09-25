# ============================================================
# SESN3 Virtual Knockout — scTenifoldKnk Analysis
# Part 2: Virtual KO of SESN3 in CD4/CD8 T cells from HC
# ============================================================
#
# References:
#   scTenifoldKnk: https://github.com/cailab-tamu/scTenifoldKnk
#   Paper: Osorio D et al., Patterns 2022
#
# ============================================================

library(Seurat)
library(scTenifoldKnk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)

set.seed(123)
n_cores <- 8

# ---- 0. Load processed Seurat object ----
data_dir <- "c:/Users/zyf/Desktop/课题一生信文章/验证部分工作/scRNA_analysis"
out_dir <- file.path(data_dir, "scTenifoldKnk_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(file.path(data_dir, "merged_seurat_annotated.rds"))

# ---- 1. Extract HC CD4 T cells ----
# (Adjust cell type labels based on your annotation)

cd4_hc_cells <- WhichCells(obj, expression = condition == "HC" & cell_type == "CD4 T cells")
cd8_hc_cells <- WhichCells(obj, expression = condition == "HC" & cell_type == "CD8 T cells")

cat("HC CD4 T cells:", length(cd4_hc_cells), "\n")
cat("HC CD8 T cells:", length(cd8_hc_cells), "\n")

# Subsample per donor (max 500 cells each) to avoid donor imbalance
subsample_donor <- function(obj, cells, max_per_donor = 500) {
  meta <- obj@meta.data[cells, ]
  meta$barcode <- rownames(meta)
  sampled <- meta %>%
    group_by(sample) %>%
    slice_sample(n = min(max_per_donor, n())) %>%
    pull(barcode)
  return(sampled)
}

cd4_hc_sampled <- subsample_donor(obj, cd4_hc_cells)
cd8_hc_sampled <- subsample_donor(obj, cd8_hc_cells)

cat("Sampled CD4:", length(cd4_hc_sampled), "| CD8:", length(cd8_hc_sampled), "\n")

# ---- 2. Prepare expression matrix ----

prepare_matrix <- function(obj, cells, n_hvg = 4000) {
  # Subset to selected cells
  sub <- subset(obj, cells = cells)

  # Refind variable features
  sub <- NormalizeData(sub)
  sub <- FindVariableFeatures(sub, nfeatures = n_hvg)

  # Keep HVGs + SESN3
  genes_use <- unique(c(VariableFeatures(sub), "SESN3"))
  genes_use <- intersect(genes_use, rownames(sub))

  # Get raw counts matrix
  counts <- GetAssayData(sub, assay = "RNA", slot = "counts")
  counts <- counts[genes_use, ]

  cat("Matrix dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")
  return(counts)
}

cat("\nPreparing CD4 T-cell matrix...\n")
counts_cd4 <- prepare_matrix(obj, cd4_hc_sampled)

cat("\nPreparing CD8 T-cell matrix...\n")
counts_cd8 <- prepare_matrix(obj, cd8_hc_sampled)

# ---- 3. Run scTenifoldKnk: SESN3 virtual KO ----

run_sesn3_ko <- function(counts, cell_type_label, seed = 123) {
  set.seed(seed)

  cat("\n========================================\n")
  cat("Running SESN3 KO in", cell_type_label, "\n")
  cat("========================================\n")

  # Ensure SESN3 is in the matrix
  if (!"SESN3" %in% rownames(counts)) {
    stop("SESN3 NOT FOUND IN MATRIX! Check gene names.")
  }

  cat("SESN3 mean expression:", round(mean(counts["SESN3", ]), 2), "\n")

  # Run virtual KO
  result <- tryCatch({
    scTenifoldKnk(
      countMatrix = counts,
      gKO = "SESN3",
      qc = TRUE,
      nCores = n_cores
    )
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })

  return(result)
}

ko_cd4 <- run_sesn3_ko(counts_cd4, "CD4_Tcells_ HC")
ko_cd8 <- run_sesn3_ko(counts_cd8, "CD8_Tcells_HC")

# ---- 4. Extract and save differential regulation results ----

process_ko_result <- function(ko_result, label) {
  if (is.null(ko_result)) {
    cat("No result for", label, "\n")
    return(NULL)
  }

  dr <- ko_result$diffRegulation
  dr <- dr[order(dr$p.adj), ]

  # Save full results
  write.csv(dr, file.path(out_dir, paste0("SESN3_KO_", label, "_diffRegulation.csv")),
            row.names = FALSE)

  cat(label, "- Top 10 perturbed genes:\n")
  print(head(dr[, c("gene", "distance", "p.adj")], 10))

  return(dr)
}

dr_cd4 <- process_ko_result(ko_cd4, "CD4_Tcells")
dr_cd8 <- process_ko_result(ko_cd8, "CD8_Tcells")

# ---- 5. Multiple iterations for stability ----

n_iterations <- 20
cat("\nRunning", n_iterations, "iterations for stability...\n")

perturbation_freq_cd4 <- list()
perturbation_freq_cd8 <- list()

for (i in 1:n_iterations) {
  cat("Iteration", i, "/", n_iterations, "...\n")

  # Resample cells
  cd4_samp <- subsample_donor(obj, cd4_hc_cells)
  cd8_samp <- subsample_donor(obj, cd8_hc_cells)

  # Prepare matrices
  cnt_cd4 <- prepare_matrix(obj, cd4_samp)
  cnt_cd8 <- prepare_matrix(obj, cd8_samp)

  # Run KO
  set.seed(i * 100 + 123)
  res_cd4 <- tryCatch(
    scTenifoldKnk(countMatrix = cnt_cd4, gKO = "SESN3", qc = TRUE, nCores = n_cores),
    error = function(e) NULL
  )
  res_cd8 <- tryCatch(
    scTenifoldKnk(countMatrix = cnt_cd8, gKO = "SESN3", qc = TRUE, nCores = n_cores),
    error = function(e) NULL
  )

  # Store significant genes (p.adj < 0.05)
  if (!is.null(res_cd4)) {
    sig <- res_cd4$diffRegulation$gene[res_cd4$diffRegulation$p.adj < 0.05]
    perturbation_freq_cd4[[i]] <- sig
  }
  if (!is.null(res_cd8)) {
    sig <- res_cd8$diffRegulation$gene[res_cd8$diffRegulation$p.adj < 0.05]
    perturbation_freq_cd8[[i]] <- sig
  }
}

# Calculate gene perturbation frequency
calc_freq <- function(freq_list, label) {
  all_genes <- unique(unlist(freq_list))
  freq_df <- data.frame(
    gene = all_genes,
    frequency = sapply(all_genes, function(g) mean(sapply(freq_list, function(x) g %in% x))),
    stringsAsFactors = FALSE
  )
  freq_df <- freq_df[order(freq_df$frequency, decreasing = TRUE), ]
  write.csv(freq_df, file.path(out_dir, paste0("SESN3_KO_", label, "_perturbation_frequency.csv")),
            row.names = FALSE)
  cat(label, "- Stable perturbed genes (freq >= 0.5):", sum(freq_df$frequency >= 0.5), "\n")
  return(freq_df)
}

freq_cd4 <- calc_freq(perturbation_freq_cd4, "CD4_Tcells")
freq_cd8 <- calc_freq(perturbation_freq_cd8, "CD8_Tcells")

# ---- 6. GO enrichment of perturbed genes ----

library(clusterProfiler)
library(org.Hs.eg.db)

run_enrichment <- function(freq_df, cell_label, min_freq = 0.5) {
  genes_stable <- freq_df$gene[freq_df$frequency >= min_freq]

  if (length(genes_stable) < 5) {
    cat("Not enough stable genes for enrichment in", cell_label, "\n")
    return(NULL)
  }

  # Convert to Entrez
  gene_ids <- bitr(genes_stable, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db, drop = TRUE)

  if (nrow(gene_ids) < 5) {
    cat("Not enough mapped genes for enrichment in", cell_label, "\n")
    return(NULL)
  }

  # GO BP enrichment
  ego_bp <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.1)

  if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
    write.csv(as.data.frame(ego_bp),
              file.path(out_dir, paste0("SESN3_KO_", cell_label, "_GO_BP.csv")),
              row.names = FALSE)

    # Dot plot
    pdf(file.path(out_dir, paste0("SESN3_KO_", cell_label, "_GO_dotplot.pdf")),
        width = 10, height = 6)
    print(dotplot(ego_bp, showCategory = 15) + ggtitle(paste("SESN3 KO —", cell_label)))
    dev.off()
  }

  return(ego_bp)
}

go_cd4 <- run_enrichment(freq_cd4, "CD4_Tcells")
go_cd8 <- run_enrichment(freq_cd8, "CD8_Tcells")

# ---- 7. MG DEG overlap ----

# Get MG vs HC DEGs in CD4/CD8 T cells
find_tcell_degs <- function(obj, ct) {
  sub <- subset(obj, subset = cell_type == ct)
  Idents(sub) <- "condition"
  degs <- FindMarkers(sub, ident.1 = "MG", ident.2 = "HC",
                      logfc.threshold = 0.25, min.pct = 0.1)
  degs$gene <- rownames(degs)
  return(degs)
}

degs_cd4 <- find_tcell_degs(obj, "CD4 T cells")
degs_cd8 <- find_tcell_degs(obj, "CD8 T cells")

# Overlap with stable perturbed genes
overlap_test <- function(perturbed_genes, degs, label) {
  perturbed <- perturbed_genes[perturbed_genes$frequency >= 0.5, "gene"]

  deg_sig <- degs$gene[degs$p_val_adj < 0.05]
  all_genes <- degs$gene

  overlap <- intersect(perturbed, deg_sig)
  cat(label, "- Overlap:", length(overlap), "genes:", paste(overlap[1:min(10, length(overlap))], collapse = ", "), "\n")

  # Hypergeometric test
  n_overlap <- length(overlap)
  n_perturbed <- length(perturbed)
  n_deg <- length(deg_sig)
  n_total <- length(all_genes)

  p_val <- phyper(n_overlap - 1, n_perturbed, n_total - n_perturbed, n_deg,
                  lower.tail = FALSE)
  cat(label, "- Enrichment p-value:", p_val, "\n")

  return(data.frame(
    cell_type = label,
    n_perturbed = n_perturbed,
    n_deg = n_deg,
    n_overlap = n_overlap,
    p_hypergeom = p_val,
    overlap_genes = paste(overlap, collapse = ";"),
    stringsAsFactors = FALSE
  ))
}

overlap_cd4 <- overlap_test(freq_cd4, degs_cd4, "CD4_Tcells")
overlap_cd8 <- overlap_test(freq_cd8, degs_cd8, "CD8_Tcells")

overlap_summary <- rbind(overlap_cd4, overlap_cd8)
write.csv(overlap_summary, file.path(out_dir, "SESN3_KO_MG_DEG_overlap.csv"), row.names = FALSE)

# ---- 8. Negative control: random gene with similar expression ----

# Pick a gene with similar mean expression to SESN3 but no MR signal
sesn3_mean_cd4 <- mean(counts_cd4["SESN3", ])
all_gene_means <- rowMeans(counts_cd4)
# Find genes with expression within 10% of SESN3
nearby_genes <- names(all_gene_means)[
  abs(all_gene_means - sesn3_mean_cd4) / sesn3_mean_cd4 < 0.1
]
# Exclude known MR hits
known_genes <- c("SESN3", "AGER", "ZSCAN12", "CFD", "DXO", "ZBTB9",
                 "HLA-DOA", "BTN2A1", "HCG23", "TYMP", "TNFRSF14")
control_candidates <- setdiff(nearby_genes, known_genes)
# Pick 2 random controls
if (length(control_candidates) >= 2) {
  neg_control_genes <- sample(control_candidates, 2)
} else {
  neg_control_genes <- sample(setdiff(rownames(counts_cd4), known_genes), 2)
}

cat("\nNegative control genes:", paste(neg_control_genes, collapse = ", "), "\n")

# Run KO for each negative control
run_neg_control_ko <- function(counts, gene_name, label) {
  cat("Running negative control KO:", gene_name, "in", label, "\n")
  result <- tryCatch({
    scTenifoldKnk(countMatrix = counts, gKO = gene_name, qc = TRUE, nCores = n_cores)
  }, error = function(e) {
    cat("ERROR for", gene_name, ":", e$message, "\n")
    return(NULL)
  })
  return(result)
}

neg_controls <- list()
for (nc_gene in neg_control_genes) {
  neg_controls[[paste0("CD4_", nc_gene)]] <- run_neg_control_ko(counts_cd4, nc_gene, "CD4")
  neg_controls[[paste0("CD8_", nc_gene)]] <- run_neg_control_ko(counts_cd8, nc_gene, "CD8")
}

# Compare perturbation magnitude
perturbation_summary <- data.frame(
  KO = c("SESN3_CD4", "SESN3_CD8",
         paste0(neg_control_genes[1], "_CD4"), paste0(neg_control_genes[1], "_CD8"),
         paste0(neg_control_genes[2], "_CD4"), paste0(neg_control_genes[2], "_CD8")),
  n_perturbed = c(
    if (!is.null(dr_cd4)) sum(dr_cd4$p.adj < 0.05) else NA,
    if (!is.null(dr_cd8)) sum(dr_cd8$p.adj < 0.05) else NA,
    if (!is.null(neg_controls[[paste0("CD4_", neg_control_genes[1])]]))
      sum(neg_controls[[paste0("CD4_", neg_control_genes[1])]]$diffRegulation$p.adj < 0.05) else NA,
    if (!is.null(neg_controls[[paste0("CD8_", neg_control_genes[1])]]))
      sum(neg_controls[[paste0("CD8_", neg_control_genes[1])]]$diffRegulation$p.adj < 0.05) else NA,
    if (!is.null(neg_controls[[paste0("CD4_", neg_control_genes[2])]]))
      sum(neg_controls[[paste0("CD4_", neg_control_genes[2])]]$diffRegulation$p.adj < 0.05) else NA,
    if (!is.null(neg_controls[[paste0("CD8_", neg_control_genes[2])]]))
      sum(neg_controls[[paste0("CD8_", neg_control_genes[2])]]$diffRegulation$p.adj < 0.05) else NA
  ),
  stringsAsFactors = FALSE
)

write.csv(perturbation_summary, file.path(out_dir, "SESN3_KO_negative_control_comparison.csv"),
          row.names = FALSE)

# ---- 9. Monocyte control (cell-type specificity) ----

# Run SESN3 KO in monocytes (should show weaker perturbation than T cells)
mono_hc_cells <- WhichCells(obj, expression = condition == "HC" &
                            cell_type %in% c("CD14 Monocytes", "FCGR3A Monocytes", "Monocytes"))

if (length(mono_hc_cells) > 50) {
  mono_sampled <- subsample_donor(obj, mono_hc_cells)
  counts_mono <- prepare_matrix(obj, mono_sampled)

  cat("\nRunning SESN3 KO in monocytes (cell-type specificity control)...\n")
  ko_mono <- run_sesn3_ko(counts_mono, "Monocytes_HC")

  if (!is.null(ko_mono)) {
    dr_mono <- process_ko_result(ko_mono, "Monocytes")
    # Compare: T cells should show more perturbation than monocytes
    n_tcell <- if (!is.null(dr_cd4)) sum(dr_cd4$p.adj < 0.05) else 0
    n_mono <- if (!is.null(dr_mono)) sum(dr_mono$p.adj < 0.05) else 0
    cat("\nCell-type specificity: CD4 T =", n_tcell, "| Monocytes =", n_mono, "\n")
  }
}

# ---- 10. Save everything ----

save.image(file.path(out_dir, "SESN3_KO_workspace.RData"))
cat("\n========================================\n")
cat("scTenifoldKnk analysis complete!\n")
cat("Output directory:", out_dir, "\n")
cat("========================================\n")
