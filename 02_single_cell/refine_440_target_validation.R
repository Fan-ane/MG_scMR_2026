suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript refine_440_target_validation.R <output_dir>")

output_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
results_dir <- file.path(output_dir, "results")
figures_dir <- file.path(output_dir, "figures")

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

obj <- readRDS(file.path(results_dir, "seurat_440_validation.rds"))
expr_mat <- GetAssayData(obj, assay = "RNA", slot = "data")

marker_sets <- list(
  B_IN = c("IGHD", "IGHM", "TCL1A", "FCER2", "IL4R", "CD24", "CD38"),
  B_Mem = c("CD27", "TNFRSF13B", "CD80", "CD86", "FCRL5", "IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGHA2"),
  Plasma = c("MZB1", "XBP1", "JCHAIN", "SDC1", "PRDM1", "TNFRSF17"),
  CD4_NC = c("CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LTB", "TCF7"),
  CD8_NC = c("CD3D", "CD3E", "CD8A", "CD8B", "CCR7", "LEF1", "TCF7"),
  CD8_ET = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "GZMH"),
  CD8_S100B = c("CD3D", "CD3E", "CD8A", "CD8B", "S100B", "GZMK")
)

score_set <- function(cells, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0 || length(cells) == 0) return(rep(0, length(cells)))
  as.numeric(Matrix::colMeans(expr_mat[genes, cells, drop = FALSE]))
}

refined <- as.character(obj$cell_type_detail)
names(refined) <- colnames(obj)

log_msg("Refining B-cell subtypes")
b_cells <- names(refined)[refined == "B_cell"]
if (length(b_cells) > 0) {
  b_scores <- data.frame(
    B_IN = score_set(b_cells, marker_sets$B_IN),
    B_Mem = score_set(b_cells, marker_sets$B_Mem),
    row.names = b_cells
  )
  b_best <- colnames(b_scores)[max.col(b_scores, ties.method = "first")]
  b_best[apply(b_scores, 1, max, na.rm = TRUE) <= 0] <- "B_cell"
  refined[b_cells] <- b_best
}

log_msg("Refining naive/cytotoxic CD8 T-cell labels")
t_cells <- names(refined)[refined %in% c("CD4_NC", "CD8_ET")]
if (length(t_cells) > 0) {
  t_scores <- data.frame(
    CD4_NC = score_set(t_cells, marker_sets$CD4_NC),
    CD8_NC = score_set(t_cells, marker_sets$CD8_NC),
    CD8_ET = score_set(t_cells, marker_sets$CD8_ET),
    CD8_S100B = score_set(t_cells, marker_sets$CD8_S100B),
    row.names = t_cells
  )
  t_best <- colnames(t_scores)[max.col(t_scores, ties.method = "first")]
  zero_score <- apply(t_scores, 1, max, na.rm = TRUE) <= 0
  t_best[zero_score] <- refined[t_cells][zero_score]
  refined[t_cells] <- t_best
}

obj$cell_type_detail_refined <- factor(refined)
write.csv(as.data.frame(table(obj$cell_type_detail, obj$cell_type_detail_refined)),
          file.path(results_dir, "celltype_refinement_crosswalk.csv"), row.names = FALSE)

prop <- as.data.frame(table(obj$sample, obj$group, obj$cell_type_detail_refined), stringsAsFactors = FALSE)
colnames(prop) <- c("sample", "group", "cell_type_detail", "cells")
prop <- prop %>%
  group_by(sample) %>%
  mutate(total_cells = sum(cells), proportion = cells / total_cells) %>%
  ungroup()
write.csv(prop, file.path(results_dir, "celltype_proportions_by_sample_refined.csv"), row.names = FALSE)

tier1_genes <- c("HCG23", "ZSCAN12", "AGER", "SCAND3", "ZNF322",
                 "SESN3", "TNFRSF14", "IL12RB2", "ORMDL3", "HFE", "CTSH")
drug_genes <- c("AGER", "CFD", "TYMP")
candidate_genes <- unique(c(tier1_genes, drug_genes))
candidate_present <- intersect(candidate_genes, rownames(obj))
candidate_missing <- setdiff(candidate_genes, rownames(obj))

targets <- data.frame(
  gene = c("HCG23", "ZSCAN12", "AGER", "SCAND3", "ZNF322", "SESN3", "SESN3",
           "TNFRSF14", "IL12RB2", "ORMDL3", "HFE", "CTSH", "CFD", "CFD", "CFD",
           "TYMP", "TYMP", "TYMP", "TYMP"),
  target_cell_type = c("B_IN", "Mono_C", "Mono_NC", "CD4_NC", "DC", "CD4_NC", "CD8_NC",
                       "CD8_ET", "NK", "B_Mem", "CD8_NC", "B_IN", "CD8_S100B", "CD8_ET", "NK",
                       "Mono_C", "B_Mem", "Mono_NC", "DC"),
  discovery_direction = c("protective", "risk", "risk", "protective", "protective", "protective", "protective",
                          "protective", "risk", "protective", "protective", "risk", "risk", "risk", "risk",
                          "protective", "protective", "protective", "protective"),
  stringsAsFactors = FALSE
)

calc_target <- function(gene, ct, discovery_direction) {
  present <- gene %in% rownames(obj)
  cells <- colnames(obj)[obj$cell_type_detail_refined == ct]
  n_hc <- sum(obj$group[cells] == "HC")
  n_mg <- sum(obj$group[cells] == "MG")
  out <- data.frame(
    gene = gene,
    target_cell_type = ct,
    discovery_direction = discovery_direction,
    expected_MG_direction = ifelse(discovery_direction == "risk", "up", "down"),
    gene_present = present,
    target_cells_HC = n_hc,
    target_cells_MG = n_mg,
    pct_expressing_HC = NA_real_,
    pct_expressing_MG = NA_real_,
    mean_expr_HC = NA_real_,
    mean_expr_MG = NA_real_,
    avg_log2FC_MG_vs_HC = NA_real_,
    p_val = NA_real_,
    p_val_adj = NA_real_,
    observed_MG_direction = NA_character_,
    expressed_in_target = FALSE,
    direction_consistent = NA,
    stringsAsFactors = FALSE
  )
  if (!present || length(cells) == 0 || n_hc < 5 || n_mg < 5) return(out)
  vals <- as.numeric(expr_mat[gene, cells])
  grp <- obj$group[cells]
  out$pct_expressing_HC <- mean(vals[grp == "HC"] > 0) * 100
  out$pct_expressing_MG <- mean(vals[grp == "MG"] > 0) * 100
  out$mean_expr_HC <- mean(vals[grp == "HC"])
  out$mean_expr_MG <- mean(vals[grp == "MG"])
  out$expressed_in_target <- max(out$pct_expressing_HC, out$pct_expressing_MG, na.rm = TRUE) >= 5
  sub_obj <- subset(obj, cells = cells)
  Idents(sub_obj) <- sub_obj$group
  fm <- tryCatch({
    FindMarkers(sub_obj, ident.1 = "MG", ident.2 = "HC", features = gene,
                test.use = "wilcox", logfc.threshold = 0, min.pct = 0, verbose = FALSE)
  }, error = function(e) NULL)
  if (!is.null(fm) && nrow(fm) > 0) {
    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(fm))[1]
    out$avg_log2FC_MG_vs_HC <- fm[1, fc_col]
    out$p_val <- fm[1, "p_val"]
    out$p_val_adj <- fm[1, "p_val_adj"]
    out$observed_MG_direction <- ifelse(out$avg_log2FC_MG_vs_HC > 0, "up",
                                        ifelse(out$avg_log2FC_MG_vs_HC < 0, "down", "flat"))
    out$direction_consistent <- out$observed_MG_direction == out$expected_MG_direction
  }
  out
}

target_validation <- bind_rows(mapply(
  calc_target,
  targets$gene,
  targets$target_cell_type,
  targets$discovery_direction,
  SIMPLIFY = FALSE
))
write.csv(target_validation, file.path(results_dir, "target_gene_validation_refined.csv"), row.names = FALSE)

log_msg("Recomputing DEG summaries for refined labels")
deg_summary <- list()
cell_types <- sort(unique(as.character(obj$cell_type_detail_refined)))
for (ct in cell_types) {
  cells <- colnames(obj)[obj$cell_type_detail_refined == ct]
  n_hc <- sum(obj$group[cells] == "HC")
  n_mg <- sum(obj$group[cells] == "MG")
  if (n_hc < 20 || n_mg < 20) next
  sub_obj <- subset(obj, cells = cells)
  Idents(sub_obj) <- sub_obj$group
  fm <- tryCatch({
    FindMarkers(sub_obj, ident.1 = "MG", ident.2 = "HC", test.use = "wilcox",
                min.pct = 0.1, logfc.threshold = 0.25, max.cells.per.ident = 5000,
                verbose = FALSE)
  }, error = function(e) NULL)
  if (is.null(fm) || nrow(fm) == 0) next
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(fm))[1]
  fm$gene <- rownames(fm)
  sig <- fm[!is.na(fm$p_val_adj) & fm$p_val_adj < 0.05, , drop = FALSE]
  write.csv(fm, file.path(results_dir, paste0("DEG_refined_", ct, "_MG_vs_HC.csv")), row.names = FALSE)
  deg_summary[[ct]] <- data.frame(
    cell_type_detail = ct,
    target_cells_HC = n_hc,
    target_cells_MG = n_mg,
    tested_genes = nrow(fm),
    sig_genes = nrow(sig),
    MG_up = sum(sig[[fc_col]] > 0),
    MG_down = sum(sig[[fc_col]] < 0),
    stringsAsFactors = FALSE
  )
}
deg_summary <- bind_rows(deg_summary)
if (nrow(deg_summary) > 0) {
  deg_summary <- deg_summary %>%
    mutate(MG_up_prop = ifelse(sig_genes > 0, MG_up / sig_genes, 0),
           MG_down_prop = ifelse(sig_genes > 0, MG_down / sig_genes, 0))
}
write.csv(deg_summary, file.path(results_dir, "celltype_deg_up_down_summary_refined.csv"), row.names = FALSE)

saveRDS(obj, file.path(results_dir, "seurat_440_validation_refined.rds"))

pC <- DimPlot(obj, reduction = "umap", group.by = "cell_type_detail_refined",
              label = TRUE, repel = TRUE, pt.size = 0.15) +
  ggtitle("C Refined immune cell types")
ggsave(file.path(figures_dir, "C_refined_celltypes_umap.pdf"), pC, width = 8, height = 6)
ggsave(file.path(figures_dir, "C_refined_celltypes_umap.png"), pC, width = 8, height = 6, dpi = 300)

pD <- ggplot(prop, aes(x = group, y = proportion, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_point(aes(shape = sample), position = position_jitter(width = 0.12, height = 0), size = 1.7) +
  facet_wrap(~ cell_type_detail, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell proportion", title = "D Refined cell-type proportions") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(figures_dir, "D_celltype_proportions_refined.pdf"), pD, width = 12, height = 9)
ggsave(file.path(figures_dir, "D_celltype_proportions_refined.png"), pD, width = 12, height = 9, dpi = 300)

if (length(candidate_present) > 0) {
  pE <- DotPlot(obj, features = candidate_present, group.by = "cell_type_detail_refined") +
    RotatedAxis() +
    labs(title = "E Candidate genes across refined immune cell types") +
    theme_bw(base_size = 10)
  ggsave(file.path(figures_dir, "E_candidate_gene_dotplot_refined.pdf"), pE, width = max(8, length(candidate_present) * 0.55), height = 5.8)
  ggsave(file.path(figures_dir, "E_candidate_gene_dotplot_refined.png"), pE, width = max(8, length(candidate_present) * 0.55), height = 5.8, dpi = 300)
}

if (nrow(deg_summary) > 0) {
  deg_long <- rbind(
    data.frame(cell_type_detail = deg_summary$cell_type_detail, direction = "MG_up", genes = deg_summary$MG_up),
    data.frame(cell_type_detail = deg_summary$cell_type_detail, direction = "MG_down", genes = -deg_summary$MG_down)
  )
  pG <- ggplot(deg_long, aes(x = reorder(cell_type_detail, genes), y = genes, fill = direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "Significant DEG count (FDR < 0.05)", title = "G Refined up/down DEG balance by cell type") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(figures_dir, "G_deg_up_down_by_celltype_refined.pdf"), pG, width = 8, height = 6)
  ggsave(file.path(figures_dir, "G_deg_up_down_by_celltype_refined.png"), pG, width = 8, height = 6, dpi = 300)
}

hdat <- target_validation %>%
  mutate(
    minus_log10_fdr = -log10(pmax(p_val_adj, 1e-300)),
    consistency = ifelse(is.na(direction_consistent), "not_tested",
                         ifelse(direction_consistent, "consistent", "opposite"))
  )
pH <- ggplot(hdat, aes(x = gene, y = target_cell_type)) +
  geom_point(aes(size = minus_log10_fdr, color = avg_log2FC_MG_vs_HC, shape = consistency), alpha = 0.9) +
  scale_color_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0,
                        na.value = "grey80", name = "log2FC\nMG vs HC") +
  scale_size_continuous(range = c(1.5, 7), name = "-log10 FDR") +
  labs(x = NULL, y = NULL, title = "H Refined target-cell validation") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
ggsave(file.path(figures_dir, "H_target_gene_validation_refined.pdf"), pH, width = 10.5, height = 6.5)
ggsave(file.path(figures_dir, "H_target_gene_validation_refined.png"), pH, width = 10.5, height = 6.5, dpi = 300)

answer1 <- target_validation %>% filter(gene_present, expressed_in_target) %>% distinct(gene, target_cell_type)
answer2 <- target_validation %>% filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  distinct(gene, target_cell_type, avg_log2FC_MG_vs_HC, p_val_adj)
answer3 <- target_validation %>% filter(!is.na(direction_consistent)) %>%
  count(direction_consistent, name = "n")

report_file <- file.path(output_dir, "report_validation_summary_refined.md")
writeLines(c(
  "# 440 refined single-cell validation summary",
  "",
  sprintf("- Cells retained after QC: %s.", format(ncol(obj), big.mark = ",")),
  sprintf("- Refined cell types: %s.", paste(sort(unique(as.character(obj$cell_type_detail_refined))), collapse = ", ")),
  sprintf("- Candidate genes present in matrix: %s.", paste(candidate_present, collapse = ", ")),
  if (length(candidate_missing) > 0) sprintf("- Candidate genes missing from matrix: %s.", paste(candidate_missing, collapse = ", ")) else "- Candidate genes missing from matrix: none.",
  "",
  "## Three document questions",
  "",
  sprintf("1. Expressed in target cell type: %d target gene-cell pairs passed the >=5%% detection threshold. See results/target_gene_validation_refined.csv.", nrow(answer1)),
  sprintf("2. Differential between MG and HC: %d target gene-cell pairs reached FDR < 0.05 by Wilcoxon test. See results/target_gene_validation_refined.csv.", nrow(answer2)),
  sprintf("3. Direction consistency: %s. Risk genes are expected to be MG-up; protective genes are expected to be MG-down.", paste(paste(answer3$direction_consistent, answer3$n, sep = "="), collapse = "; ")),
  "",
  "## Core refined figures",
  "",
  "- C_refined_celltypes_umap",
  "- D_celltype_proportions_refined",
  "- E_candidate_gene_dotplot_refined",
  "- G_deg_up_down_by_celltype_refined",
  "- H_target_gene_validation_refined"
), report_file)

log_msg("Done refining outputs")
