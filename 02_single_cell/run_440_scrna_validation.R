suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_440_scrna_validation.R <input_dir> <output_dir>")
}

input_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_dir <- file.path(output_dir, "raw_10x")
results_dir <- file.path(output_dir, "results")
figures_dir <- file.path(output_dir, "figures")
logs_dir <- file.path(output_dir, "logs")
for (d in c(raw_dir, results_dir, figures_dir, logs_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

set.seed(440)
theme_set(theme_bw(base_size = 10))

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
}

zip_files <- list.files(input_dir, pattern = "matrix.*\\.zip$", full.names = TRUE)
if (length(zip_files) == 0) stop("No matrix zip files found in: ", input_dir)

sample_from_zip <- function(path) {
  x <- basename(path)
  x <- sub("_matrix_10X\\.zip$", "", x)
  x <- sub("_matrix\\.zip$", "", x)
  sub("\\.zip$", "", x)
}

samples <- data.frame(
  sample = vapply(zip_files, sample_from_zip, character(1)),
  zip = zip_files,
  stringsAsFactors = FALSE
)
samples$group <- ifelse(grepl("^H", samples$sample, ignore.case = TRUE), "HC", "MG")
samples <- samples[order(samples$group, samples$sample), ]
write.csv(samples, file.path(results_dir, "sample_metadata.csv"), row.names = FALSE)

find_10x_dir <- function(parent) {
  mtx <- list.files(parent, pattern = "^matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(mtx) == 0) stop("No matrix.mtx found under: ", parent)
  dirname(mtx[1])
}

read_one_sample <- function(sample, zip, group) {
  out_dir <- file.path(raw_dir, sample)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(out_dir, pattern = "^matrix\\.mtx(\\.gz)?$", recursive = TRUE)) == 0) {
    log_msg("Unzipping", sample)
    utils::unzip(zip, exdir = out_dir)
  }
  data_dir <- find_10x_dir(out_dir)
  log_msg("Reading", sample, "from", data_dir)
  counts <- Read10X(data.dir = data_dir, gene.column = 2, unique.features = TRUE)
  if (is.list(counts)) {
    counts <- counts[[grep("Gene Expression|RNA", names(counts), ignore.case = TRUE)[1]]]
  }
  obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
  obj$sample <- sample
  obj$group <- group
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  before <- ncol(obj)
  obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 20)
  after <- ncol(obj)
  attr(obj, "qc") <- data.frame(sample = sample, group = group, cells_before = before, cells_after = after)
  obj
}

objs <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  objs[[i]] <- read_one_sample(samples$sample[i], samples$zip[i], samples$group[i])
}
names(objs) <- samples$sample
qc <- bind_rows(lapply(objs, attr, "qc"))
write.csv(qc, file.path(results_dir, "sample_qc.csv"), row.names = FALSE)

log_msg("Merging samples")
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "MG_HC_440")
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj <- JoinLayers(obj)
}
obj$group <- factor(obj$group, levels = c("HC", "MG"))

log_msg("Normalizing and reducing dimensions")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

marker_sets <- list(
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "HLA-DRA"),
  B_IN = c("IGHD", "IGHM", "TCL1A", "FCER2", "IL4R", "CD24", "CD38"),
  B_Mem = c("CD27", "TNFRSF13B", "CD80", "CD86", "FCRL5", "IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGHA2"),
  Plasma = c("MZB1", "XBP1", "JCHAIN", "SDC1", "PRDM1", "TNFRSF17"),
  CD4_NC = c("CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LTB", "TCF7"),
  CD8_NC = c("CD3D", "CD3E", "CD8A", "CD8B", "CCR7", "LEF1", "TCF7"),
  CD8_ET = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "GNLY", "GZMH"),
  CD8_S100B = c("CD3D", "CD3E", "CD8A", "CD8B", "S100B", "GZMK"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "NCAM1", "PRF1"),
  Mono_C = c("LYZ", "S100A8", "S100A9", "FCN1", "VCAN", "CD14", "LST1"),
  Mono_NC = c("LYZ", "FCGR3A", "MS4A7", "LST1", "LILRB2", "IFITM3"),
  DC = c("FCER1A", "CLEC10A", "CST3", "HLA-DRA", "ITGAX", "LILRA4", "GZMB"),
  Platelet = c("PPBP", "PF4", "GP9", "SDPR")
)

expr_mat <- GetAssayData(obj, assay = "RNA", slot = "data")
clusters <- levels(obj$seurat_clusters)
score_rows <- list()
for (cl in clusters) {
  cells <- WhichCells(obj, idents = cl)
  for (nm in names(marker_sets)) {
    genes <- intersect(marker_sets[[nm]], rownames(expr_mat))
    score <- if (length(genes) == 0) NA_real_ else mean(Matrix::rowMeans(expr_mat[genes, cells, drop = FALSE]))
    score_rows[[length(score_rows) + 1]] <- data.frame(cluster = cl, cell_type = nm, score = score)
  }
}
cluster_scores <- bind_rows(score_rows)
cluster_anno <- cluster_scores %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(cell_type_detail = ifelse(is.na(score) | score <= 0, "Other", cell_type))
write.csv(cluster_scores, file.path(results_dir, "cluster_marker_scores.csv"), row.names = FALSE)
write.csv(cluster_anno, file.path(results_dir, "cluster_annotation.csv"), row.names = FALSE)

anno_map <- setNames(cluster_anno$cell_type_detail, cluster_anno$cluster)
obj$cell_type_detail <- unname(anno_map[as.character(obj$seurat_clusters)])
obj$cell_type_detail[is.na(obj$cell_type_detail)] <- "Other"
obj$cell_type_detail <- factor(obj$cell_type_detail)

prop <- as.data.frame(table(obj$sample, obj$group, obj$cell_type_detail), stringsAsFactors = FALSE)
colnames(prop) <- c("sample", "group", "cell_type_detail", "cells")
prop <- prop %>%
  group_by(sample) %>%
  mutate(total_cells = sum(cells), proportion = cells / total_cells) %>%
  ungroup()
write.csv(prop, file.path(results_dir, "celltype_proportions_by_sample.csv"), row.names = FALSE)

tier1_genes <- c("HCG23", "ZSCAN12", "AGER", "SCAND3", "ZNF322",
                 "SESN3", "TNFRSF14", "IL12RB2", "ORMDL3", "HFE", "CTSH")
drug_genes <- c("AGER", "CFD", "TYMP")
candidate_genes <- unique(c(tier1_genes, drug_genes))
candidate_present <- intersect(candidate_genes, rownames(obj))
candidate_missing <- setdiff(candidate_genes, rownames(obj))
writeLines(candidate_missing, file.path(results_dir, "candidate_genes_missing.txt"))

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
  cells <- WhichCells(obj, expression = cell_type_detail == ct)
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
write.csv(target_validation, file.path(results_dir, "target_gene_validation.csv"), row.names = FALSE)

log_msg("Running cell-type-level DEG summaries")
deg_summary <- list()
cell_types <- sort(unique(as.character(obj$cell_type_detail)))
for (ct in cell_types) {
  cells <- WhichCells(obj, expression = cell_type_detail == ct)
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
  fm$direction <- ifelse(fm[[fc_col]] > 0, "MG_up", "MG_down")
  sig <- fm[!is.na(fm$p_val_adj) & fm$p_val_adj < 0.05, , drop = FALSE]
  write.csv(fm, file.path(results_dir, paste0("DEG_", ct, "_MG_vs_HC.csv")), row.names = FALSE)
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
write.csv(deg_summary, file.path(results_dir, "celltype_deg_up_down_summary.csv"), row.names = FALSE)

saveRDS(obj, file.path(results_dir, "seurat_440_validation.rds"))

log_msg("Plotting")
pA1 <- DimPlot(obj, reduction = "umap", group.by = "group", pt.size = 0.15) + ggtitle("A1 Group")
pA2 <- DimPlot(obj, reduction = "umap", group.by = "sample", pt.size = 0.15) + ggtitle("A2 Sample")
pB <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.15) + ggtitle("B Clusters")
pC <- DimPlot(obj, reduction = "umap", group.by = "cell_type_detail", label = TRUE, repel = TRUE, pt.size = 0.15) + ggtitle("C Cell types")
ggsave(file.path(figures_dir, "A_B_C_umap_overview.pdf"), (pA1 + pA2) / (pB + pC), width = 13, height = 10)
ggsave(file.path(figures_dir, "A_B_C_umap_overview.png"), (pA1 + pA2) / (pB + pC), width = 13, height = 10, dpi = 300)

pD <- ggplot(prop, aes(x = group, y = proportion, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_point(aes(shape = sample), position = position_jitter(width = 0.12, height = 0), size = 1.7) +
  facet_wrap(~ cell_type_detail, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell proportion", title = "D Cell-type proportions") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(figures_dir, "D_celltype_proportions.pdf"), pD, width = 12, height = 9)
ggsave(file.path(figures_dir, "D_celltype_proportions.png"), pD, width = 12, height = 9, dpi = 300)

if (length(candidate_present) > 0) {
  pE <- DotPlot(obj, features = candidate_present, group.by = "cell_type_detail") +
    RotatedAxis() +
    labs(title = "E Candidate genes across annotated immune cell types")
  ggsave(file.path(figures_dir, "E_candidate_gene_dotplot.pdf"), pE, width = max(8, length(candidate_present) * 0.55), height = 5.5)
  ggsave(file.path(figures_dir, "E_candidate_gene_dotplot.png"), pE, width = max(8, length(candidate_present) * 0.55), height = 5.5, dpi = 300)
}

feature_genes <- intersect(c("AGER", "ZSCAN12", "SESN3", "CFD", "HCG23"), rownames(obj))
if (length(feature_genes) > 0) {
  pF <- FeaturePlot(obj, features = feature_genes, reduction = "umap", pt.size = 0.12,
                    order = TRUE, ncol = min(3, length(feature_genes))) &
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(figures_dir, "F_key_gene_featureplots.pdf"), pF, width = 11, height = ceiling(length(feature_genes) / 3) * 3.2)
  ggsave(file.path(figures_dir, "F_key_gene_featureplots.png"), pF, width = 11, height = ceiling(length(feature_genes) / 3) * 3.2, dpi = 300)
}

if (nrow(deg_summary) > 0) {
  deg_long <- rbind(
    data.frame(cell_type_detail = deg_summary$cell_type_detail, direction = "MG_up", genes = deg_summary$MG_up),
    data.frame(cell_type_detail = deg_summary$cell_type_detail, direction = "MG_down", genes = -deg_summary$MG_down)
  )
  pG <- ggplot(deg_long, aes(x = reorder(cell_type_detail, genes), y = genes, fill = direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "Significant DEG count (FDR < 0.05)", title = "G Up/down DEG balance by cell type") +
    theme(legend.position = "bottom")
  ggsave(file.path(figures_dir, "G_deg_up_down_by_celltype.pdf"), pG, width = 8, height = 6)
  ggsave(file.path(figures_dir, "G_deg_up_down_by_celltype.png"), pG, width = 8, height = 6, dpi = 300)
}

hdat <- target_validation %>%
  mutate(
    minus_log10_fdr = -log10(pmax(p_val_adj, 1e-300)),
    consistency = ifelse(is.na(direction_consistent), "not_tested",
                         ifelse(direction_consistent, "consistent", "opposite")),
    label = paste(gene, target_cell_type, sep = " / ")
  )
if (nrow(hdat) > 0) {
  pH <- ggplot(hdat, aes(x = gene, y = target_cell_type)) +
    geom_point(aes(size = minus_log10_fdr, color = avg_log2FC_MG_vs_HC, shape = consistency), alpha = 0.9) +
    scale_color_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0,
                          na.value = "grey80", name = "log2FC\nMG vs HC") +
    scale_size_continuous(range = c(1.5, 7), name = "-log10 FDR") +
    labs(x = NULL, y = NULL, title = "H Target-cell validation of candidate genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  ggsave(file.path(figures_dir, "H_target_gene_validation.pdf"), pH, width = 10.5, height = 6.5)
  ggsave(file.path(figures_dir, "H_target_gene_validation.png"), pH, width = 10.5, height = 6.5, dpi = 300)
}

answer1 <- target_validation %>%
  filter(gene_present, expressed_in_target) %>%
  distinct(gene, target_cell_type)
answer2 <- target_validation %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  distinct(gene, target_cell_type, avg_log2FC_MG_vs_HC, p_val_adj)
answer3 <- target_validation %>%
  filter(!is.na(direction_consistent)) %>%
  count(direction_consistent, name = "n")

report_file <- file.path(output_dir, "report_validation_summary.md")
writeLines(c(
  "# 440 single-cell validation summary",
  "",
  sprintf("- Samples: %d HC and %d MG.", sum(samples$group == "HC"), sum(samples$group == "MG")),
  sprintf("- Cells retained after QC: %s.", format(ncol(obj), big.mark = ",")),
  sprintf("- Candidate genes present in matrix: %s.", paste(candidate_present, collapse = ", ")),
  if (length(candidate_missing) > 0) sprintf("- Candidate genes missing from matrix: %s.", paste(candidate_missing, collapse = ", ")) else "- Candidate genes missing from matrix: none.",
  "",
  "## Three document questions",
  "",
  sprintf("1. Expressed in target cell type: %d target gene-cell pairs passed the >=5%% detection threshold. See results/target_gene_validation.csv.", nrow(answer1)),
  sprintf("2. Differential between MG and HC: %d target gene-cell pairs reached FDR < 0.05 by Wilcoxon test. See results/target_gene_validation.csv.", nrow(answer2)),
  sprintf("3. Direction consistency: %s. Risk genes are expected to be MG-up; protective genes are expected to be MG-down.", paste(paste(answer3$direction_consistent, answer3$n, sep = "="), collapse = "; ")),
  "",
  "## Core figures",
  "",
  "- A_B_C_umap_overview: group/sample UMAP, clusters, and cell-type annotation.",
  "- D_celltype_proportions: HC vs MG cell-type proportions by sample.",
  "- E_candidate_gene_dotplot: candidate-gene expression across cell types.",
  "- F_key_gene_featureplots: key candidate FeaturePlots.",
  "- G_deg_up_down_by_celltype: MG-up and MG-down DEG counts by cell type.",
  "- H_target_gene_validation: target cell-type log2FC, FDR, and direction consistency."
), report_file)

log_msg("Done. Outputs written to", output_dir)
