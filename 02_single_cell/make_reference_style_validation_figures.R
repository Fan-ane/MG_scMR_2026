suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript make_reference_style_validation_figures.R <output_dir>")

output_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
results_dir <- file.path(output_dir, "results")
figures_dir <- file.path(output_dir, "figures")
ref_dir <- file.path(figures_dir, "reference_style")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

save_ref_plot <- function(filename_base, plot, width, height) {
  ggsave(file.path(ref_dir, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 450)
  ggsave(file.path(ref_dir, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, device = grDevices::cairo_pdf)
}

batch_rds <- file.path(results_dir, "seurat_440_validation_refined_batch_corrected.rds")
base_rds <- file.path(results_dir, "seurat_440_validation_refined.rds")
obj <- readRDS(if (file.exists(batch_rds)) batch_rds else base_rds)
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
umap_reduction <- if ("umap_harmony" %in% Reductions(obj)) "umap_harmony" else "umap"
cluster_field <- if ("harmony_clusters" %in% colnames(obj[[]])) "harmony_clusters" else "seurat_clusters"
emb <- Embeddings(obj, umap_reduction)

obj$sample_label <- factor(obj$sample,
  levels = c("H1", "H2", "H3", "M0707-01", "M0709-02", "M0709-03", "M210715-01"),
  labels = c("HC1", "HC2", "HC3", "MG1", "MG2", "MG3", "MG4")
)
obj$cell_type_plot <- factor(obj$cell_type_detail_refined,
  levels = c("CD8_ET", "CD8_NC", "CD8_S100B", "CD4_NC", "NK", "Mono_C", "Mono_NC",
             "B_IN", "B_Mem", "B_cell", "Plasma", "Platelet", "DC"),
  labels = c("CD8 ET", "CD8 NC", "CD8 S100B", "CD4 T", "NK", "Mono C", "Mono NC",
             "B IN", "B Mem", "B cell", "Plasma", "Platelet", "DC")
)
major_map <- c(
  CD8_ET = "CD8_T", CD8_NC = "CD8_T", CD8_S100B = "CD8_T",
  CD4_NC = "CD4_T", NK = "NK", Mono_C = "Mono C", Mono_NC = "Mono NC",
  B_IN = "B IN", B_Mem = "B Mem", B_cell = "B IN", Plasma = "Plasma",
  Platelet = "Platelet", DC = "DC"
)
obj$cell_type_major_plot <- factor(
  unname(major_map[as.character(obj$cell_type_detail_refined)]),
  levels = c("CD8_T", "CD4_T", "NK", "Mono C", "Mono NC", "B IN", "Platelet", "B Mem", "Plasma", "DC")
)

sample_cols <- c(
  HC1 = "#1F78B4", HC2 = "#A6CEE3", HC3 = "#33A02C",
  MG1 = "#FF7F00", MG2 = "#FDBF6F", MG3 = "#31A354", MG4 = "#90D889"
)
cell_cols <- c(
  "CD8_T" = "#1F78B4",
  "CD4_T" = "#A6CEE3",
  "NK" = "#FF7F00",
  "Mono C" = "#FDBF6F",
  "Mono NC" = "#33A02C",
  "B IN" = "#90D889",
  "Platelet" = "#E31A1C",
  "B Mem" = "#FB9A99",
  "Plasma" = "#CAB2D6",
  "DC" = "#9467BD"
)
cluster_raw <- as.character(obj[[cluster_field]][, 1])
cluster_levels <- as.character(sort(unique(as.integer(cluster_raw))))
cluster_values <- factor(cluster_raw, levels = cluster_levels)
cluster_cols <- grDevices::hcl.colors(length(levels(cluster_values)), palette = "Dark 3")
names(cluster_cols) <- levels(cluster_values)

umap_df <- data.frame(
  UMAP1 = emb[, 1],
  UMAP2 = emb[, 2],
  sample_label = obj$sample_label,
  cluster = cluster_values,
  cell_type_plot = obj$cell_type_plot,
  cell_type_major_plot = obj$cell_type_major_plot,
  stringsAsFactors = FALSE
)
set.seed(440)
umap_df <- umap_df[sample(seq_len(nrow(umap_df))), ]

label_centers <- function(var) {
  aggregate(umap_df[, c("UMAP1", "UMAP2")], list(label = umap_df[[var]]), median)
}

theme_ref_umap <- theme_classic(base_size = 11) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.55),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.title = element_text(size = 11),
    axis.text = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(4, 4, 4, 4)
  )

panel_tag <- function(label) {
  plot_annotation(tag_levels = list(label)) &
    theme(plot.tag = element_text(face = "bold", size = 22),
          plot.tag.position = c(0, 1))
}

pA <- ggplot(umap_df, aes(UMAP1, UMAP2, color = sample_label)) +
  geom_point(size = 0.18, alpha = 1, stroke = 0) +
  scale_color_manual(values = sample_cols, name = "Sample", drop = FALSE) +
  coord_equal() +
  labs(x = "UMAP1", y = "UMAP2") +
  annotate("text", x = min(emb[, 1]) + 0.8, y = min(emb[, 2]) + 0.8,
           label = "7 Samples", fontface = "bold.italic", hjust = 0, size = 3.4) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_ref_umap

pB <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.16, alpha = 1, stroke = 0) +
  geom_text(data = label_centers("cluster"), aes(x = UMAP1, y = UMAP2, label = label), inherit.aes = FALSE,
            size = 2.3, color = "black") +
  scale_color_manual(values = cluster_cols, name = "Cluster", drop = FALSE) +
  coord_equal() +
  labs(x = "UMAP1", y = "UMAP2") +
  annotate("text", x = min(emb[, 1]) + 0.8, y = min(emb[, 2]) + 0.8,
           label = paste0(format(ncol(obj), big.mark = ","), " Cells"),
           fontface = "bold.italic", hjust = 0, size = 3.4) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_ref_umap

pC <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cell_type_major_plot)) +
  geom_point(size = 0.18, alpha = 1, stroke = 0) +
  geom_text(data = label_centers("cell_type_major_plot"), aes(x = UMAP1, y = UMAP2, label = label), inherit.aes = FALSE,
            size = 2.4, color = "black") +
  scale_color_manual(values = cell_cols, name = "Cell type", drop = FALSE) +
  coord_equal() +
  labs(x = "UMAP1", y = "UMAP2") +
  annotate("text", x = min(emb[, 1]) + 0.8, y = min(emb[, 2]) + 0.8,
           label = paste0(length(unique(obj$cell_type_major_plot)), " Cell types"),
           fontface = "bold.italic", hjust = 0, size = 3.4) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_ref_umap

save_ref_plot("A_sample_umap_reference_style", pA, width = 5.4, height = 4.6)
save_ref_plot("B_cluster_umap_reference_style", pB, width = 6.6, height = 4.8)
save_ref_plot("C_celltype_umap_reference_style", pC, width = 6.2, height = 4.8)

prop_group <- as.data.frame(table(obj$group, obj$cell_type_major_plot), stringsAsFactors = FALSE)
colnames(prop_group) <- c("group", "cell_type", "cells")
prop_group$cell_type <- factor(prop_group$cell_type, levels = levels(obj$cell_type_major_plot))
prop_group <- prop_group %>%
  group_by(group) %>%
  mutate(prop = cells / sum(cells), ymax = cumsum(prop), ymin = lag(ymax, default = 0),
         label_pos = (ymax + ymin) / 2,
         label = ifelse(prop >= 0.01, paste0(round(prop * 100), "%"), ""))

donut_one <- function(g) {
  dat <- prop_group[prop_group$group == g, ]
  total <- sum(dat$cells)
  ggplot(dat, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.45, fill = cell_type)) +
    geom_rect(color = "white", linewidth = 0.45) +
    geom_text(aes(x = 4.25, y = label_pos, label = label), size = 3) +
    coord_polar(theta = "y") +
    xlim(c(0.6, 4.5)) +
    scale_fill_manual(values = cell_cols, drop = FALSE) +
    annotate("text", x = 0.6, y = 0.5, label = paste0(g, "\n", format(total, big.mark = ","), " Cells"),
             fontface = "bold.italic", size = 5) +
    theme_void(base_size = 11) +
    theme(legend.position = "none")
}

pD_hc <- donut_one("HC")
pD_mg <- donut_one("MG")
pD_leg <- ggplot(prop_group, aes(cell_type, cells, fill = cell_type)) +
  geom_col(alpha = 0) +
  scale_fill_manual(values = cell_cols, drop = FALSE, name = "Cell type") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))
pD <- pD_hc + pD_mg + pD_leg + plot_layout(widths = c(1, 1, 0.55), guides = "collect")
save_ref_plot("D_celltype_donut_reference_style", pD, width = 10.6, height = 4.6)

dot_genes <- intersect(c("HCG23", "ZSCAN12", "AGER", "ZNF322", "SESN3", "TNFRSF14",
                         "IL12RB2", "ORMDL3", "HFE", "CTSH", "CFD", "TYMP"), rownames(obj))
pE <- DotPlot(obj, features = dot_genes, group.by = "cell_type_major_plot", cols = c("#2C7BB6", "#F7F7F7", "#D7191C")) +
  RotatedAxis() +
  scale_color_gradient2(low = "#2C7BB6", mid = "#F7F7F7", high = "#D7191C",
                        midpoint = 0, name = "Relative\nExpression") +
  scale_size(range = c(0.5, 7), name = "Percentage\nExpression") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right"
  )
save_ref_plot("E_candidate_dotplot_reference_style", pE, width = 8.6, height = 5.2)

feature_genes <- intersect(c("AGER", "ZSCAN12", "SESN3", "CFD", "HCG23", "TYMP"), rownames(obj))
make_feature_plot <- function(gene) {
  gene_expr <- as.numeric(expr[gene, colnames(obj)])
  max_cutoff <- if (any(gene_expr > 0)) {
    unname(stats::quantile(gene_expr[gene_expr > 0], 0.98))
  } else {
    1
  }
  max_cutoff <- max(max_cutoff, 1)
  dat <- data.frame(
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    expr = pmin(gene_expr, max_cutoff)
  )
  dat_zero <- dat[dat$expr <= 0, ]
  dat_pos <- dat[dat$expr > 0, ]
  dat_pos <- dat_pos[order(dat_pos$expr), ]
  legend_breaks <- pretty(c(0, max_cutoff), n = 4)
  legend_breaks <- legend_breaks[legend_breaks >= 0 & legend_breaks <= max_cutoff]
  if (!0 %in% legend_breaks) legend_breaks <- c(0, legend_breaks)

  ggplot() +
    geom_point(data = dat_zero, aes(UMAP1, UMAP2), color = "#EFEFEF", size = 0.22, alpha = 1, stroke = 0) +
    geom_point(data = dat_pos, aes(UMAP1, UMAP2, color = expr), size = 0.3, alpha = 1, stroke = 0) +
    scale_color_gradientn(
      colors = c("#EFEFEF", "#FEE5D9", "#FC9272", "#DE2D26"),
      limits = c(0, max_cutoff),
      breaks = legend_breaks,
      oob = scales::squish,
      name = NULL
    ) +
    coord_equal() +
    labs(x = "UMAP1", y = "UMAP2") +
    ggtitle(gene) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 13, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.35, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}
feature_plots <- lapply(feature_genes, function(gene) {
  make_feature_plot(gene)
})
for (i in seq_along(feature_plots)) {
  save_ref_plot(paste0("F_featureplot_", feature_genes[i], "_reference_style"),
                feature_plots[[i]], width = 4.3, height = 3.9)
}
pF <- wrap_plots(feature_plots, ncol = 2)
save_ref_plot("F_key_featureplots_reference_style", pF, width = 6.2, height = 9.8)

deg_summary <- read.csv(file.path(results_dir, "celltype_deg_up_down_summary_refined.csv"), stringsAsFactors = FALSE)
deg_summary$cell_type_major <- unname(major_map[deg_summary$cell_type_detail])
deg_summary <- deg_summary %>%
  group_by(cell_type_major) %>%
  summarise(sig_genes = sum(sig_genes), MG_up = sum(MG_up), MG_down = sum(MG_down), .groups = "drop")
deg_summary$cell_type <- factor(deg_summary$cell_type_major,
  levels = rev(c("CD8_T", "CD4_T", "NK", "Mono C", "Mono NC", "B IN", "Platelet", "B Mem", "Plasma", "DC"))
)
deg_long <- rbind(
  data.frame(cell_type = deg_summary$cell_type, group = "Up", prop = deg_summary$MG_up / pmax(deg_summary$sig_genes, 1)),
  data.frame(cell_type = deg_summary$cell_type, group = "Down", prop = deg_summary$MG_down / pmax(deg_summary$sig_genes, 1))
)
pG <- ggplot(deg_long, aes(x = cell_type, y = prop, fill = group)) +
  geom_col(width = 0.88, color = "white", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c(Up = "#E31A1C", Down = "#1F78B4"), name = "Group") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = 11) +
  theme(legend.title = element_text(face = "bold"), axis.text.y = element_text(size = 10))
save_ref_plot("G_deg_up_down_reference_style", pG, width = 4.8, height = 5.2)

target <- read.csv(file.path(results_dir, "target_gene_validation_refined.csv"), stringsAsFactors = FALSE)
target <- target[!is.na(target$avg_log2FC_MG_vs_HC), ]
target$plot_label <- paste(target$gene, target$target_cell_type, sep = " / ")
target$plot_label <- factor(target$plot_label, levels = rev(target$plot_label[order(target$target_cell_type, target$gene)]))
target$consistency <- ifelse(target$direction_consistent, "Consistent", "Opposite")
target$fdr_label <- ifelse(is.na(target$p_val_adj), "",
                           ifelse(target$p_val_adj < 0.001, "p<0.001",
                                  paste0("p=", signif(target$p_val_adj, 2))))
x_range <- range(target$avg_log2FC_MG_vs_HC, na.rm = TRUE)
x_pad <- max(0.25, diff(x_range) * 0.28)
target$label_x <- target$avg_log2FC_MG_vs_HC +
  ifelse(target$avg_log2FC_MG_vs_HC >= 0, x_pad * 0.18, -x_pad * 0.18)
target$label_hjust <- ifelse(target$avg_log2FC_MG_vs_HC >= 0, 0, 1)
pH <- ggplot(target, aes(y = plot_label, x = avg_log2FC_MG_vs_HC)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey35") +
  geom_segment(aes(x = 0, xend = avg_log2FC_MG_vs_HC, yend = plot_label), color = "grey35", linewidth = 0.35) +
  geom_point(aes(fill = consistency, size = -log10(pmax(p_val_adj, 1e-300))), shape = 21, color = "black", alpha = 0.9) +
  geom_text(aes(x = label_x, label = fdr_label, hjust = label_hjust),
            size = 2.6, color = ifelse(target$p_val_adj < 0.05, "#D7191C", "black")) +
  scale_fill_manual(values = c(Consistent = "#4DBBD5", Opposite = "#E64B35")) +
  scale_size(range = c(2, 8), name = "-log10 FDR") +
  coord_cartesian(xlim = c(x_range[1] - x_pad, x_range[2] + x_pad), clip = "off") +
  labs(x = "log2 Fold Change", y = NULL, fill = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        legend.position = "right",
        plot.margin = margin(5, 25, 5, 5))
save_ref_plot("H_target_validation_reference_style", pH, width = 7, height = 6.2)

message("Reference-style single-panel figures were written to: ", ref_dir)
