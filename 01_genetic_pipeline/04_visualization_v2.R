#!/usr/bin/env Rscript
# =============================================================================
# Step 4 (v2): Revised Visualization
# Changes from v1:
#   - Fig 2: 3-panel nSNP distribution (nSNP=1/2/3), reference paper style
#   - Fig 3: Manhattan plot, unified palette, no title in figure
#   - Fig 4: Forest plot, red/blue for risk/protective, compact layout
#   - Fig 5: PheWAS gene heatmap, unified style, no title
#   - Fig 6: SNP-level PheWAS, unified style, no title
# =============================================================================

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)

RESULTS_DIR <- "results"
FIG_DIR <- "results/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ct_order <- c("B IN", "B Mem", "CD4 ET", "CD4 NC", "CD4 SOX4",
              "CD8 ET", "CD8 NC", "CD8 S100B", "DC",
              "Mono C", "Mono NC", "NK", "NK R", "Plasma")

# Unified cell type colors (consistent muted academic palette)
ct_colors <- c(
  "B IN" = "#F781BF", "B Mem" = "#999999",
  "CD4 ET" = "#FF7F00", "CD4 NC" = "#E41A1C", "CD4 SOX4" = "#FFFF33",
  "CD8 ET" = "#4DAF4A", "CD8 NC" = "#377EB8", "CD8 S100B" = "#00CED1",
  "DC" = "#8DA0CB",
  "Mono C" = "#66C2A5", "Mono NC" = "#FC8D62",
  "NK" = "#984EA3", "NK R" = "#A65628",
  "Plasma" = "#E78AC3"
)

# Shared theme: clean academic style, no titles in image
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )
}

# Sample sizes (OneK1K)
sample_sizes <- data.table(
  cell_type = ct_order,
  n_samples = c(982L, 982L, 982L, 982L, 857L,
                982L, 982L, 981L, 968L,
                969L, 932L, 928L, 967L, 643L)
)

# Load data from CSVs
mr_all <- fread(file.path(RESULTS_DIR, "Supplementary_Table_S2_MR_all.csv"))
mr_sig <- fread(file.path(RESULTS_DIR, "Supplementary_Table_S2_MR_significant.csv"))
coloc_merged <- fread(file.path(RESULTS_DIR, "Supplementary_Table_Coloc_merged.csv"))
phewas_gene <- fread(file.path(RESULTS_DIR, "Supplementary_Table_PheWAS_gene.csv"))
phewas_snp <- fread(file.path(RESULTS_DIR, "Supplementary_Table_PheWAS_snp.csv"))

# =============================================================================
# Fig 2: Instrument nSNP Distribution (3-panel, reference style)
# =============================================================================
cat("Generating Fig 2...\n")

nsnp_data <- mr_all[, .(nSNP = nsnp[1]), by = .(gene, cell_type)]
nsnp_summary <- nsnp_data[, .(Count = .N), by = .(cell_type, nSNP)]

all_combos <- CJ(cell_type = ct_order, nSNP = 1:3)
nsnp_summary <- merge(all_combos, nsnp_summary, by = c("cell_type", "nSNP"), all.x = TRUE)
nsnp_summary[is.na(Count), Count := 0L]
nsnp_summary[, cell_type := factor(cell_type, levels = ct_order)]

make_nsnp_panel <- function(dt, n, fill_color) {
  sub <- dt[nSNP == n]
  max_count <- max(sub$Count)
  p <- ggplot(sub, aes(x = cell_type, y = Count)) +
    geom_col(fill = fill_color, width = 0.7) +
    geom_text(aes(label = Count), vjust = -0.3, size = 3.2) +
    labs(x = NULL, y = "Count", subtitle = paste0("nSNP = ", n)) +
    theme_pub() +
    theme(
      plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )
  if (max_count == 0) {
    p <- p +
      scale_y_continuous(limits = c(0, 5), breaks = 0:5) +
      annotate("text", x = 7.5, y = 2.5,
               label = "n = 0 for all cell types",
               color = "grey50", size = 4, fontface = "italic")
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  }
  p
}

p2a <- make_nsnp_panel(nsnp_summary, 1, "#F8766D")
p2b <- make_nsnp_panel(nsnp_summary, 2, "#619CFF")
p2c <- make_nsnp_panel(nsnp_summary, 3, "#00BA38")

fig2 <- p2a / p2b / p2c + plot_layout(heights = c(3, 1, 1))
ggsave(file.path(FIG_DIR, "Fig2_eGene_distribution.pdf"), fig2, width = 10, height = 12)
ggsave(file.path(FIG_DIR, "Fig2_eGene_distribution.png"), fig2, width = 10, height = 12, dpi = 300)
cat("  Fig 2 done\n")

# =============================================================================
# Fig 3: Manhattan plot (by cell type, unified colors, no title)
# =============================================================================
cat("Generating Fig 3...\n")

mr_all[, cell_type_f := factor(cell_type, levels = ct_order)]
mr_all[, neglogp := -log10(pval)]

fdr_threshold <- -log10(max(mr_sig$pval))

top_genes <- mr_all[fdr < 0.05][order(pval)][1:min(50, sum(mr_all$fdr < 0.05, na.rm = TRUE))]

p3 <- ggplot(mr_all, aes(x = cell_type_f, y = neglogp, color = cell_type_f)) +
  geom_jitter(width = 0.3, alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = 14, y = fdr_threshold + 0.3, label = "FDR = 0.05",
           color = "red", hjust = 1, size = 3) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 2.5,
                  max.overlaps = 30, segment.size = 0.2) +
  scale_color_manual(values = ct_colors) +
  labs(x = NULL, y = expression(-log[10](P))) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(FIG_DIR, "Fig3_Manhattan_MR.pdf"), p3, width = 14, height = 7)
ggsave(file.path(FIG_DIR, "Fig3_Manhattan_MR.png"), p3, width = 14, height = 7, dpi = 300)
cat("  Fig 3 done\n")

# =============================================================================
# Fig 4: Forest plot (compact, red/blue, reference style with table)
# =============================================================================
cat("Generating Fig 4...\n")

tier1 <- coloc_merged[tier == "Tier 1 (Strong)"]

if (nrow(tier1) > 30) {
  tier1 <- tier1[order(-abs(log(or)))][1:30]
}

tier1[, direction := ifelse(b < 0, "Protective", "Risk")]
tier1[, label := paste0(gene, " (", cell_type, ")")]

# Get MR method from MR sig data (nsnps in coloc_merged is coloc overlap, not MR instruments)
mr_method <- mr_sig[, .(gene, cell_type, mr_method = method, mr_nsnp = nsnp)]
tier1 <- merge(tier1, mr_method, by = c("gene", "cell_type"), all.x = TRUE)
tier1[is.na(mr_method), mr_method := "Wald ratio"]
tier1 <- tier1[order(or)]
tier1[, label := factor(label, levels = label)]
tier1[, or_ci_text := sprintf("%.2f (%.2f-%.2f)", or, or_lci, or_uci)]
tier1[, pp_text := sprintf("%.1f%%", PP.H4 * 100)]

n_genes <- nrow(tier1)
plot_height <- max(6, n_genes * 0.25 + 1.5)

# Left table panel: cell_type, gene, method, FDR
p4_table <- ggplot(tier1, aes(y = label)) +
  geom_text(aes(x = 0, label = cell_type), size = 2.2, hjust = 0) +
  geom_text(aes(x = 1, label = gene), size = 2.2, hjust = 0, fontface = "italic") +
  geom_text(aes(x = 2, label = mr_method), size = 2.0, hjust = 0) +
  geom_text(aes(x = 3, label = sprintf("%.1e", fdr)), size = 2.0, hjust = 0) +
  annotate("text", x = 0, y = n_genes + 0.8, label = "Cell type",
           size = 2.5, fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = n_genes + 0.8, label = "eGene",
           size = 2.5, fontface = "bold", hjust = 0) +
  annotate("text", x = 2, y = n_genes + 0.8, label = "Method",
           size = 2.5, fontface = "bold", hjust = 0) +
  annotate("text", x = 3, y = n_genes + 0.8, label = "FDR",
           size = 2.5, fontface = "bold", hjust = 0) +
  scale_x_continuous(limits = c(-0.1, 4)) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

# Center: forest plot
p4_forest <- ggplot(tier1, aes(x = or, y = label, color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5) +
  geom_errorbar(aes(xmin = or_lci, xmax = or_uci), width = 0.25,
                orientation = "y") +
  scale_color_manual(values = c("Protective" = "#2166AC", "Risk" = "#B2182B"),
                     name = NULL) +
  scale_x_log10() +
  labs(x = "Odds Ratio (95% CI, log scale)", y = NULL) +
  annotate("text", x = 0.15, y = 0.3, label = "Protective",
           color = "#2166AC", fontface = "bold", size = 3.5) +
  annotate("text", x = 8, y = 0.3, label = "Risk",
           color = "#B2182B", fontface = "bold", size = 3.5) +
  coord_cartesian(clip = "off") +
  theme_pub() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey95"))

# Right panel: OR (95% CI) and PP.H4
p4_right <- ggplot(tier1, aes(y = label)) +
  geom_text(aes(x = 0, label = or_ci_text), size = 2.0, hjust = 0) +
  geom_text(aes(x = 1, label = pp_text), size = 2.0, hjust = 0) +
  annotate("text", x = 0, y = n_genes + 0.8, label = "OR (95% CI)",
           size = 2.5, fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = n_genes + 0.8, label = "PP.H4",
           size = 2.5, fontface = "bold", hjust = 0) +
  scale_x_continuous(limits = c(-0.1, 1.5)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 0))

fig4 <- p4_table + p4_forest + p4_right + plot_layout(widths = c(2.5, 3, 1.5))
ggsave(file.path(FIG_DIR, "Fig4_Coloc_Forest.pdf"), fig4, width = 16, height = plot_height)
ggsave(file.path(FIG_DIR, "Fig4_Coloc_Forest.png"), fig4, width = 16, height = plot_height, dpi = 300)
cat("  Fig 4 done\n")

# =============================================================================
# Fig 5: Gene-level PheWAS heatmap (no title, unified style)
# =============================================================================
cat("Generating Fig 5...\n")

phewas_plot <- phewas_gene[overall_score > 0.3]
phewas_plot <- phewas_plot[, .SD[order(-overall_score)][1:min(5, .N)], by = gene]

if (nrow(phewas_plot) > 0) {
  p5_height <- min(40, max(8, uniqueN(phewas_plot$disease) * 0.35))
  p5_width <- min(40, max(10, uniqueN(phewas_plot$gene) * 0.5))

  p5 <- ggplot(phewas_plot, aes(x = gene, y = disease, fill = overall_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_grid(therapeutic_area ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "lightyellow", high = "darkorange",
                        name = "Association\nScore") +
    labs(x = "Colocalized eGene", y = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
          axis.text.y = element_text(size = 6),
          strip.text.y = element_text(angle = 0, size = 7))

  ggsave(file.path(FIG_DIR, "Fig5_PheWAS_gene.pdf"), p5, width = p5_width, height = p5_height)
  ggsave(file.path(FIG_DIR, "Fig5_PheWAS_gene.png"), p5, width = p5_width, height = p5_height, dpi = 300)
}
cat("  Fig 5 done\n")

# =============================================================================
# Fig 6: SNP-level PheWAS (no title, unified style)
# =============================================================================
cat("Generating Fig 6...\n")

# Reconstruct trait_category from trait text
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

phewas_snp[, neglogp := -log10(as.numeric(pval))]
phewas_snp[, trait_category := sapply(trait, categorize_trait)]
phewas_snp[, trait_category_f := factor(trait_category)]

top_snp_assoc <- phewas_snp[order(-neglogp)][1:min(30, nrow(phewas_snp))]

p6 <- ggplot(phewas_snp, aes(x = trait_category_f, y = neglogp, color = trait_category_f)) +
  geom_jitter(width = 0.3, alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  geom_text_repel(data = top_snp_assoc,
                  aes(label = paste(gene, trait, sep = "\n")),
                  size = 2, max.overlaps = 20, segment.size = 0.2) +
  labs(x = NULL, y = expression(-log[10](P))) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(title = "Trait Category"))

ggsave(file.path(FIG_DIR, "Fig6_PheWAS_snp.pdf"), p6, width = 14, height = 8)
ggsave(file.path(FIG_DIR, "Fig6_PheWAS_snp.png"), p6, width = 14, height = 8, dpi = 300)
cat("  Fig 6 done\n")

cat("All 5 figures generated.\n")
