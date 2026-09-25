#!/usr/bin/env Rscript
# =============================================================================
# Step 4 (v3): Publication-quality Visualization — Unified style
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

# --- Unified color palette ---
ct_colors <- c(
  "B IN" = "#F781BF", "B Mem" = "#999999",
  "CD4 ET" = "#FF7F00", "CD4 NC" = "#E41A1C", "CD4 SOX4" = "#FFFF33",
  "CD8 ET" = "#4DAF4A", "CD8 NC" = "#377EB8", "CD8 S100B" = "#00CED1",
  "DC" = "#8DA0CB",
  "Mono C" = "#66C2A5", "Mono NC" = "#FC8D62",
  "NK" = "#984EA3", "NK R" = "#A65628",
  "Plasma" = "#E78AC3"
)

dir_colors <- c("Protective" = "#2166AC", "Risk" = "#B2182B")
tier_fills <- c("Tier 1 (PP.H4>80%)" = "#B2182B",
                "Tier 2 (PP.H4 50-80%)" = "#FC8D62",
                "MR Significant" = "#FDDBC7")

theme_pub <- function(base_sz = 12) {
  theme_minimal(base_size = base_sz) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

# Load data
mr_all <- fread(file.path(RESULTS_DIR, "Supplementary_Table_S2_MR_all.csv"))
mr_sig <- fread(file.path(RESULTS_DIR, "Supplementary_Table_S2_MR_significant.csv"))
coloc_merged <- fread(file.path(RESULTS_DIR, "Supplementary_Table_Coloc_merged.csv"))
phewas_gene <- fread(file.path(RESULTS_DIR, "Supplementary_Table_PheWAS_gene.csv"))
phewas_snp <- fread(file.path(RESULTS_DIR, "Supplementary_Table_PheWAS_snp.csv"))

colocalized <- coloc_merged[PP.H4 > 0.5]

# =============================================================================
# Fig 2: eGene distribution — 3-layer stacked (Sig / Tier2 / Tier1)
# =============================================================================
cat("Generating Fig 2...\n")

# Count per cell type at each level
sig_counts <- mr_sig[, .N, by = cell_type]
setnames(sig_counts, "N", "n_sig")

coloc_counts <- colocalized[, .(
  n_tier1 = sum(tier == "Tier 1 (Strong)"),
  n_tier2 = sum(tier == "Tier 2 (Moderate)")
), by = cell_type]

all_ct <- data.table(cell_type = ct_order)
fig2_data <- merge(all_ct, sig_counts, by = "cell_type", all.x = TRUE)
fig2_data <- merge(fig2_data, coloc_counts, by = "cell_type", all.x = TRUE)
fig2_data[is.na(n_sig), n_sig := 0L]
fig2_data[is.na(n_tier1), n_tier1 := 0L]
fig2_data[is.na(n_tier2), n_tier2 := 0L]
fig2_data[, n_sig_only := n_sig - n_tier1 - n_tier2]
fig2_data[, cell_type := factor(cell_type, levels = rev(ct_order))]

fig2_long <- melt(fig2_data, id.vars = "cell_type",
                  measure.vars = c("n_tier1", "n_tier2", "n_sig_only"),
                  variable.name = "layer", value.name = "count")
fig2_long[, layer := factor(layer,
  levels = c("n_sig_only", "n_tier2", "n_tier1"),
  labels = c("MR Significant", "Tier 2 (PP.H4 50-80%)", "Tier 1 (PP.H4>80%)"))]

p2 <- ggplot(fig2_long, aes(x = count, y = cell_type, fill = layer)) +
  geom_col(width = 0.7) +
  geom_text(data = fig2_data,
            aes(x = n_sig, y = cell_type, label = n_sig, fill = NULL),
            hjust = -0.3, size = 3.8) +
  scale_fill_manual(values = tier_fills, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of eGene–Cell Type Pairs", y = NULL) +
  theme_pub(13) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 11))

ggsave(file.path(FIG_DIR, "Fig2_eGene_distribution.pdf"), p2, width = 9, height = 7)
ggsave(file.path(FIG_DIR, "Fig2_eGene_distribution.png"), p2, width = 9, height = 7, dpi = 300)

# =============================================================================
# Fig 3: Manhattan — colocalized genes, top 15 labels
# =============================================================================
cat("Generating Fig 3...\n")

coloc_ids <- paste(colocalized$gene, colocalized$cell_type)
mr_coloc <- mr_all[paste(gene, cell_type) %in% coloc_ids]
mr_coloc[, cell_type_f := factor(cell_type, levels = ct_order)]
mr_coloc[, neglogp := -log10(pval)]
mr_coloc <- merge(mr_coloc, colocalized[, .(gene, cell_type, PP.H4, tier)],
                  by = c("gene", "cell_type"), all.x = TRUE)

set.seed(42)
mr_coloc[, x_jit := as.numeric(cell_type_f) + runif(.N, -0.3, 0.3)]
top_labels <- mr_coloc[order(pval)][!duplicated(gene)][1:15]

p3 <- ggplot(mr_coloc, aes(x = x_jit, y = neglogp, color = cell_type_f)) +
  geom_point(aes(size = ifelse(tier == "Tier 1 (Strong)", "Tier 1", "Tier 2")),
             alpha = 0.7) +
  scale_size_manual(values = c("Tier 1" = 3.5, "Tier 2" = 1.2), name = NULL) +
  geom_text_repel(data = top_labels, aes(label = gene), size = 3.5,
                  max.overlaps = 15, segment.size = 0.3,
                  fontface = "italic", box.padding = 0.5) +
  scale_color_manual(values = ct_colors, guide = "none") +
  scale_x_continuous(breaks = seq_along(ct_order), labels = ct_order) +
  labs(x = NULL, y = expression(-log[10](P))) +
  theme_pub(13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        legend.position = "top")

ggsave(file.path(FIG_DIR, "Fig3_Manhattan_MR.pdf"), p3, width = 12, height = 7)
ggsave(file.path(FIG_DIR, "Fig3_Manhattan_MR.png"), p3, width = 12, height = 7, dpi = 300)

# =============================================================================
# Fig 4: Forest plot — Tier 1 + Tier 2 mixed, compact
# =============================================================================
cat("Generating Fig 4...\n")

tier1 <- coloc_merged[tier == "Tier 1 (Strong)"]
tier2 <- coloc_merged[tier == "Tier 2 (Moderate)"]

# Top 5 Tier 1 by |log(OR)| + top 5 Tier 2
tier1_top <- tier1[order(-abs(log(or)))][1:min(8, nrow(tier1))]
tier2_top <- tier2[order(-abs(log(or)))][1:min(4, nrow(tier2))]
forest_data <- rbind(tier1_top, tier2_top)
forest_data[, direction := ifelse(b < 0, "Protective", "Risk")]
forest_data[, tier_label := ifelse(tier == "Tier 1 (Strong)", "Tier 1", "Tier 2")]

forest_data[, label := paste0(gene, " (", cell_type, ")")]
forest_data <- forest_data[order(or)]
forest_data[, label := factor(label, levels = label)]
forest_data[, or_text := sprintf("%.2f (%.2f–%.2f)", or, or_lci, or_uci)]
forest_data[, pp_text := sprintf("%.1f%%", PP.H4 * 100)]

n_genes <- nrow(forest_data)
x_or <- max(forest_data$or_uci) * 1.3
x_pp <- max(forest_data$or_uci) * 4.5

p4 <- ggplot(forest_data, aes(x = or, y = label, color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_point(aes(shape = tier_label), size = 4) +
  geom_errorbar(aes(xmin = or_lci, xmax = or_uci), width = 0.3, linewidth = 0.8,
                orientation = "y") +
  geom_text(aes(x = x_or, label = or_text), size = 3.5, hjust = 0, color = "grey20") +
  geom_text(aes(x = x_pp, label = pp_text), size = 3.5, hjust = 0, color = "grey20") +
  annotate("text", x = x_or, y = n_genes + 0.7,
           label = "OR (95% CI)", size = 3.5, fontface = "bold", hjust = 0) +
  annotate("text", x = x_pp, y = n_genes + 0.7,
           label = "PP.H4", size = 3.5, fontface = "bold", hjust = 0) +
  scale_color_manual(values = dir_colors, name = "Direction") +
  scale_shape_manual(values = c("Tier 1" = 16, "Tier 2" = 17), name = "Evidence") +
  scale_x_log10() +
  coord_cartesian(clip = "off") +
  labs(x = "Odds Ratio (95% CI, log scale)", y = NULL) +
  theme_pub(13) +
  theme(axis.text.y = element_text(size = 11),
        legend.position = "top",
        legend.box = "horizontal",
        panel.grid.major.y = element_line(color = "grey95"),
        plot.margin = margin(10, 140, 10, 10))

ggsave(file.path(FIG_DIR, "Fig4_Coloc_Forest.pdf"), p4, width = 13, height = 7)
ggsave(file.path(FIG_DIR, "Fig4_Coloc_Forest.png"), p4, width = 13, height = 7, dpi = 300)

# =============================================================================
# Fig 5: PheWAS gene — lollipop chart
# =============================================================================
cat("Generating Fig 5...\n")

phewas_filt <- phewas_gene[overall_score > 0.5]
phewas_filt <- phewas_filt[gene %in% colocalized$gene]

gene_max <- phewas_filt[, .(max_score = max(overall_score)), by = gene]
top_genes_phewas <- gene_max[order(-max_score)][1:min(15, nrow(gene_max))]$gene
phewas_filt <- phewas_filt[gene %in% top_genes_phewas]
phewas_filt <- phewas_filt[, .SD[order(-overall_score)][1:min(3, .N)], by = gene]

if (nrow(phewas_filt) > 0) {
  phewas_filt[, ta_short := gsub(";.*", "", therapeutic_area)]
  phewas_filt[, ta_short := trimws(ta_short)]

  ta_map <- c(
    "genetic, familial or congenital disease" = "Genetic/Congenital",
    "gastrointestinal disease" = "Gastrointestinal",
    "endocrine system disease" = "Endocrine",
    "nervous system disease" = "Neurological",
    "immune system disease" = "Immune",
    "cancer or benign tumor" = "Neoplasm",
    "cardiovascular disease" = "Cardiovascular",
    "hematologic disease" = "Hematological",
    "urinary system disease" = "Urinary",
    "psychiatric disorder" = "Psychiatric",
    "musculoskeletal or connective tissue disease" = "Musculoskeletal",
    "nutritional or metabolic disease" = "Metabolic",
    "phenotype" = "Other",
    "respiratory or thoracic disease" = "Respiratory"
  )
  phewas_filt[, ta_label := ifelse(ta_short %in% names(ta_map), ta_map[ta_short], "Other")]

  wrap_label <- function(x, width = 35) {
    sapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"))
  }
  phewas_filt[, disease_wrap := wrap_label(disease, 35)]
  phewas_filt <- phewas_filt[order(gene, -overall_score)]
  phewas_filt[, row_label := paste0(disease_wrap, "  (", gene, ")")]
  phewas_filt[, row_label := factor(row_label, levels = rev(row_label))]

  ta_colors <- c(
    "Genetic/Congenital" = "#B2182B", "Gastrointestinal" = "#FF7F00",
    "Endocrine" = "#FFCC00", "Neurological" = "#4DAF4A",
    "Immune" = "#2166AC", "Neoplasm" = "#984EA3",
    "Cardiovascular" = "#E78AC3", "Hematological" = "#A65628",
    "Urinary" = "#66C2A5", "Psychiatric" = "#8DA0CB",
    "Metabolic" = "#FC8D62", "Musculoskeletal" = "#B3B3B3",
    "Respiratory" = "#00CED1", "Other" = "#999999"
  )

  n_rows <- nrow(phewas_filt)
  p5_height <- max(7, n_rows * 0.28 + 2)

  p5 <- ggplot(phewas_filt, aes(x = overall_score, y = row_label, color = ta_label)) +
    geom_segment(aes(x = 0, xend = overall_score, yend = row_label),
                 linewidth = 0.6, color = "grey80") +
    geom_point(size = 3.5) +
    geom_text(aes(label = sprintf("%.2f", overall_score)),
              hjust = -0.4, size = 3, color = "grey30") +
    scale_color_manual(values = ta_colors, name = "Therapeutic Area") +
    scale_x_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Overall Association Score (Open Targets)", y = NULL) +
    theme_pub(12) +
    theme(axis.text.y = element_text(size = 8.5),
          legend.position = "right",
          legend.text = element_text(size = 9),
          panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3))

  ggsave(file.path(FIG_DIR, "Fig5_PheWAS_gene.pdf"), p5, width = 12, height = p5_height)
  ggsave(file.path(FIG_DIR, "Fig5_PheWAS_gene.png"), p5, width = 12, height = p5_height, dpi = 300)
}

# =============================================================================
# Fig 6: SNP PheWAS — genome-wide significant, top 10 labels
# =============================================================================
cat("Generating Fig 6...\n")

categorize_trait <- function(trait) {
  tl <- tolower(trait)
  if (grepl("immun|autoimmun|inflamm|allerg|asthma|arthrit|lupus|sclerosis", tl)) "Immune/AutoImmune"
  else if (grepl("blood|hemoglobin|platelet|corpuscul|erythro|leuko|iron|ferrit", tl)) "Hematological"
  else if (grepl("lipid|cholesterol|triglycerid|hdl|ldl|fatty", tl)) "Lipid"
  else if (grepl("cancer|carcinoma|neoplasm|tumor|glioma|melanoma", tl)) "Neoplasm"
  else if (grepl("diabet|glucose|insulin|bmi|body mass|obes|metabol|weight|height", tl)) "Metabolic"
  else if (grepl("protein|biomark|complement|cytokine|interleukin", tl)) "Biomarker"
  else if (grepl("stature|anthropom|hip|waist", tl)) "Anthropometric"
  else "Other"
}

phewas_snp[, neglogp_raw := -log10(as.numeric(pval))]
phewas_snp[is.infinite(neglogp_raw), neglogp_raw := 320]
phewas_snp[, trait_category := sapply(trait, categorize_trait)]

phewas_gw <- phewas_snp[as.numeric(pval) < 5e-8]
phewas_gw[, trait_category_f := factor(trait_category)]

pval_cap <- 310
phewas_gw[, neglogp := pmin(neglogp_raw, pval_cap)]

set.seed(123)
phewas_gw[, x_jit := as.numeric(trait_category_f) + runif(.N, -0.25, 0.25)]

top10 <- phewas_gw[order(-neglogp_raw)]
top10 <- top10[!duplicated(paste(gene, trait))][1:min(10, nrow(top10))]

# Trait category colors — harmonized with overall palette
tc_colors <- c(
  "Anthropometric" = "#FC8D62", "Biomarker" = "#FFCC00",
  "Hematological" = "#4DAF4A", "Immune/AutoImmune" = "#2166AC",
  "Lipid" = "#66C2A5", "Metabolic" = "#8DA0CB",
  "Neoplasm" = "#984EA3", "Other" = "#E78AC3"
)

tc_levels <- levels(phewas_gw$trait_category_f)

p6 <- ggplot(phewas_gw, aes(x = x_jit, y = neglogp, color = trait_category_f)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = top10,
                  aes(label = paste0(gene, ": ", trait)),
                  size = 3.2, max.overlaps = 15, segment.size = 0.3,
                  box.padding = 0.6, nudge_y = 5,
                  min.segment.length = 0) +
  scale_x_continuous(breaks = seq_along(tc_levels), labels = tc_levels) +
  scale_color_manual(values = tc_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = expression(-log[10](P))) +
  theme_pub(13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
        legend.position = "none",
        plot.margin = margin(25, 15, 10, 10))

ggsave(file.path(FIG_DIR, "Fig6_PheWAS_snp.pdf"), p6, width = 12, height = 8)
ggsave(file.path(FIG_DIR, "Fig6_PheWAS_snp.png"), p6, width = 12, height = 8, dpi = 300)

cat("All 5 figures generated.\n")
