#!/usr/bin/env Rscript

# ==============================================================================
# MG / OneK1K V4 COMPLETE COMPOSITE, PANEL AND TABLE DELIVERY
# ------------------------------------------------------------------------------
# Purpose:
#   Redraw completed locked-V2 results into a coherent publication package.
#   This script does NOT rerun the 20-GB eQTL analysis, colocalisation, doublet
#   detection, Harmony, edgeR models, fgsea, LODO, or exact permutations.
#
# Output:
#   01_main_figures          four manuscript-oriented figure groups
#   02_supplementary_figures reviewer diagnostics and all-50 Hallmark display
#   03_source_tables         exact plotted source tables
#   04_logs                  status, sessionInfo and source audit
#   MG_V4_publication_figures_tables_DRAFT_HGNC_PENDING.zip
#
# Run in RStudio: open this entire file and click Source.
# ============================================================================== 

options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------ EDIT HERE ------------------------------------

PROJECT_ROOT <- "E:/jd/20260524/py/proteomics/oneK1K"
V2_ROOT <- file.path(PROJECT_ROOT, "yellow_revision_results_v2_locked")
OUTPUT_ROOT <- file.path(V2_ROOT, "18_client_figures_and_tables_v4")

PNG_DPI <- 600
UMAP_MAX_CELLS <- 70000L
VOLCANO_LABEL_N <- 6L
FOREST_N <- 14L  # main figure; the complete ranked table remains in the source archive
SEED <- 20260920L
HGNC_MAPPING_FILE <- file.path(PROJECT_ROOT, "hgnc_gene_symbol_review.csv")
# The 7,852/396 reanalysis baseline is not reconciled to the submitted
# manuscript's 7,598/365. Keep the unreconciled audit in the supplement.
MHC_CHR <- "6"
MHC_START_GRCH37 <- 25726063L
MHC_END_GRCH37 <- 33400644L

# Project threshold requested by the user. FDR is retained and displayed but
# is not substituted for this differential-expression display threshold.
DE_P_THRESHOLD <- 0.05
DE_ABS_LOG2FC_THRESHOLD <- 0.5

# ------------------------------ PACKAGES -------------------------------------

required <- c("data.table", "ggplot2", "patchwork", "scales", "ggrepel", "png")
missing_pkg <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) {
  stop("Install required packages first: ", paste(missing_pkg, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

if (!dir.exists(V2_ROOT)) stop("V2_ROOT does not exist: ", V2_ROOT)

DIR_MAIN <- file.path(OUTPUT_ROOT, "01_main_figures")
DIR_SUPP <- file.path(OUTPUT_ROOT, "02_supplementary_figures")
DIR_TAB  <- file.path(OUTPUT_ROOT, "03_source_tables")
DIR_LOG  <- file.path(OUTPUT_ROOT, "04_logs")
DIR_CODE <- file.path(OUTPUT_ROOT, "05_code")
DIR_PANEL <- file.path(OUTPUT_ROOT, "06_individual_panels")
DIR_TABLE <- file.path(OUTPUT_ROOT, "07_analysis_tables")

# This is a new dedicated output directory. Clearing only its generated
# subfolders prevents stale .pdf.1/.pdf.2 files from surviving a rerun.
for (d in c(DIR_MAIN, DIR_SUPP, DIR_TAB, DIR_LOG, DIR_CODE, DIR_PANEL, DIR_TABLE)) {
  if (dir.exists(d)) unlink(d, recursive = TRUE, force = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
old_archives <- list.files(OUTPUT_ROOT,
  pattern = "^MG_V4_publication_figures_tables.*\\.zip$", full.names = TRUE)
if (length(old_archives)) unlink(old_archives, force = TRUE)

set.seed(SEED)

# ------------------------------ STYLE ----------------------------------------

COL_HC <- "#2C7FB8"
COL_MG <- "#D7301F"
COL_CD4 <- "#1B9E77"
COL_CD8 <- "#6A51A3"
COL_POS <- "#C51B7D"
COL_NEG <- "#2B8CBE"
COL_GOLD <- "#D99A2B"
COL_PURPLE <- "#6A3D9A"
COL_GREY <- "#BDBDBD"
COL_DARK <- "#202020"
COL_GRID <- "#E8E8E8"

GROUP_COLS <- c(HC = COL_HC, MG = COL_MG)
TARGET_COLS <- c(CD4_NC = COL_CD4, CD8_NC = COL_CD8)

theme_pub <- function(base_size = 11) {
  base_size <- max(base_size, 10.8)
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(face = "bold", colour = COL_DARK, size = rel(1.02)),
      axis.text = element_text(colour = COL_DARK, size = rel(0.98)),
      axis.line = element_line(colour = COL_DARK, linewidth = 0.35),
      axis.ticks = element_line(colour = COL_DARK, linewidth = 0.30),
      plot.title = element_text(face = "bold", colour = COL_DARK,
                                size = rel(1.08), hjust = 0),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = COL_DARK),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.96)),
      legend.key.height = grid::unit(0.42, "cm"),
      panel.spacing = grid::unit(0.72, "lines"),
      plot.margin = margin(14, 20, 14, 16)
    )
}

theme_heat <- function(base_size = 9) {
  base_size <- max(base_size, 10.5)
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = COL_DARK, size = rel(0.96)),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(12, 18, 12, 12)
    )
}

save_dual <- function(name, plot, width, height, directory = DIR_MAIN) {
  pdf_file <- file.path(directory, paste0(name, ".pdf"))
  png_file <- file.path(directory, paste0(name, ".png"))
  ggsave(pdf_file, plot, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf, bg = "white")
  ggsave(png_file, plot, width = width, height = height, units = "in",
         dpi = PNG_DPI, limitsize = FALSE, bg = "white")
  invisible(c(pdf_file, png_file))
}

panel_tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 15, colour = COL_DARK),
  plot.tag.position = c(0.01, 0.99)
)

clean_cell <- function(x) {
  x <- as.character(x)
  x <- gsub("_NC$", " NC", x)
  x <- gsub("_", " ", x)
  x
}

clean_pathway <- function(x) {
  x <- sub("^HALLMARK_", "", as.character(x))
  tools::toTitleCase(tolower(gsub("_", " ", x)))
}

# Display-only symbol translation. Never change statistical join keys or
# collapse rows whose old names happen to map to the same approved symbol.
hgnc_mapping <- data.table(old_symbol = character(), approved_symbol = character(),
                           hgnc_id = character())
if (file.exists(HGNC_MAPPING_FILE)) {
  hgnc_mapping <- fread(HGNC_MAPPING_FILE)
  need_columns <- c("old_symbol", "approved_symbol", "hgnc_id")
  if (!all(need_columns %in% names(hgnc_mapping)))
    stop("HGNC mapping requires columns: ", paste(need_columns, collapse = ", "))
  hgnc_mapping[, (need_columns) := lapply(.SD, function(z) trimws(as.character(z))),
               .SDcols = need_columns]
  if (anyNA(hgnc_mapping[, ..need_columns]))
    stop("HGNC mapping has missing fields; fill or remove incomplete rows")
  hgnc_mapping <- hgnc_mapping[nzchar(old_symbol) & nzchar(approved_symbol)]
  if (!nrow(hgnc_mapping))
    warning("HGNC mapping CSV is empty; old symbols remain in draft figures")
  if (anyDuplicated(hgnc_mapping$old_symbol)) stop("Duplicate old_symbol in HGNC mapping")
  if (any(!grepl("^HGNC:[0-9]+$", hgnc_mapping$hgnc_id)))
    stop("Each mapping must include a verified HGNC ID (HGNC:digits)")
  fwrite(hgnc_mapping, file.path(DIR_LOG, "HGNC_display_mapping_used.csv"))
} else {
  warning("HGNC mapping not supplied; old symbols remain in figures: ",
          HGNC_MAPPING_FILE)
}
display_gene <- function(x) {
  z <- as.character(x)
  i <- match(z, hgnc_mapping$old_symbol)
  z[!is.na(i)] <- hgnc_mapping$approved_symbol[i[!is.na(i)]]
  z
}

# -------------------------- INPUT DISCOVERY ----------------------------------

norm_slash <- function(x) gsub("\\\\", "/", normalizePath(x, winslash = "/",
                                                        mustWork = FALSE))

find_input <- function(pattern, required = TRUE, prefer = NULL) {
  x <- list.files(V2_ROOT, pattern = pattern, recursive = TRUE,
                  full.names = TRUE, ignore.case = TRUE)
  # Always read locked analysis stages, never recent visualisation exports.
  rendered_roots <- c(file.path(V2_ROOT, "17_client_requested_figures"), OUTPUT_ROOT)
  for (rendered in rendered_roots)
    x <- x[!startsWith(norm_slash(x), paste0(norm_slash(rendered), "/"))]
  relative_input <- substring(norm_slash(x), nchar(norm_slash(V2_ROOT)) + 2L)
  x <- x[!grepl("^(1[3-8]_|99_)", relative_input)]
  if (length(x)) {
    out_norm <- norm_slash(OUTPUT_ROOT)
    x <- x[!startsWith(norm_slash(x), paste0(out_norm, "/"))]
  }
  if (length(x) && !is.null(prefer)) {
    priority <- ifelse(grepl(prefer, norm_slash(x), ignore.case = TRUE), 0L, 1L)
    # Prefer the requested analysis stage, then the most recently modified
    # copy within that stage. Ascending mtime silently selected stale tables.
    x <- x[order(priority, -as.numeric(file.info(x)$mtime))]
  } else if (length(x)) {
    x <- x[order(file.info(x)$mtime, decreasing = TRUE)]
  }
  if (!length(x)) {
    if (required) stop("Required input not found; pattern: ", pattern)
    return(NA_character_)
  }
  x[1L]
}

read_input <- function(pattern, required = TRUE, prefer = NULL) {
  f <- find_input(pattern, required, prefer)
  if (is.na(f)) return(NULL)
  x <- tryCatch(fread(f), error = function(e) {
    if (required) stop("Could not read ", f, ": ", e$message)
    warning("Could not read optional input ", f, ": ", e$message)
    NULL
  })
  if (!is.null(x)) attr(x, "source_file") <- f
  x
}

need <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, " lacks columns: ", paste(missing, collapse = ", "))
}

copy_plotted_table <- function(x, filename) {
  fwrite(x, file.path(DIR_TAB, filename))
}

source_audit <- data.table(
  object = character(), pattern = character(), source_file = character(),
  status = character()
)

register_source <- function(object, pattern, required = TRUE, prefer = NULL) {
  f <- find_input(pattern, required = FALSE, prefer = prefer)
  source_audit <<- rbind(
    source_audit,
    data.table(object = object, pattern = pattern,
               source_file = ifelse(is.na(f), "", f),
               status = ifelse(is.na(f), ifelse(required, "MISSING_REQUIRED", "MISSING_OPTIONAL"),
                               "FOUND"))
  )
  if (required && is.na(f)) stop("Required input missing for ", object, ": ", pattern)
  f
}

status <- data.table(
  figure = character(), status = character(), details = character()
)

mark <- function(figure, state, details) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", figure,
          " | ", state, " | ", details)
  status <<- rbind(status, data.table(figure = figure, status = state,
                                     details = details))
}
mark("HGNC gene-symbol input", if (nrow(hgnc_mapping)) "INPUT_PROVIDED" else
       "DRAFT_PENDING_MAPPING",
     if (nrow(hgnc_mapping)) paste(nrow(hgnc_mapping), "rows; review HGNC IDs and rendered names")
     else "No client HGNC table supplied; figures with old symbols cannot be final")

# ----------------------- REQUIRED TABLES -------------------------------------

mhc <- read_input("^MHC_exclusion_before_after_counts\\.csv$",
                  prefer = "07_reviewer_completion")
post_counts <- read_input("^Pair_variant_locus_discovery_counts\\.csv$",
                          prefer = "07_reviewer_completion")
rescale <- read_input("^eQTL_rho_to_per_allele_rescaling\\.csv\\.gz$",
                      prefer = "07_reviewer_completion")
joint <- read_input("^MR_coloc_joint_evidence_by_gene_cell\\.csv\\.gz$",
                    prefer = "07_reviewer_completion")
coloc <- read_input("^Coloc_with_conditional_shared_signal_probability\\.csv\\.gz$",
                    prefer = "09_final_evidence_consolidation")
inst <- read_input("^All_postMHC_instrument_strength_harmonisation\\.csv\\.gz$",
                   prefer = "09_final_evidence_consolidation")
expr <- read_input("^Corrected_SESN3_logCPM_by_donor\\.csv$",
                   prefer = "06_annotation_corrected_sensitivity")
lodo <- read_input("^Corrected_SESN3_full_and_LODO\\.csv$",
                   prefer = "06_annotation_corrected_sensitivity")
perm <- read_input("^Genomewide_exact_permutation_statistics\\.csv\\.gz$",
                   prefer = "11_exact_permutation_calibration")
sesn3_exact <- read_input("^SESN3_genomewide_exact_permutation_summary\\.csv$",
                          prefer = "11_exact_permutation_calibration")
de <- read_input("^Corrected_DE_both_target_celltypes\\.csv\\.gz$",
                 prefer = "06_annotation_corrected_sensitivity")
de_conc <- read_input("^CD4NC_CD8NC_DE_concordance_all_genes\\.csv\\.gz$",
                      prefer = "08_extended_sensitivity")
hallmark <- read_input("^(Corrected_Hallmark_50x2_complete|Hallmark_all50_both_celltypes_FDR)\\.csv$",
                       prefer = "06_annotation_corrected_sensitivity")
hallmark_lodo <- read_input("^Hallmark_full_and_LODO_all_results\\.csv\\.gz$",
                            prefer = "10_genomewide_LODO_composition")
lib_qc <- read_input("^Pseudobulk_library_QC\\.csv$",
                     prefer = "10_genomewide_LODO_composition")
comp <- read_input("^Single_cell_composition_plot_data\\.csv$", required = FALSE)
sample_qc <- read_input("^Single_cell_per_sample_QC_by_group\\.csv$", required = FALSE)
pca_comp <- read_input("^Celltype_composition_CLR_PCA\\.csv$", required = FALSE)
pb_pca <- read_input("^Pseudobulk_PCA_coordinates\\.csv$", required = FALSE)
hier <- read_input("^Hierarchical_multiple_testing_summary\\.csv$", required = FALSE)
robust <- read_input("^Genomewide_LODO_robustness_counts\\.csv$", required = FALSE)
exact_all <- read_input("^Genomewide_exact_permutation_gene_summary\\.csv\\.gz$", required = FALSE)
phe <- read_input("^PheWAS_descriptive_categories\\.csv$", required = FALSE)
dice <- read_input("^DICE_historical_claim_audit_summary\\.csv$", required = FALSE)
evidence <- read_input("^Final_cross_modality_gene_cell_evidence\\.csv\\.gz$", required = FALSE)

input_objects <- list(
  MHC = mhc, postMHC_discovery_counts = post_counts,
  eQTL_rescaling = rescale, MR_coloc = joint,
  conditional_coloc = coloc, instruments = inst, SESN3_expression = expr,
  SESN3_LODO = lodo, exact_permutations = perm, SESN3_exact = sesn3_exact,
  pseudobulk_DE = de, DE_concordance = de_conc, Hallmark = hallmark,
  Hallmark_LODO = hallmark_lodo, pseudobulk_QC = lib_qc,
  cell_composition = comp, sample_QC = sample_qc, composition_PCA = pca_comp,
  pseudobulk_PCA = pb_pca, hierarchical_testing = hier,
  genomewide_LODO = robust, exact_gene_summary = exact_all,
  PheWAS = phe, DICE = dice, cross_modality = evidence
)
source_audit <- rbindlist(lapply(names(input_objects), function(nm) {
  x <- input_objects[[nm]]
  data.table(
    object = nm,
    pattern = "",
    source_file = if (is.null(x)) "" else as.character(attr(x, "source_file")),
    status = if (is.null(x)) "MISSING_OPTIONAL" else "FOUND"
  )
}), fill = TRUE)

# Drug annotation is a review queue, never a treatment recommendation. Group
# repeated cell-type rows by approved gene label and preserve the original
# source strings so each indication and mechanism can be checked manually.
drug_source <- read_input("^Final_prioritised_genetic_evidence\\.csv$",
                          required = FALSE, prefer = "09_final_evidence_consolidation")
if (!is.null(drug_source)) {
  need(drug_source, c("gene", "cell_type", "rsid", "MR_FDR_global",
                      "PP_H4", "historical_n_drugs", "approved_drugs"),
       "Final genetic drug review source")
  drug_source[, approved_display_gene := display_gene(gene)]
  drug_source[, historical_n_drugs := suppressWarnings(as.numeric(historical_n_drugs))]
  drug_queue <- drug_source[is.finite(MR_FDR_global) & MR_FDR_global < 0.05 &
    is.finite(PP_H4) & PP_H4 > 0.8]
  drug_queue[, `:=`(source_is_postMHC = TRUE,
    indication_for_MG_verified = FALSE,
    mechanism_direction_verified = FALSE,
    drug_claim_status = "Manual pharmacology and indication verification required")]
  fwrite(drug_queue, file.path(DIR_TAB,
    "PostMHC_strong_gene_drug_manual_review_queue.csv"))
  historical_flagged <- drug_queue[historical_n_drugs > 0]
  fwrite(historical_flagged, file.path(DIR_TAB,
    "Historical_drug_mentions_in_strong_postMHC_pairs_NOT_MG_validation.csv"))
  mark("Post-MHC drug annotation", "REVIEW_REQUIRED",
       paste(uniqueN(drug_queue$gene), "strong-evidence genes;",
             uniqueN(historical_flagged$gene), "historically drug-annotated genes"))
} else {
  mark("Post-MHC drug annotation", "BLOCKED", "Final prioritised evidence table unavailable")
}

# ---------------- MAIN FIGURE 1: GENETIC REANALYSIS --------------------------

need(post_counts, c("discovery_set", "gene_cell_pairs", "distinct_lead_variants",
                    "position_based_1Mb_loci"), "Post-MHC discovery counts")
primary_counts <- post_counts[discovery_set == "MR_global_FDR005"]
if (nrow(primary_counts) != 1L) stop("Expected one primary post-MHC discovery row")
post_display <- data.table(
  label = factor(c("Gene-cell pairs", "Distinct lead variants", "Physical 1-Mb groups"),
                 levels = rev(c("Gene-cell pairs", "Distinct lead variants", "Physical 1-Mb groups"))),
  n = as.numeric(unlist(primary_counts[, .(gene_cell_pairs,
                distinct_lead_variants, position_based_1Mb_loci)]))
)
p_post <- ggplot(post_display, aes(n, label)) +
  geom_col(fill = "#3E78A8", width = 0.60, colour = COL_DARK, linewidth = 0.32) +
  geom_text(aes(label = n), hjust = -0.22, size = 4.1, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Post-MHC genetic findings", subtitle = "Global BH FDR < 0.05",
       x = "Count", y = NULL,
       caption = "Physical 1-Mb groups are not LD-independent loci.") +
  theme_pub(11.4) + theme(axis.text.y = element_text(size = 10.6))

need(mhc, c("stage", "metric", "value"), "MHC count table")
mhc[, stage := factor(stage, levels = c("Before MHC exclusion", "After MHC exclusion"))]
mhc[, metric_label := fifelse(
  grepl("tested", metric, ignore.case = TRUE), "Tested gene-cell pairs",
  fifelse(grepl("FDR", metric, ignore.case = TRUE), "Global FDR < 0.05", metric)
)]
mhc_display <- mhc[metric %in% c("tested_gene_cell_pairs", "FDR_significant_pairs")]
if (nrow(mhc_display) != 4L) stop("Before/after MHC audit needs two comparable metrics in both stages")
post_sig_from_audit <- as.numeric(mhc_display[
  stage == "After MHC exclusion" & metric == "FDR_significant_pairs", value])
if (length(post_sig_from_audit) != 1L || !is.finite(post_sig_from_audit) ||
    post_sig_from_audit != primary_counts$gene_cell_pairs)
  stop("Source conflict: post-MHC significant gene-cell count disagrees between MHC audit and discovery table")
mhc_display[, metric_label := factor(metric_label,
  levels = c("Tested gene-cell pairs", "Global FDR < 0.05"))]
p_mhc <- ggplot(mhc_display, aes(stage, value, fill = stage)) +
  geom_col(width = 0.62, colour = COL_DARK, linewidth = 0.32) +
  geom_text(aes(label = scales::comma(value)), vjust = -0.35, size = 3.2,
            fontface = "bold") +
  facet_wrap(~metric_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Before MHC exclusion" = "#9E9E9E",
                               "After MHC exclusion" = "#2C7FB8")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "MHC exclusion in the current analysis pipeline",
       subtitle = "Current pipeline baseline differs from the submitted manuscript; reconcile original filters before citation.",
       x = NULL, y = "Count") +
  theme_pub(12) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10.5),
        strip.text = element_text(size = 11.5))

need(rescale, c("rho", "beta_eqtl", "rescaling_factor"), "eQTL rescaling table")
rescale[, `:=`(rho = as.numeric(rho), beta_eqtl = as.numeric(beta_eqtl),
               rescaling_factor = as.numeric(rescaling_factor))]
rescale_plot <- rescale[is.finite(rho) & is.finite(beta_eqtl) &
                          is.finite(rescaling_factor)]
if (nrow(rescale_plot) > 30000L) rescale_plot <- rescale_plot[sample(.N, 30000L)]
p_rescale <- ggplot(rescale_plot, aes(rho, beta_eqtl, colour = rescaling_factor)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 0.65, alpha = 0.38) +
  scale_colour_viridis_c(option = "C", begin = 0.10, end = 0.92) +
  labs(title = "OneK1K effect-size rescaling",
       subtitle = "Published Spearman rho converted to approximate per-allele effect",
       x = "Published Spearman rho", y = "Corrected eQTL effect",
       colour = "Scale factor") +
  theme_pub(10.5)

need(joint, c("gene", "cell_type", "FDR_global_recomputed", "PP.H4"),
     "MR-coloc evidence")
joint[, `:=`(
  FDR_global_recomputed = as.numeric(FDR_global_recomputed),
  PP.H4 = as.numeric(PP.H4)
)]
joint_plot <- joint[is.finite(FDR_global_recomputed) & is.finite(PP.H4)]
joint_plot[, `:=`(
  neglog10_FDR = pmin(-log10(pmax(FDR_global_recomputed, 1e-300)), 10),
  evidence = fcase(
    FDR_global_recomputed < 0.05 & PP.H4 > 0.80, "Strong convergent",
    FDR_global_recomputed < 0.05 & PP.H4 > 0.50, "Moderate convergent",
    FDR_global_recomputed < 0.05, "MR only",
    default = "Other"
  )
)]
joint_plot[, evidence := factor(evidence,
                                levels = c("Other", "MR only", "Moderate convergent",
                                           "Strong convergent"))]
joint_lab <- joint_plot[evidence == "Strong convergent"][order(FDR_global_recomputed,
                                                                 -PP.H4)]
joint_lab <- joint_lab[!duplicated(display_gene(gene))][1:min(.N, 5L)]
joint_lab[, display_label := display_gene(gene)]
p_joint <- ggplot(joint_plot, aes(neglog10_FDR, PP.H4, colour = evidence)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "#666666") +
  geom_hline(yintercept = c(0.50, 0.80), linetype = c(3, 2),
             colour = c("#999999", COL_PURPLE)) +
  geom_point(size = 2.2, alpha = 0.88) +
  ggrepel::geom_text_repel(data = joint_lab, aes(label = display_label),
                           colour = COL_DARK, size = 3.8,
                           box.padding = 0.7, point.padding = 0.3,
                           force = 5, force_pull = 0.45,
                           max.iter = 50000, max.overlaps = Inf,
                           min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = c("Other" = COL_GREY, "MR only" = COL_GOLD,
                                 "Moderate convergent" = "#E08214",
                                 "Strong convergent" = COL_PURPLE)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(title = "MR-colocalisation convergence",
       subtitle = "Only variants with evaluable colocalisation are shown",
       x = expression(-log[10](global~BH~FDR)), y = "PP.H4", colour = NULL) +
  theme_pub(11.3) + theme(legend.position = "bottom")

need(coloc, c("gene", "cell_type", "PP.H4", "or", "or_lci", "or_uci", "pval"),
     "Conditional coloc table")
forest <- coloc[is.finite(PP.H4) & PP.H4 > 0.80 & is.finite(or) &
                  is.finite(or_lci) & is.finite(or_uci)]
setorder(forest, pval, -PP.H4)
forest <- forest[1:min(.N, FOREST_N)]
forest[, row_label := paste0(display_gene(gene), " (", clean_cell(cell_type), ")")]
# If two old symbols share an approved display name, retain the old symbol
# as a disambiguator and do not merge their estimates.
forest[, plot_label := fifelse(duplicated(row_label) | duplicated(row_label, fromLast = TRUE),
                              paste0(row_label, " [", gene, "]"), row_label)]
if (anyDuplicated(forest$plot_label)) stop("Ambiguous duplicate forest labels")
forest[, plot_label := factor(plot_label, levels = rev(plot_label))]
p_forest <- ggplot(forest, aes(or, plot_label, colour = cell_type)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#666666") +
  geom_errorbarh(aes(xmin = or_lci, xmax = or_uci), height = 0.16,
                 linewidth = 0.50) +
  geom_point(aes(size = PP.H4), shape = 16) +
  scale_x_log10() +
  scale_colour_manual(values = TARGET_COLS, labels = clean_cell) +
  scale_size_continuous(range = c(2.0, 3.4), limits = c(0.8, 1)) +
  labs(title = "Strong colocalised effects",
       subtitle = "Top results ranked by corrected MR P value",
       x = "Odds ratio per corrected exposure unit (95% CI)", y = NULL,
       colour = NULL, size = "PP.H4") +
  theme_pub(11.3) + theme(legend.position = "bottom",
                          axis.text.y = element_text(size = 10.6))

fig1 <- (p_post | p_rescale) / (p_joint | p_forest) +
  plot_layout(heights = c(0.78, 1.22)) +
  plot_annotation(title = "Corrected post-MHC genetic analysis",
                  tag_levels = "A") & panel_tag_theme
save_dual("Figure1_genetic_reanalysis", fig1, 18.0, 12.2)
save_dual("Supplementary_MHC_reanalysis_baseline_UNRECONCILED", p_mhc,
          12.5, 5.2, DIR_SUPP)
copy_plotted_table(mhc, "Supplementary_MHC_counts_source.csv")
copy_plotted_table(post_display, "Figure1A_postMHC_discovery_counts.csv")
copy_plotted_table(rescale, "Figure1B_eQTL_rescaling.csv.gz")
copy_plotted_table(joint_plot, "Figure1C_MR_coloc_plot_data.csv")
copy_plotted_table(forest, "Figure1D_strong_coloc_forest.csv")
mark("Figure 1", "DONE", "Effect correction, MR-coloc and forest; unreconciled MHC baseline in supplement")

# Exploratory (post hoc) boundary sensitivity of the POST-MHC MR test family.
# Report BOTH original-BH retention and BH recalculated over each smaller
# family. These are not reruns of eQTL, colocalisation or GWAS analyses.
mr_boundary <- read_input("^MR_postMHC_global_and_within_cell_BH\\.csv\\.gz$",
                          required = FALSE, prefer = "07_reviewer_completion")
if (!is.null(mr_boundary)) {
  need(mr_boundary, c("chr", "pos", "gene", "cell_type", "rsid",
                      "FDR_global_recomputed"), "MHC boundary source")
  mr_boundary[, chr_clean := gsub("^chr", "", as.character(chr), ignore.case = TRUE)]
  mr_boundary[, pos_clean := suppressWarnings(as.numeric(pos))]
  mr_boundary[, q_original := suppressWarnings(as.numeric(FDR_global_recomputed))]
  if (nrow(mr_boundary) != 7328L) {
    warning("MR test family has ", nrow(mr_boundary),
            " rows, not the reviewed 7,328; inspect source before comparing counts")
  }
  if (any(is.na(mr_boundary$chr_clean) | !is.finite(mr_boundary$pos_clean)))
    stop("MR boundary table has unknown chromosome/position; cannot calculate sensitivity")
  boundary_rows <- list()
  for (flank in c(0L, 250000L, 500000L, 1000000L)) {
    keep <- !(mr_boundary$chr_clean == MHC_CHR &
      mr_boundary$pos_clean >= MHC_START_GRCH37 - flank &
      mr_boundary$pos_clean <= MHC_END_GRCH37 + flank)
    x <- mr_boundary[keep]
    if (!"pval" %in% names(x)) stop("MR boundary audit requires raw pval to recalculate BH")
    q_retested <- p.adjust(suppressWarnings(as.numeric(x$pval)), method = "BH")
    boundary_rows[[as.character(flank)]] <- data.table(
      added_flank_kb = flank / 1000,
      tested_postMHC_pairs = nrow(x),
      pairs_with_original_global_BH_q_below_005 = sum(x$q_original < 0.05, na.rm = TRUE),
      pairs_with_recalculated_global_BH_q_below_005 = sum(q_retested < 0.05, na.rm = TRUE),
      distinct_significant_lead_rsids = uniqueN(x[q_original < 0.05 &
        !is.na(rsid) & nzchar(as.character(rsid))]$rsid)
    )
  }
  boundary_summary <- rbindlist(boundary_rows)
  fwrite(boundary_summary, file.path(DIR_TAB, "MHC_boundary_exploratory_sensitivity.csv"))
  near_boundary <- mr_boundary[chr_clean == MHC_CHR & q_original < 0.05 &
    pos_clean >= MHC_START_GRCH37 - 500000 & pos_clean < MHC_START_GRCH37]
  near_boundary[, distance_below_MHC_start_bp := MHC_START_GRCH37 - pos_clean]
  setorder(near_boundary, distance_below_MHC_start_bp)
  fwrite(near_boundary[, .(gene, cell_type, rsid, chr = chr_clean, pos = pos_clean,
     distance_below_MHC_start_bp, q_original)],
    file.path(DIR_TAB, "MHC_boundary_significant_pairs_within_500kb.csv"))
  mark("MHC boundary sensitivity", "EXPLORATORY",
       paste("Post hoc 0/250/500/1000 kb extensions;", nrow(near_boundary),
             "significant pairs within 500 kb of lower boundary"))
} else {
  mark("MHC boundary sensitivity", "BLOCKED", "Post-MHC MR per-pair table missing")
}

# ---------------- MAIN FIGURE 2: SESN3 ---------------------------------------

need(expr, c("donor_id", "group", "cell_type", "SESN3_logCPM"), "SESN3 expression")
expr[, group := factor(group, levels = c("HC", "MG"))]
expr[, cell_label := clean_cell(cell_type)]
p_expr <- ggplot(expr, aes(group, SESN3_logCPM, colour = group, fill = group)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, alpha = 0.14,
               colour = COL_DARK, linewidth = 0.48) +
  geom_point(size = 2.8, position = position_jitter(width = 0.045)) +
  ggrepel::geom_text_repel(aes(label = donor_id), size = 3.0,
                           direction = "y", max.overlaps = Inf,
                           box.padding = 0.20, point.padding = 0.15,
                           min.segment.length = 0, show.legend = FALSE) +
  facet_wrap(~cell_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = GROUP_COLS) +
  scale_fill_manual(values = GROUP_COLS) +
  labs(title = "SESN3 donor-level expression",
       subtitle = "TMM-normalised pseudobulk; MAIT-like cluster excluded",
       x = NULL, y = expression(log[2]~CPM)) +
  theme_pub(9.5) + theme(legend.position = "none")

need(lodo, c("cell_type", "scenario", "logFC"), "SESN3 LODO table")
lodo[, scenario_label := fifelse(scenario == "Full", "Full model",
                                  paste("Omit", sub("^Leave_out_", "", scenario)))]
scenario_order <- c("Full model", paste("Omit", c("H1", "H2", "H3", "M0707_01",
                                                    "M0709_02", "M0709_03", "M210715_01")))
lodo[, scenario_label := factor(scenario_label, levels = rev(scenario_order))]
lodo[, model_group := fcase(
  scenario == "Full", "Full model",
  omitted_group == "HC", "Omitted HC",
  omitted_group == "MG", "Omitted MG",
  default = "LODO"
)]
lodo[, cell_label := clean_cell(cell_type)]
p_lodo <- ggplot(lodo, aes(logFC, scenario_label, colour = model_group)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#666666") +
  geom_segment(aes(x = 0, xend = logFC, yend = scenario_label),
               linewidth = 0.38, colour = "#D7D7D7") +
  geom_point(size = 2.7) +
  facet_wrap(~cell_label, nrow = 1) +
  scale_colour_manual(values = c("Full model" = COL_DARK, "Omitted HC" = COL_HC,
                                 "Omitted MG" = COL_MG, "LODO" = "#777777")) +
  labs(title = "Leave-one-donor-out stability",
       x = "log2 fold change (MG versus HC)", y = NULL, colour = NULL) +
  theme_pub(9.0) + theme(legend.position = "bottom")

need(perm, c("gene", "cell_type", "logFC", "F", "is_observed"), "Permutation table")
need(sesn3_exact, c("cell_type", "exact_P_F", "exact_P_abs_logFC",
                    "exact_BH_within_celltype"),
     "SESN3 exact permutation summary")
sesn3_perm <- perm[gene == "SESN3"]
if (!setequal(unique(sesn3_perm$cell_type), c("CD4_NC", "CD8_NC")) ||
    any(!is.finite(sesn3_perm$F)) ||
    any(sesn3_perm[, .N, by = cell_type]$N != 35L) ||
    nrow(sesn3_perm[is_observed == TRUE, .N, by = cell_type]) != 2L ||
    any(sesn3_perm[is_observed == TRUE, .N, by = cell_type]$N != 1L))
  stop("SESN3 permutation panel requires 35 allocations and one observed per cell type")
for (ct in unique(sesn3_perm$cell_type)) {
  d <- sesn3_perm[cell_type == ct]
  f_observed <- d[is_observed == TRUE]$F
  logfc_observed <- abs(d[is_observed == TRUE]$logFC)
  p_f <- sum(d$F >= f_observed - 1e-9) / nrow(d)
  p_abs <- sum(abs(d$logFC) >= logfc_observed - 1e-9) / nrow(d)
  summary_row <- sesn3_exact[cell_type == ct]
  if (nrow(summary_row) != 1L ||
      abs(p_f - summary_row$exact_P_F) > 1e-6 ||
      abs(p_abs - summary_row$exact_P_abs_logFC) > 1e-6)
    stop("Fig2C permutation table and exact P summary disagree for ", ct)
}
sesn3_perm[, cell_label := clean_cell(cell_type)]
exact_anno <- copy(sesn3_exact)
exact_anno[, cell_label := clean_cell(cell_type)]
exact_anno[, label := sprintf("Exact P (F) = %.4f\nExact BH q = %.4f",
                              exact_P_F, exact_BH_within_celltype)]
p_perm <- ggplot(sesn3_perm, aes(F)) +
  geom_histogram(bins = 14,
                 fill = "#9ECAE1", colour = "white", linewidth = 0.25) +
  geom_vline(data = sesn3_perm[is_observed == TRUE], aes(xintercept = F),
             colour = COL_MG, linewidth = 0.9) +
  geom_text(data = exact_anno, aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.03, vjust = 1.15, size = 3.35,
            lineheight = 1.15) +
  facet_wrap(~cell_label, nrow = 1, scales = "free_y") +
  labs(title = "Exact donor-label permutation",
       subtitle = "All 35 allocations; red = observed F; exact tail probability uses F",
       x = "edgeR F statistic", y = "Allocations",
       caption = "Absolute-log2FC exact P is a separate statistic (0.0571 in both cell types).") +
  theme_pub(11.2)

# Locate one canonical regional figure per target population. Existing corrected
# regional plots are retained because the V2 archive lacks a complete LD panel.
copy_regional <- function(cell_type) {
  canonical <- paste0("SESN3_regional_coloc_", cell_type)
  pdfs <- list.files(V2_ROOT, pattern = paste0("(^|/)", canonical, "\\.pdf$"),
                     recursive = TRUE, full.names = TRUE)
  pdfs <- pdfs[!startsWith(norm_slash(pdfs), paste0(norm_slash(OUTPUT_ROOT), "/"))]
  pdfs <- pdfs[!startsWith(norm_slash(pdfs),
    paste0(norm_slash(file.path(V2_ROOT, "17_client_requested_figures")), "/"))]
  if (!length(pdfs)) {
    pdfs <- list.files(V2_ROOT,
      pattern = paste0("Y7_", canonical, "\\.pdf$"),
      recursive = TRUE, full.names = TRUE)
  }
  if (!length(pdfs)) return(NULL)
  src_pdf <- pdfs[order(nchar(norm_slash(pdfs)))][1L]
  out_pdf <- file.path(DIR_MAIN, paste0("Figure2_regional_", cell_type, ".pdf"))
  file.copy(src_pdf, out_pdf, overwrite = TRUE)

  pngs <- list.files(V2_ROOT, pattern = paste0("(^|/)", canonical, "\\.png$"),
                     recursive = TRUE, full.names = TRUE)
  pngs <- pngs[!startsWith(norm_slash(pngs), paste0(norm_slash(OUTPUT_ROOT), "/"))]
  out_png <- file.path(DIR_MAIN, paste0("Figure2_regional_", cell_type, ".png"))
  pngs <- pngs[!startsWith(norm_slash(pngs),
    paste0(norm_slash(file.path(V2_ROOT, "17_client_requested_figures")), "/"))]
  if (length(pngs)) file.copy(pngs[order(nchar(norm_slash(pngs)))][1L], out_png,
                              overwrite = TRUE)
  list(pdf = out_pdf, png = if (file.exists(out_png)) out_png else NA_character_)
}

reg4 <- copy_regional("CD4_NC")
reg8 <- copy_regional("CD8_NC")

fig2 <- p_expr / (p_lodo | p_perm) +
  plot_layout(heights = c(0.82, 1.18)) +
  plot_annotation(title = "Donor-level validation and robustness of SESN3",
                  tag_levels = "A") & panel_tag_theme
save_dual("Figure2_SESN3_donor_robustness", fig2, 16.2, 11.4)
copy_plotted_table(expr, "Figure2A_SESN3_expression.csv")
copy_plotted_table(lodo, "Figure2B_SESN3_LODO.csv")
copy_plotted_table(sesn3_perm, "Figure2C_SESN3_exact_permutation.csv.gz")
copy_plotted_table(exact_anno[, .(cell_type, exact_P_F, exact_P_abs_logFC,
                                  exact_BH_within_celltype)],
                   "Figure2C_exact_P_statistic_audit.csv")
mark("Figure 2", "DONE",
     paste("donor expression, LODO, exact permutation; regional PDFs:",
           sum(vapply(list(reg4, reg8), Negate(is.null), logical(1))), "/2"))

# ---------------- MAIN FIGURE 3: SINGLE-CELL VALIDATION ----------------------

p_umap_cell <- p_umap_group <- p_marker <- p_pb_pca <- NULL
seu_file <- find_input("MG_7donor_QC_Harmony_previous440_reference_v2\\.rds$",
                       required = FALSE)
reviewed_meta <- read_input("^Reviewed_cell_metadata\\.csv\\.gz$", required = FALSE)

if (!is.na(seu_file) && !is.null(reviewed_meta) &&
    requireNamespace("Seurat", quietly = TRUE) &&
    requireNamespace("SeuratObject", quietly = TRUE)) {
  message("Loading Seurat checkpoint for V3 UMAP and marker plots: ", seu_file)
  seu <- readRDS(seu_file)
  need(reviewed_meta, c("cell", "donor_id", "group", "reviewed_cell_type"),
       "Reviewed cell metadata")
  idx <- match(colnames(seu), reviewed_meta$cell)
  if (sum(!is.na(idx)) < 0.90 * ncol(seu)) {
    stop("Fewer than 90% of Seurat cells match Reviewed_cell_metadata.csv.gz")
  }
  seu$v3_celltype <- reviewed_meta$reviewed_cell_type[idx]
  seu$v3_group <- reviewed_meta$group[idx]
  seu$v3_donor <- reviewed_meta$donor_id[idx]

  reductions <- names(seu@reductions)
  reduction_use <- if ("umap" %in% reductions) "umap" else reductions[grepl("umap", reductions)][1L]
  if (is.na(reduction_use) || !nzchar(reduction_use)) stop("UMAP reduction not found in Seurat object")
  emb <- as.data.table(SeuratObject::Embeddings(seu, reduction = reduction_use),
                       keep.rownames = "cell")
  setnames(emb, names(emb)[2:3], c("UMAP_1", "UMAP_2"))
  meta_plot <- reviewed_meta[match(emb$cell, cell)]
  emb[, `:=`(cell_type = meta_plot$reviewed_cell_type,
             donor_id = meta_plot$donor_id, group = meta_plot$group)]
  if (nrow(emb) > UMAP_MAX_CELLS) emb_plot <- emb[sample(.N, UMAP_MAX_CELLS)] else emb_plot <- emb

  cell_levels <- sort(unique(na.omit(emb_plot$cell_type)))
  cell_cols <- setNames(grDevices::hcl.colors(length(cell_levels), "Dark 3"), cell_levels)
  centers <- emb_plot[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)),
                      by = cell_type]
  p_umap_cell <- ggplot(emb_plot, aes(UMAP_1, UMAP_2, colour = cell_type)) +
    geom_point(size = 0.12, alpha = 0.70) +
    ggrepel::geom_text_repel(data = centers, aes(label = clean_cell(cell_type)),
                             size = 3.1, fontface = "bold", colour = COL_DARK,
                             max.overlaps = Inf, min.segment.length = 0,
                             box.padding = 0.28, show.legend = FALSE) +
    scale_colour_manual(values = cell_cols) +
    labs(title = "Reviewed cell-type annotation", x = "UMAP 1", y = "UMAP 2") +
    theme_pub(9.0) + theme(legend.position = "none")

  p_umap_group <- ggplot(emb_plot, aes(UMAP_1, UMAP_2, colour = group)) +
    geom_point(size = 0.14, alpha = 0.46) +
    scale_colour_manual(values = GROUP_COLS, na.value = COL_GREY) +
    labs(title = "Disease-group distribution", x = "UMAP 1", y = "UMAP 2",
         colour = NULL) +
    theme_pub(9.0) + theme(legend.position = "bottom")

  marker_panel <- unique(c(
    "IL7R", "CCR7", "LTB", "CD8A", "CD8B", "CCL5", "NKG7", "GNLY",
    "TRAC", "TRDC", "KLRB1", "SLC4A10", "MS4A1", "CD79A", "CD37",
    "MZB1", "JCHAIN", "LYZ", "LST1", "FCGR3A", "PPBP", "PF4",
    "HDC", "MS4A2", "MKI67", "TOP2A"
  ))
  marker_panel <- marker_panel[marker_panel %in% rownames(seu)]
  if (length(marker_panel) >= 8L) {
    assay_use <- if ("RNA" %in% SeuratObject::Assays(seu)) "RNA" else
      SeuratObject::DefaultAssay(seu)
    SeuratObject::DefaultAssay(seu) <- assay_use
    p_marker <- Seurat::DotPlot(
      seu, features = marker_panel, group.by = "v3_celltype",
      assay = assay_use, dot.scale = 5.5
    ) +
      scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7",
                             high = "#B2182B", midpoint = 0) +
      labs(title = "Canonical marker validation", x = NULL, y = NULL,
           colour = "Scaled\nexpression", size = "% cells") +
      theme_pub(10.2) +
      theme(axis.text.x = element_text(angle = 48, hjust = 1,
                                       face = "italic", size = 9.4),
            axis.text.y = element_text(size = 9.5))
  }
  fwrite(emb, file.path(DIR_TAB, "Figure3_UMAP_coordinates.csv.gz"))
  rm(seu); invisible(gc())
}

if (!is.null(pb_pca)) {
  if (all(c("PC1", "PC2", "group", "cell_type", "donor_id") %in% names(pb_pca))) {
    pb_pca[, cell_label := clean_cell(cell_type)]
    p_pb_pca <- ggplot(pb_pca, aes(PC1, PC2, colour = group, label = donor_id)) +
      geom_point(size = 3.0) +
      ggrepel::geom_text_repel(size = 2.7, max.overlaps = Inf,
                               min.segment.length = 0) +
      facet_wrap(~cell_label, scales = "free", nrow = 1) +
      scale_colour_manual(values = GROUP_COLS) +
      labs(title = "Donor-level pseudobulk PCA", x = "PC1", y = "PC2",
           colour = NULL) +
      theme_pub(9.0) + theme(legend.position = "bottom")
  }
}

if (is.null(p_pb_pca) && !is.null(pca_comp)) {
  need(pca_comp, c("donor_id", "group", "PC1", "PC2"), "Composition PCA")
  p_pb_pca <- ggplot(pca_comp, aes(PC1, PC2, colour = group, label = donor_id)) +
    geom_point(size = 3.1) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = Inf,
                             min.segment.length = 0) +
    scale_colour_manual(values = GROUP_COLS) +
    labs(title = "Cell-composition CLR PCA", x = "PC1", y = "PC2",
         colour = NULL) +
    theme_pub(9.0) + theme(legend.position = "bottom")
}

sc_panels <- Filter(Negate(is.null), list(p_umap_cell, p_umap_group, p_marker, p_pb_pca))
if (length(sc_panels) >= 2L) {
  fig3 <- wrap_plots(sc_panels, ncol = 2, guides = "collect") +
    plot_annotation(title = "Single-cell annotation and donor-level quality assessment",
                    tag_levels = "A") & panel_tag_theme
  save_dual("Figure3_single_cell_validation", fig3, 18.4,
            ifelse(length(sc_panels) >= 4L, 13.0, 7.6))
  mark("Figure 3", "DONE", paste(length(sc_panels), "panels"))
} else {
  mark("Figure 3", "PARTIAL", "Seurat RDS/metadata or PCA coordinates unavailable")
}

# ---------------- MAIN FIGURE 4: TRANSCRIPTOMIC VALIDATION -------------------

need(de, c("gene", "logFC", "PValue", "cell_type"), "Corrected DE table")
de_full <- de[scenario == "Full" | is.na(scenario)]
de_full[, class := fcase(
  PValue < DE_P_THRESHOLD & logFC >= DE_ABS_LOG2FC_THRESHOLD, "Higher in MG",
  PValue < DE_P_THRESHOLD & logFC <= -DE_ABS_LOG2FC_THRESHOLD, "Lower in MG",
  default = "Not primary"
)]
de_full[, neglog10P := -log10(pmax(PValue, 1e-300))]
de_full[, cell_label := clean_cell(cell_type)]
lab_de <- de_full[class != "Not primary", .SD[order(PValue)][1:min(.N, VOLCANO_LABEL_N)],
                  by = cell_type]
lab_de[, display_label := display_gene(gene)]
p_volcano <- ggplot(de_full, aes(logFC, neglog10P, colour = class)) +
  geom_vline(xintercept = c(-DE_ABS_LOG2FC_THRESHOLD, DE_ABS_LOG2FC_THRESHOLD),
             linetype = 2, colour = "#777777") +
  geom_hline(yintercept = -log10(DE_P_THRESHOLD), linetype = 2,
             colour = "#777777") +
  geom_point(size = 0.55, alpha = 0.52) +
  ggrepel::geom_text_repel(data = lab_de, aes(label = display_label), size = 3.3,
                           colour = COL_DARK, max.overlaps = Inf,
                           box.padding = 0.45, point.padding = 0.25,
                           force = 2.5, min.segment.length = 0,
                           show.legend = FALSE) +
  facet_wrap(~cell_label, nrow = 1, scales = "free") +
  scale_colour_manual(values = c("Lower in MG" = COL_NEG,
                                 "Not primary" = "#C8C8C8",
                                 "Higher in MG" = COL_POS)) +
  labs(title = "Donor-level pseudobulk differential expression",
       subtitle = "Primary display threshold: raw P < 0.05 and |log2FC| >= 0.5",
       x = "log2 fold change (MG versus HC)", y = expression(-log[10](P)),
       colour = NULL) +
  theme_pub(10.2) + theme(legend.position = "bottom")

need(de_conc, c("gene", "logFC_CD4_NC", "logFC_CD8_NC", "evidence_class"),
     "DE concordance table")
conc <- de_conc[is.finite(logFC_CD4_NC) & is.finite(logFC_CD8_NC)]
rho_de <- cor(conc$logFC_CD4_NC, conc$logFC_CD8_NC, method = "spearman")
conc[, highlight := fifelse(gene == "SESN3", "SESN3",
                            fifelse(grepl("shared", evidence_class, ignore.case = TRUE),
                                    "Shared primary", "Other"))]
conc_lab <- conc[highlight != "Other"][order(-abs(logFC_CD4_NC) - abs(logFC_CD8_NC))]
conc_lab <- unique(rbind(conc[gene == "SESN3"], conc_lab[1:min(.N, 10L)]), by = "gene")
conc_lab[, display_label := display_gene(gene)]
p_conc <- ggplot(conc, aes(logFC_CD4_NC, logFC_CD8_NC, colour = highlight)) +
  geom_hline(yintercept = 0, colour = "#DDDDDD") +
  geom_vline(xintercept = 0, colour = "#DDDDDD") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 0.70, alpha = 0.56) +
  ggrepel::geom_text_repel(data = conc_lab, aes(label = display_label), colour = COL_DARK,
                           size = 3.1, max.overlaps = Inf, min.segment.length = 0,
                           box.padding = 0.45, point.padding = 0.22,
                           show.legend = FALSE) +
  scale_colour_manual(values = c("Other" = "#C7C7C7", "SESN3" = COL_MG,
                                 "Shared primary" = COL_PURPLE)) +
  labs(title = "Cross-cell-type effect concordance",
       subtitle = sprintf("All-gene Spearman rho = %.3f", rho_de),
       x = "CD4 NC log2FC", y = "CD8 NC log2FC", colour = NULL) +
  coord_equal() + theme_pub(10.2) + theme(legend.position = "bottom")

need(hallmark, c("cell_type", "pathway", "NES", "FDR_within_celltype"),
     "Hallmark table")
hm <- hallmark[test_status == "Evaluated" | is.na(test_status)]
hm[, pathway_label := clean_pathway(pathway)]
hm[, significance := FDR_within_celltype < 0.05]
hm_score <- hm[, .(rank_score = max(abs(NES), na.rm = TRUE),
                   any_fdr = any(significance, na.rm = TRUE)), by = pathway_label]
setorder(hm_score, -any_fdr, -rank_score)
main_pathways <- hm_score[1:min(.N, 20L)]$pathway_label
hm_main <- hm[pathway_label %in% main_pathways]
hm_main[, pathway_label := factor(pathway_label, levels = rev(main_pathways))]
hm_main[, cell_label := clean_cell(cell_type)]
p_hm_main <- ggplot(hm_main, aes(cell_label, pathway_label)) +
  geom_point(aes(size = -log10(pmax(FDR_within_celltype, 1e-6)), fill = NES,
                 shape = significance), colour = COL_DARK, stroke = 0.32) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 22)) +
  scale_size_continuous(range = c(1.6, 5.7)) +
  labs(title = "Hallmark pathway activity",
       subtitle = "Top 20; squares indicate FDR < 0.05",
       x = NULL, y = NULL, fill = "NES", size = expression(-log[10](FDR)),
       shape = "FDR < 0.05") +
  guides(shape = "none", size = "none") +
  theme_heat(10.0) + theme(axis.text.y = element_text(size = 9.0),
                           legend.position = "bottom")

need(hallmark_lodo, c("pathway", "cell_type", "scenario", "NES"),
     "Hallmark LODO table")
hm_lodo <- hallmark_lodo[scenario != "Full" & is.finite(NES)]
hm_full <- hallmark_lodo[scenario == "Full" & is.finite(NES),
                         .(pathway, cell_type, full_NES = NES,
                           full_FDR = fifelse(is.finite(FDR_model), FDR_model, padj))]
hm_stab <- merge(
  hm_lodo[, .(median_LODO_NES = median(NES),
              min_LODO_NES = min(NES), max_LODO_NES = max(NES),
              direction_concordance = NA_real_), by = .(pathway, cell_type)],
  hm_full, by = c("pathway", "cell_type"), all.x = TRUE
)
hm_dir <- merge(hm_lodo[, .(pathway, cell_type, NES)], hm_full,
                by = c("pathway", "cell_type"), all.x = TRUE)
hm_dir[, concordant := sign(NES) == sign(full_NES)]
hm_dir_sum <- hm_dir[, .(direction_concordance = mean(concordant, na.rm = TRUE)),
                     by = .(pathway, cell_type)]
hm_stab[, direction_concordance := NULL]
hm_stab <- merge(hm_stab, hm_dir_sum, by = c("pathway", "cell_type"), all.x = TRUE)
hm_stab[, pathway_label := clean_pathway(pathway)]
hm_stab_main <- hm_stab[pathway_label %in% main_pathways]
hm_stab_main[, pathway_label := factor(pathway_label, levels = rev(main_pathways))]
hm_stab_main[, cell_label := clean_cell(cell_type)]
p_hm_stab <- ggplot(hm_stab_main, aes(cell_label, pathway_label,
                                      fill = direction_concordance)) +
  geom_tile(colour = "black", linewidth = 0.28, width = 0.92, height = 0.92) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * direction_concordance)), size = 2.9) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1)) +
  labs(title = "Pathway stability",
       subtitle = "Direction retained in seven donor omissions",
       x = NULL, y = NULL, fill = "Concordance") +
  theme_heat(10.0) + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           # Keep patchwork's D tag away from this short title.
                           plot.title = element_text(face = "bold",
                             margin = margin(l = 44, b = 8)),
                           legend.position = "bottom")

fig4 <- (p_volcano | p_conc) /
  (p_hm_main | p_hm_stab) +
  plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(title = "Donor-level transcriptomic and pathway validation",
                  tag_levels = "A") & panel_tag_theme
save_dual("Figure4_transcriptomic_pathway_validation", fig4, 20.0, 12.2)
copy_plotted_table(de_full, "Figure4A_full_model_DE.csv.gz")
copy_plotted_table(conc, "Figure4B_DE_concordance.csv.gz")
copy_plotted_table(hm, "Figure4C_Hallmark_all50.csv")
copy_plotted_table(hm_stab, "Figure4D_Hallmark_LODO_stability.csv")
mark("Figure 4", "DONE", "DE, concordance, Hallmark and LODO stability")

# ---------------- SUPPLEMENTARY 1: ALL 50 HALLMARK ---------------------------

hm_order <- hm[, .(mean_NES = mean(NES, na.rm = TRUE)), by = pathway_label]
setorder(hm_order, mean_NES)
if (nrow(hm_order) != 50L || uniqueN(hm$cell_type) != 2L ||
    any(hm[, .N, by = cell_type]$N != 50L))
  stop("Supplementary Figure 1 requires exactly 50 pathways for each target cell type")
hm[, cell_label := clean_cell(cell_type)]
hm_limits <- max(abs(hm$NES), na.rm = TRUE)
hallmark_parts <- list()
for (part in c("A", "B")) {
  ids <- if (part == "A") seq_len(25L) else 26:50
  these <- hm_order$pathway_label[ids]
  block <- copy(hm[pathway_label %in% these])
  block[, pathway_factor := factor(pathway_label, levels = these)]
  p_hm25 <- ggplot(block, aes(cell_label, pathway_factor)) +
    geom_tile(aes(fill = NES), colour = "#454545", linewidth = 0.24,
              width = 0.94, height = 0.94) +
    geom_point(data = block[significance == TRUE], shape = 8, size = 2.1,
               colour = "black", stroke = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, limits = c(-hm_limits, hm_limits)) +
    labs(title = paste0("Hallmark pathways (", part, "; 25 of 50)"),
         subtitle = "Asterisk: within-cell-type FDR < 0.05; identical NES scale across A/B",
         x = NULL, y = NULL, fill = "NES") +
    theme_heat(10.5) +
    theme(axis.text.y = element_text(size = 10.4),
          axis.text.x = element_text(size = 10.8))
  hallmark_parts[[part]] <- p_hm25
  save_dual(paste0("Supplementary_Figure1", part, "_Hallmark_25"),
            p_hm25, 11.2, 9.0, DIR_SUPP)
  copy_plotted_table(block, paste0("Supplementary_Figure1", part,
                                   "_Hallmark_25_source.csv"))
}
# A/B remain full-size pages; a side-by-side composite is supplied for overview.
hallmark_overview <- wrap_plots(hallmark_parts, ncol = 2) +
  plot_annotation(title = "Hallmark pathways across both target cell types",
                  tag_levels = "A") & panel_tag_theme
save_dual("Supplementary_Figure1_AB_combined_overview",
          hallmark_overview, 23.5, 9.4, DIR_SUPP)
mark("Supplementary Figure 1", "DONE", "A/B full-size pages plus combined overview; 50 x 2 evaluated")

# ---------------- SUPPLEMENTARY 2: SAMPLE QC ---------------------------------

qc_panels <- list()
if (!is.null(sample_qc)) {
  need(sample_qc, c("donor_id", "group", "cells_before", "cells_after_adaptive_QC",
                    "singlet", "doublet_rate", "retention_percent"), "Sample QC")
  qc_long <- melt(sample_qc,
                  id.vars = c("donor_id", "group"),
                  measure.vars = c("cells_before", "cells_after_adaptive_QC", "singlet"),
                  variable.name = "stage", value.name = "cells")
  qc_long[, stage := factor(stage,
                            levels = c("cells_before", "cells_after_adaptive_QC", "singlet"),
                            labels = c("Raw", "After adaptive QC", "After doublet removal"))]
  p_qc_cells <- ggplot(qc_long, aes(donor_id, cells / 1000, fill = group)) +
    geom_col(width = 0.68, colour = COL_DARK, linewidth = 0.25) +
    facet_wrap(~stage, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = GROUP_COLS) +
    labs(title = "Cells retained through QC", x = NULL, y = "Cells (thousands)",
         fill = NULL) +
    theme_pub(10.8) + theme(axis.text.x = element_text(angle = 35, hjust = 1,
                                                       size = 10),
                           legend.position = "bottom")
  p_qc_rate <- ggplot(sample_qc, aes(retention_percent, donor_id, colour = group)) +
    geom_segment(aes(x = 0, xend = retention_percent, yend = donor_id),
                 colour = "#DDDDDD", linewidth = 0.48) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("%.1f%%", retention_percent)), hjust = -0.25,
              colour = COL_DARK, size = 3.2) +
    scale_colour_manual(values = GROUP_COLS) +
    scale_x_continuous(limits = c(0, max(sample_qc$retention_percent) * 1.12)) +
    labs(title = "Adaptive-QC retention", x = "Retained cells (%)", y = NULL,
         colour = NULL) +
    theme_pub(8.8) + theme(legend.position = "bottom")
  p_doublet <- ggplot(sample_qc, aes(doublet_rate * 100, donor_id, colour = group)) +
    geom_segment(aes(x = 0, xend = doublet_rate * 100, yend = donor_id),
                 colour = "#DDDDDD", linewidth = 0.48) +
    geom_point(size = 3) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * doublet_rate)), hjust = -0.25,
              colour = COL_DARK, size = 3.2) +
    scale_colour_manual(values = GROUP_COLS) +
    scale_x_continuous(limits = c(0, max(sample_qc$doublet_rate * 100) * 1.18)) +
    labs(title = "scDblFinder calls", x = "Doublet rate", y = NULL,
         colour = NULL) +
    theme_pub(8.8) + theme(legend.position = "bottom")
  qc_panels <- c(qc_panels, list(p_qc_cells, p_qc_rate, p_doublet))
}

if (!is.null(lib_qc)) {
  need(lib_qc, c("donor_id", "cell_type", "group", "library_size", "detected_genes"),
       "Pseudobulk library QC")
  lib_qc[, cell_label := clean_cell(cell_type)]
  p_lib <- ggplot(lib_qc, aes(library_size / 1e6, detected_genes,
                              colour = group, label = donor_id)) +
    geom_point(size = 2.8) +
    ggrepel::geom_text_repel(size = 3.0, max.overlaps = Inf,
                             min.segment.length = 0) +
    facet_wrap(~cell_label, nrow = 1, scales = "free") +
    scale_colour_manual(values = GROUP_COLS) +
    labs(title = "Pseudobulk library quality", x = "Library size (millions)",
         y = "Detected genes", colour = NULL) +
    theme_pub(8.8) + theme(legend.position = "bottom")
  qc_panels <- c(qc_panels, list(p_lib))
}

if (length(qc_panels)) {
  s2 <- wrap_plots(qc_panels, ncol = 2, guides = "collect") +
    plot_annotation(title = "Per-donor single-cell and pseudobulk quality control",
                    tag_levels = "A") & panel_tag_theme
  save_dual("Supplementary_Figure2_sample_QC", s2, 17.0,
            ifelse(length(qc_panels) >= 4, 11.5, 7.2), DIR_SUPP)
  mark("Supplementary Figure 2", "DONE", paste(length(qc_panels), "QC panels"))
} else {
  mark("Supplementary Figure 2", "SKIPPED", "No sample-level QC tables")
}

# ---------------- SUPPLEMENTARY 3: COMPOSITION -------------------------------

if (!is.null(comp)) {
  need(comp, c("donor_id", "cell_type", "group", "proportion"), "Composition table")
  comp[, donor_id := factor(donor_id,
                            levels = c("H1", "H2", "H3", "M0707_01", "M0709_02",
                                       "M0709_03", "M210715_01"))]
  ct_levels <- comp[, .(total = sum(proportion)), by = cell_type][order(-total)]$cell_type
  comp[, cell_type := factor(cell_type, levels = ct_levels)]
  ct_cols <- setNames(grDevices::hcl.colors(length(ct_levels), "Spectral", rev = TRUE),
                      ct_levels)
  p_comp <- ggplot(comp, aes(donor_id, proportion * 100, fill = cell_type)) +
    geom_col(width = 0.76, colour = "white", linewidth = 0.18) +
    scale_fill_manual(values = ct_cols, labels = clean_cell) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(title = "Cell-type composition across donors",
         subtitle = "Descriptive composition after QC and doublet removal",
         x = NULL, y = "Cells (%)", fill = "Cell type") +
    theme_pub(9.0) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.3),
          legend.position = "right")
  save_dual("Supplementary_Figure3_cell_composition", p_comp, 14.2, 7.8, DIR_SUPP)
  mark("Supplementary Figure 3", "DONE", "Per-donor cell-type composition")
} else {
  mark("Supplementary Figure 3", "SKIPPED", "Composition table unavailable")
}

# ---------------- SUPPLEMENTARY 4: INSTRUMENT AUDIT --------------------------

need(inst, c("f_stat", "allele_compatible", "alignment"), "Instrument audit")
inst[, f_stat := as.numeric(f_stat)]
p_f <- ggplot(inst[is.finite(f_stat)], aes(f_stat)) +
  geom_histogram(bins = 65, fill = "#4C78A8", colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = 10, linetype = 2, colour = COL_MG, linewidth = 0.65) +
  scale_x_log10() +
  labs(title = "Instrument-strength distribution", x = "F statistic (log10 scale)",
       y = "Gene-cell tests") + theme_pub(9.0)
align <- inst[, .N, by = .(allele_compatible, alignment)]
align[, label := paste0("Compatible: ", allele_compatible, "\n", alignment)]
p_align <- ggplot(align, aes(reorder(label, N), N, fill = allele_compatible)) +
  geom_col(width = 0.65, colour = COL_DARK, linewidth = 0.25) +
  geom_text(aes(label = scales::comma(N)), hjust = -0.15, size = 3) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(`TRUE` = "#4DAF4A", `FALSE` = COL_MG)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Allele-alignment audit", x = NULL, y = "Rows", fill = NULL) +
  theme_pub(9.0) + theme(legend.position = "none")
s4 <- p_f | p_align + plot_annotation(tag_levels = "A") & panel_tag_theme
save_dual("Supplementary_Figure4_instrument_alignment", s4, 16.0, 6.5, DIR_SUPP)
mark("Supplementary Figure 4", "DONE", "Instrument strength and allele alignment")

# ---------------- SUPPLEMENTARY 5: ROBUSTNESS/CALIBRATION --------------------

robust_panels <- list()
if (!is.null(robust)) {
  rlong <- melt(robust,
    id.vars = c("cell_type", "full_primary_genes"),
    measure.vars = c("direction_stable_fraction", "primary_5of7_fraction",
                     "primary_all7_fraction"),
    variable.name = "criterion", value.name = "fraction")
  rlong[, criterion := factor(criterion,
    levels = c("direction_stable_fraction", "primary_5of7_fraction", "primary_all7_fraction"),
    labels = c("Direction stable 7/7", "Primary >=5/7", "Primary 7/7"))]
  rlong[, cell_label := clean_cell(cell_type)]
  p_robust <- ggplot(rlong, aes(criterion, fraction * 100, fill = cell_type)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62,
             colour = COL_DARK, linewidth = 0.24) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * fraction)),
              position = position_dodge(width = 0.72), vjust = -0.35, size = 3.3) +
    scale_fill_manual(values = TARGET_COLS, labels = clean_cell) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Genome-wide donor-omission robustness", x = NULL,
         y = "Full-model primary genes (%)", fill = NULL) +
    theme_pub(10.8) + theme(axis.text.x = element_text(angle = 22, hjust = 1,
                                                       size = 10),
                           legend.position = "bottom")
  robust_panels <- c(robust_panels, list(p_robust))
}
if (!is.null(exact_all)) {
  exact_plot <- exact_all[is.finite(parametric_PValue) & is.finite(exact_P_abs_logFC)]
  exact_plot[, cell_label := clean_cell(cell_type)]
  p_exact <- ggplot(exact_plot,
                    aes(-log10(pmax(parametric_PValue, 1e-300)),
                        -log10(pmax(exact_P_abs_logFC, 1e-300)), colour = cell_type)) +
    geom_point(size = 0.55, alpha = 0.42, position = position_jitter(height = 0.015)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
    scale_colour_manual(values = TARGET_COLS, labels = clean_cell) +
    labs(title = "Parametric versus exact donor-label evidence",
         x = expression(-log[10](parametric~P)),
         y = expression(-log[10](exact~P)), colour = NULL) +
    theme_pub(8.8) + theme(legend.position = "bottom")
  robust_panels <- c(robust_panels, list(p_exact))
}
if (length(robust_panels)) {
  s5 <- wrap_plots(robust_panels, ncol = 2) +
    plot_annotation(tag_levels = "A") & panel_tag_theme
  save_dual("Supplementary_Figure5_donor_robustness_calibration", s5,
            16.0, 6.4, DIR_SUPP)
  mark("Supplementary Figure 5", "DONE", paste(length(robust_panels), "robustness panels"))
} else {
  mark("Supplementary Figure 5", "SKIPPED", "LODO/exact tables unavailable")
}

# ---------------- SUPPLEMENTARY 6: HIERARCHY + CLAIM AUDITS ------------------

audit_panels <- list()
if (!is.null(hier)) {
  need(hier, c("analysis_unit", "tested_units", "FDR005_units", "discovery_fraction"),
       "Hierarchical testing")
  hier[, label := paste0(FDR005_units, "/", tested_units)]
  p_hier <- ggplot(hier, aes(discovery_fraction * 100,
                             reorder(analysis_unit, discovery_fraction),
                             fill = analysis_unit)) +
    geom_col(width = 0.62, colour = COL_DARK, linewidth = 0.25) +
    geom_text(aes(label = label), hjust = -0.15, size = 3.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = c("#C95B5B", "#76539A", "#3E7DAA")) +
    labs(title = "Discoveries by analysis unit",
         subtitle = "1-Mb groups are physical intervals",
         x = "FDR < 0.05 discovery fraction (%)", y = NULL) +
    theme_pub(8.8) + theme(legend.position = "none")
  audit_panels <- c(audit_panels, list(p_hier))
}
if (!is.null(phe)) {
  p_phe <- ggplot(phe, aes(distinct_SNP_trait_pairs,
                           reorder(category, distinct_SNP_trait_pairs))) +
    geom_col(width = 0.64, fill = "#3E78A8", colour = COL_DARK, linewidth = 0.24) +
    geom_text(aes(label = distinct_SNP_trait_pairs), hjust = -0.18, size = 3.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(title = "PheWAS descriptive coverage", subtitle = "Not an enrichment test",
         x = "Distinct SNP-trait pairs", y = NULL) + theme_pub(8.8)
  audit_panels <- c(audit_panels, list(p_phe))
}
if (!is.null(dice)) {
  p_dice <- ggplot(dice, aes(N, reorder(audit_class, N))) +
    geom_col(width = 0.64, fill = COL_GOLD, colour = COL_DARK, linewidth = 0.24) +
    geom_text(aes(label = N), hjust = -0.18, size = 3.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Historical DICE claim audit",
         subtitle = "Selected rows; denominator unavailable",
         x = "Historical rows", y = NULL) + theme_pub(8.8)
  audit_panels <- c(audit_panels, list(p_dice))
}
if (length(audit_panels)) {
  s6 <- wrap_plots(audit_panels, ncol = min(3L, length(audit_panels))) +
    plot_annotation(tag_levels = "A") & panel_tag_theme
  save_dual("Supplementary_Figure6_dependency_and_claim_audits", s6,
            18.0, 7.2, DIR_SUPP)
  mark("Supplementary Figure 6", "DONE", paste(length(audit_panels), "audit panels"))
} else {
  mark("Supplementary Figure 6", "SKIPPED", "Audit tables unavailable")
}

# ---------------- SUPPLEMENTARY 7: CROSS-MODALITY MATRIX ---------------------

if (!is.null(evidence)) {
  need(evidence, c("gene", "cell_type", "MR_global_FDR005", "coloc_strong",
                   "pseudobulk_primary"), "Cross-modality evidence")
  ev <- evidence[MR_global_FDR005 == TRUE | coloc_strong == TRUE |
                   pseudobulk_primary == TRUE]
  ev[, support_n := rowSums(cbind(
    fifelse(is.na(MR_global_FDR005), FALSE, MR_global_FDR005),
    fifelse(is.na(coloc_strong), FALSE, coloc_strong),
    fifelse(is.na(pseudobulk_primary), FALSE, pseudobulk_primary)
  ))]
  setorder(ev, -support_n, MR_FDR_global, -PP_H4)
  ev <- ev[1:min(.N, 35L)]
  ev[, row_label := paste0(display_gene(gene), " (", clean_cell(cell_type), ")")]
  if (anyDuplicated(ev$row_label)) {
    ev[, row_label := fifelse(duplicated(row_label) |
      duplicated(row_label, fromLast = TRUE),
      paste0(row_label, " [", gene, "]"), row_label)]
  }
  if (anyDuplicated(ev$row_label)) stop("Repeated gene-cell labels in evidence matrix")
  ev_long <- melt(ev,
    id.vars = c("row_label", "support_n"),
    measure.vars = c("MR_global_FDR005", "coloc_strong", "pseudobulk_primary"),
    variable.name = "evidence", value.name = "supported")
  ev_long[is.na(supported), supported := FALSE]
  ev_long[, evidence := factor(evidence,
    levels = c("MR_global_FDR005", "coloc_strong", "pseudobulk_primary"),
    labels = c("MR global FDR < 0.05", "Strong colocalisation", "Pseudobulk primary"))]
  label_order <- unique(ev$row_label)
  parts <- list(A = label_order[seq_len(min(length(label_order), 18L))],
                B = if (length(label_order) > 18L)
                  label_order[19:length(label_order)] else character())
  evidence_parts <- list()
  for (part in names(parts)) {
    these <- parts[[part]]
    if (!length(these)) next
    dd <- copy(ev_long[as.character(row_label) %in% these])
    dd[, row_label := factor(as.character(row_label), levels = rev(these))]
    p_ev <- ggplot(dd, aes(evidence, row_label, fill = supported)) +
      geom_tile(colour = "black", linewidth = 0.28, width = 0.94, height = 0.94) +
      scale_fill_manual(values = c(`FALSE` = "#F1F1F1", `TRUE` = COL_PURPLE)) +
      coord_fixed(ratio = 0.25) +
      labs(title = paste0("Cross-modality evidence (", part, ")"),
           subtitle = "Gene-cell pairs ranked by convergent evidence",
           x = NULL, y = NULL, fill = "Supported") +
      theme_heat(11.0) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10.5),
            axis.text.y = element_text(size = 10.0))
    evidence_parts[[part]] <- p_ev
    save_dual(paste0("Supplementary_Figure7", part,
                     "_cross_modality_evidence"), p_ev, 8.8,
              max(7.3, 0.36 * length(these) + 1.7), DIR_SUPP)
  }
  if (length(evidence_parts) == 2L) {
    evidence_overview <- wrap_plots(evidence_parts, ncol = 2) +
      plot_annotation(title = "Cross-modality evidence across gene-cell pairs",
                      tag_levels = "A") & panel_tag_theme
    save_dual("Supplementary_Figure7_AB_combined_overview",
              evidence_overview, 19.0, 9.2, DIR_SUPP)
  }
  copy_plotted_table(ev, "Supplementary_Figure7_plot_data.csv")
  mark("Supplementary Figure 7", "DONE", paste(nrow(ev),
      "gene-cell pairs split across A/B"))
} else {
  mark("Supplementary Figure 7", "SKIPPED", "Cross-modality evidence unavailable")
}

# ---------------- INDIVIDUAL PANELS AND REGIONAL OVERVIEW -------------------

# Each multi-panel figure above is kept intact in 01_main_figures or
# 02_supplementary_figures. Export its component plots without panel letters.
panel_specs <- list(
  list("Figure1_A_postMHC_counts", p_post, 8.5, 5.5),
  list("Figure1_B_eQTL_rescaling", p_rescale, 8.5, 5.5),
  list("Figure1_C_MR_coloc", p_joint, 8.5, 6.7),
  list("Figure1_D_forest", p_forest, 10.0, 7.3),
  list("Figure2_A_SESN3_donor_expression", p_expr, 12.0, 5.4),
  list("Figure2_B_SESN3_LODO", p_lodo, 11.0, 6.6),
  list("Figure2_C_exact_permutation", p_perm, 11.0, 6.6),
  list("Figure3_A_reviewed_UMAP", p_umap_cell, 9.0, 6.8),
  list("Figure3_B_disease_UMAP", p_umap_group, 9.0, 6.8),
  list("Figure3_C_marker_dotplot", p_marker, 12.5, 8.0),
  list("Figure3_D_donor_PCA", p_pb_pca, 10.0, 6.8),
  list("Figure4_A_pseudobulk_volcano", p_volcano, 13.0, 6.0),
  list("Figure4_B_DE_concordance", p_conc, 8.0, 6.0),
  list("Figure4_C_Hallmark_activity", p_hm_main, 9.5, 9.3),
  list("Figure4_D_Hallmark_stability", p_hm_stab, 9.5, 9.3),
  list("Supplementary_MHC_UNRECONCILED", p_mhc, 12.5, 5.2),
  list("Supplementary1_A_Hallmark_25", hallmark_parts[["A"]], 11.2, 9.0),
  list("Supplementary1_B_Hallmark_25", hallmark_parts[["B"]], 11.2, 9.0),
  list("Supplementary2_A_retained_cells", if (exists("p_qc_cells")) p_qc_cells else NULL, 9.0, 7.5),
  list("Supplementary2_B_adaptive_QC", if (exists("p_qc_rate")) p_qc_rate else NULL, 8.0, 5.5),
  list("Supplementary2_C_doublets", if (exists("p_doublet")) p_doublet else NULL, 8.0, 5.5),
  list("Supplementary2_D_library_QC", if (exists("p_lib")) p_lib else NULL, 9.5, 5.5),
  list("Supplementary3_donor_composition", if (exists("p_comp")) p_comp else NULL, 14.2, 7.8),
  list("Supplementary4_A_instrument_strength", p_f, 8.0, 6.0),
  list("Supplementary4_B_allele_alignment", p_align, 8.0, 6.0),
  list("Supplementary5_A_donor_omission", if (exists("p_robust")) p_robust else NULL, 9.0, 6.0),
  list("Supplementary5_B_exact_calibration", if (exists("p_exact")) p_exact else NULL, 9.0, 6.0),
  list("Supplementary6_A_discoveries", if (exists("p_hier")) p_hier else NULL, 9.0, 5.5),
  list("Supplementary6_B_PheWAS", if (exists("p_phe")) p_phe else NULL, 9.0, 5.5),
  list("Supplementary6_C_DICE", if (exists("p_dice")) p_dice else NULL, 9.0, 5.5),
  list("Supplementary7_A_cross_modality", if (exists("evidence_parts")) evidence_parts[["A"]] else NULL, 8.8, 8.0),
  list("Supplementary7_B_cross_modality", if (exists("evidence_parts")) evidence_parts[["B"]] else NULL, 8.8, 8.0)
)
panel_inventory <- rbindlist(lapply(panel_specs, function(spec) {
  ok <- inherits(spec[[2L]], "ggplot")
  if (ok) save_dual(spec[[1L]], spec[[2L]] +
       labs(subtitle = NULL, caption = NULL), spec[[3L]], spec[[4L]], DIR_PANEL)
  data.table(panel = spec[[1L]], exported = ok,
             size_in = paste(spec[[3L]], spec[[4L]], sep = " x "))
}))
fwrite(panel_inventory, file.path(OUTPUT_ROOT, "V4_INDIVIDUAL_PANEL_INDEX.csv"))

# The corrected core supplied only finished regional graphics. PNGs can be
# assembled for an overview, while their original PDFs remain separate vector
# originals. The combined PDF embeds PNG and is identified as raster below.
if (!is.null(reg4) && !is.null(reg8) && file.exists(reg4$png) &&
    file.exists(reg8$png)) {
  region_image <- function(path) {
    ggplot() + annotation_custom(grid::rasterGrob(
      png::readPNG(path), width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"), interpolate = TRUE),
      xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  }
  regional_overview <- region_image(reg4$png) | region_image(reg8$png)
  save_dual("Figure2_regional_CD4NC_CD8NC_combined_RASTER_overview",
            regional_overview, 21.0, 7.2, DIR_MAIN)
}

# ---------------- CURATED TABLES AND FIGURE-TO-DATA INDEX -------------------

table_sources <- list(
  Table01_postMHC_discovery_counts = post_counts,
  Table02_MHC_before_after_UNRECONCILED = mhc,
  Table03_MR_coloc_evidence = joint,
  Table04_SESN3_donor_expression = expr,
  Table05_SESN3_LODO = lodo,
  Table06_SESN3_exact_inference = sesn3_exact,
  Table07_pseudobulk_DE = de,
  Table08_Hallmark_50x2 = hallmark,
  Table09_sample_QC = sample_qc,
  Table10_per_donor_composition = comp,
  Table11_hierarchical_sensitivity = hier,
  Table12_PheWAS_descriptive_only = phe,
  Table13_DICE_historical_only = dice,
  Table14_instrument_alignment = inst,
  Table15_cross_modality_evidence = evidence
)
table_index <- rbindlist(lapply(names(table_sources), function(nm) {
  obj <- table_sources[[nm]]
  if (is.null(obj)) return(data.table(table = nm, status = "MISSING_OPTIONAL",
                                      file = "", rows = NA_integer_))
  target <- paste0(nm, if (nrow(obj) > 50000L) ".csv.gz" else ".csv")
  fwrite(obj, file.path(DIR_TABLE, target))
  data.table(table = nm, status = "COPIED_FROM_LOCKED_V2",
             file = paste0("07_analysis_tables/", target), rows = nrow(obj))
}), fill = TRUE)
fwrite(table_index, file.path(OUTPUT_ROOT, "V4_ANALYSIS_TABLE_INDEX.csv"))
figure_to_data <- data.table(
  figure = c("Figure1", "Figure2", "Figure2 regional", "Figure3",
             "Figure4", "Supplementary1", "Supplementary2", "Supplementary3",
             "Supplementary4", "Supplementary5", "Supplementary6",
             "Supplementary7", "Supplementary MHC audit"),
  source_prefix = c("Figure1A-D", "Figure2A-C", "Legacy regional originals",
                    "Figure3_UMAP", "Figure4A-D", "Supplementary_Figure1A/B",
                    "Table09_sample_QC", "Table10_per_donor_composition",
                    "Table14_instrument_alignment", "Table05_SESN3_LODO",
                    "Table12_PheWAS and Table13_DICE", "Supplementary_Figure7",
                    "Table02_MHC_before_after_UNRECONCILED"),
  caveat = c("Physical 1-Mb groups are descriptive", "Exact P(F) differs by cell type",
             "Regional source-data table unavailable; composed overview is raster",
             "The UMAP coordinates are integrated; DE uses raw RNA counts",
             "The two cell types share donors", "All 50 sets per cell type; shared colour scale",
             "Seven donors", "Descriptive, not a statistical composition test",
             "F statistic and allele alignment", "LODO is sensitivity, not replication",
             "No PheWAS enrichment or full DICE replication denominator",
             "Gene names require HGNC review", "Submitted baseline has not been reconciled")
)
fwrite(figure_to_data, file.path(OUTPUT_ROOT, "V4_FIGURE_TO_TABLE_MAP.csv"))

# ---------------- FIGURE INDEX AND REVIEWER COVERAGE -------------------------

used_symbols <- unique(c(as.character(joint_lab$gene), as.character(forest$gene),
                         as.character(lab_de$gene), as.character(conc_lab$gene),
                         if (exists("ev", inherits = FALSE)) as.character(ev$gene)
                         else character()))
used_symbols <- used_symbols[!is.na(used_symbols) & nzchar(used_symbols)]
gene_label_audit <- data.table(old_symbol = sort(used_symbols))
gene_label_audit[, `:=`(display_symbol = display_gene(old_symbol),
  hgnc_id = hgnc_mapping$hgnc_id[match(old_symbol, hgnc_mapping$old_symbol)],
  label_status = fifelse(old_symbol %in% hgnc_mapping$old_symbol,
                        "CLIENT_MAPPING_USED", "NO_CLIENT_MAPPING_FOR_THIS_SYMBOL"))]
fwrite(gene_label_audit, file.path(DIR_TAB, "HGNC_all_plotted_gene_labels_audit.csv"))

figure_index <- data.table(
  figure = c(
    "Figure 1", "Figure 2", "Figure 2 regional CD4 NC", "Figure 2 regional CD8 NC",
    "Figure 3", "Figure 4", "Supplementary Figure 1",
    "Supplementary Figure 2", "Supplementary Figure 3",
    "Supplementary Figure 4", "Supplementary Figure 5",
    "Supplementary Figure 6", "Supplementary Figure 7"
  ),
  content = c(
    "A: post-MHC discoveries; B: eQTL effect correction; C: MR-coloc; D: compact forest; pre-MHC audit in supplement",
    "SESN3 donor expression; LODO; exact permutation",
    "SESN3 CD4 NC regional colocalisation",
    "SESN3 CD8 NC regional colocalisation",
    "Reviewed UMAP; group distribution; marker validation; donor PCA",
    "Pseudobulk DE; cross-cell concordance; Hallmark; Hallmark LODO",
    "All 50 Hallmark pathways in two readable 25-pathway parts (A/B)",
    "Per-donor single-cell and pseudobulk QC",
    "Cell-type composition across donors",
    "Instrument strength and allele alignment",
    "Genome-wide LODO and exact donor-label calibration",
    "Locus, PheWAS and DICE claim audits",
    "Cross-modality evidence matrix in two readable parts (A/B)"
  ),
  recommended_location = c(
    "Main", "Main", "Main or supplementary", "Main or supplementary",
    "Main", "Main", rep("Supplementary", 7)
  )
)
fwrite(figure_index, file.path(OUTPUT_ROOT, "V4_FIGURE_INDEX.csv"))

coverage <- data.table(
  reviewer_requirement = c(
    "Extended-MHC exclusion", "Corrected eQTL effect size", "MR and colocalisation",
    "Independent SNP/locus reporting", "SESN3 regional colocalisation",
    "SESN3 donor-level expression", "Single-cell pseudobulk and sample QC",
    "All 50 Hallmark pathways and FDR", "PheWAS supplementary analysis",
    "DICE complete denominator", "qPCR composition adjustment",
    "Strict ancestry-matched LD clumping"
  ),
  V3_location = c(
    "Supplementary MHC baseline audit", "Figure 1B", "Figure 1C-D",
    "Supplementary Figure 6A; physical loci only",
    "Figure 2 regional CD4 NC and CD8 NC", "Figure 2",
    "Figure 3 and Supplementary Figure 2", "Figure 4 and Supplementary Figure 1",
    "Supplementary Figure 6B; descriptive only",
    "Supplementary Figure 6C; incomplete historical denominator",
    "Not completed: participant-level composition data unavailable",
    "Not completed: ancestry-matched LD reference unavailable"
  ),
  status = c(
    rep("COVERED", 9), "LIMITED_BY_INPUT", "BLOCKED_MISSING_DATA",
    "BLOCKED_MISSING_DATA"
  )
)
fwrite(coverage, file.path(OUTPUT_ROOT, "V4_REVIEWER_REQUIREMENT_COVERAGE.csv"))

# Every letter on a combined figure has a dedicated explanation suitable for
# a manuscript legend. These legends are a draft: confirm numbering and terms
# against the submitted manuscript before publication.
figure_legends <- c(
  "# Draft panel legends for the redrawn figures",
  "",
  "## Figure 1. Corrected post-MHC genetic analysis",
  "A. Numbers of global-BH-significant gene-cell pairs, distinct lead variants and descriptive 1-Mb physical groups in the post-MHC analysis. Physical groups are not LD-independent loci.",
  "B. Published OneK1K Spearman rho and the corresponding rescaled per-allele eQTL effect; point colour represents the conversion factor.",
  "C. Corrected MR global BH q-values and regional colocalisation posterior probabilities (PP.H4) for evaluable gene-cell pairs; dashed lines mark the indicated thresholds. Each selected label appears once per gene.",
  "D. Corrected MR odds ratios and 95% confidence intervals for the top strong-colocalisation results; full results remain in the source table.",
  "",
  "## Figure 2. SESN3 donor-level robustness",
  "A. Normalised pseudobulk SESN3 expression for each donor in CD4 NC and CD8 NC populations, after the annotation correction.",
  "B. Full and seven leave-one-donor-out edgeR model estimates for SESN3; all points share donors from the same cohort.",
  "C. Exact permutation distribution of the edgeR F statistic over all 35 allocations of three HC and four MG donors. F-statistic exact P is 0.0571 for CD4 NC and 0.0286 for CD8 NC. The distinct absolute-log2FC test yields 0.0571 for both; within-cell-type BH-adjusted exact P values are displayed separately in each panel.",
  "",
  "## Figure 3. Single-cell review",
  "A. Reviewed cell-type annotation of the integrated UMAP. B. HC/MG distribution on the same embedding. C. Canonical marker expression by reviewed cell type. D. Donor-level pseudobulk PCA (or cell-composition PCA when pseudobulk coordinates are unavailable).",
  "",
  "## Figure 4. Transcriptomic and pathway validation",
  "A. Donor-level pseudobulk differential expression with raw P < 0.05 and |log2FC| >= 0.5 used only as the displayed threshold; FDR remains available in source tables.",
  "B. Gene-level CD4 NC and CD8 NC log2FC concordance in the same donors; this is not an independent replication cohort.",
  "C. Top Hallmark pathways by the within-cell-type FDR and NES; all 50 evaluated pathways per target population are shown in Supplementary Figure 1A/B.",
  "D. Proportion of seven donor-omission models retaining the full-model NES direction for the selected pathways.",
  "",
  "## Supplementary Figure 1",
  "A and B. Twenty-five distinct Hallmark pathways per page. All 50 pathways are represented for both populations on the same NES scale; stars mark within-cell-type FDR < 0.05.",
  "",
  "## Supplementary Figure 2. Per-donor QC",
  "A. Cell counts before filtering, after adaptive QC and after doublet removal. B. Adaptive-QC retention. C. scDblFinder doublet rate. D. Pseudobulk library size and detected genes when the library QC source is available.",
  "",
  "## Supplementary Figure 3. Cell composition",
  "Per-donor percentages for reviewed cell types after QC and doublet removal; this is a descriptive composition chart.",
  "",
  "## Supplementary Figure 4. Instrument audit",
  "A. Distribution of corrected instrument F statistics, with F = 10 as a guide. B. Allele-compatibility and alignment row counts.",
  "",
  "## Supplementary Figure 5. Donor sensitivity",
  "A. Genome-wide direction and primary-threshold stability across leave-one-donor-out models. B. Parametric versus exact donor-label P values.",
  "",
  "## Supplementary Figure 6. Evidence and claim audits",
  "A. BH discoveries at distinct units of analysis; physical 1-Mb intervals do not establish LD independence. B. Descriptive PheWAS SNP-trait category counts without enrichment inference. C. Historical DICE records from a selected table; replication denominator unavailable.",
  "",
  "## Supplementary Figure 7. Cross-modality evidence",
  "A and B. Support matrix for the selected gene-cell pairs split across two readable pages; rows are pairs and columns represent MR, colocalisation and donor-level pseudobulk evidence."
)
writeLines(figure_legends, file.path(OUTPUT_ROOT, "V4_DRAFT_PANEL_LEGENDS.md"))

visual_qa <- data.table(
  figure = c("Figure 1", "Figure 2", "SESN3 regional CD4 NC",
             "SESN3 regional CD8 NC", "Figure 3", "Figure 4",
             "Supplementary Figure 1A", "Supplementary Figure 1B",
             paste0("Supplementary Figure ", 2:6),
             "Supplementary Figure 7A", "Supplementary Figure 7B",
             "Supplementary MHC audit"),
  qa_state = "PENDING_REVIEW_OF_NEWLY_RENDERED_PDF",
  examine = c(
    "A-D panel titles; MR gene labels; all forest labels and legend",
    "Donor IDs and exact P labels; 0.0571/0.0286 match the F-statistic source table",
    "Retained existing core-pipeline PDF: check legend/title/axes; source data not supplied for a full redraw",
    "Retained existing core-pipeline PDF: check legend/title/axes; source data not supplied for a full redraw",
    "UMAP annotation and 26 italic marker names; no axis or legend clipping",
    "Volcano labels; pathway labels; all four panel headings and legends",
    "First 25 Hallmark pathways; shared NES scale and legible gene-set labels",
    "Remaining 25 Hallmark pathways; same NES scale as page A",
    "All seven donor IDs, QC fractions, and pseudobulk labels",
    "Seven donor IDs and all cell-type legend entries",
    "F-statistic and allele-alignment count labels",
    "Robustness bar headings, percentages, and exact-P axes",
    "Three panel headings, long phenotype labels, and DICE denominator caveat",
    "First 18 gene-cell labels; three column headings; square tiles",
    "Remaining gene-cell labels; identical palette and square tiles",
    "Only comparable before/after metrics; explicit baseline-reconciliation caveat"
  )
)
fwrite(visual_qa, file.path(OUTPUT_ROOT, "V4_REQUIRED_VISUAL_QA.csv"))

# ---------------- FINAL VALIDATION AND PACKAGE -------------------------------

fwrite(status, file.path(OUTPUT_ROOT, "V4_VISUALISATION_STATUS.csv"))
fwrite(source_audit, file.path(DIR_LOG, "V4_SOURCE_AUDIT.csv"))
capture.output(sessionInfo(), file = file.path(DIR_LOG, "sessionInfo.txt"))

# Preserve the exact script used when it can be resolved in Rscript or RStudio.
script_path <- NA_character_
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- sub("^--file=", "", script_arg[1L])
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                          error = function(e) NA_character_)
}
if (!is.na(script_path) && nzchar(script_path) && file.exists(script_path)) {
  file.copy(script_path,
            file.path(DIR_CODE, "MG_V4_complete_figures_and_tables.R"),
            overwrite = TRUE)
}

all_outputs <- list.files(OUTPUT_ROOT, recursive = TRUE, full.names = TRUE,
                          include.dirs = FALSE)
bad_names <- all_outputs[grepl("\\.(pdf|png)\\.[0-9]+$", all_outputs) |
                           grepl("/Y7_", norm_slash(all_outputs))]
if (length(bad_names)) stop("Duplicate/legacy figure names detected: ",
                            paste(basename(bad_names), collapse = ", "))

figure_files <- all_outputs[grepl("\\.(pdf|png)$", all_outputs, ignore.case = TRUE)]
if (!length(figure_files)) stop("No figures were produced")

manifest_file <- file.path(OUTPUT_ROOT, "V4_FILE_MANIFEST.csv")
zip_file <- file.path(OUTPUT_ROOT, if (nrow(hgnc_mapping))
  "MG_V4_publication_figures_tables_PENDING_VISUAL_QA.zip" else
  "MG_V4_publication_figures_tables_DRAFT_HGNC_PENDING.zip")
all_outputs <- all_outputs[
  norm_slash(all_outputs) != norm_slash(manifest_file) &
    norm_slash(all_outputs) != norm_slash(zip_file)
]
root_norm <- norm_slash(OUTPUT_ROOT)
manifest <- data.table(
  relative_path = substring(norm_slash(all_outputs), nchar(root_norm) + 2L),
  size_bytes = file.info(all_outputs)$size,
  md5 = unname(tools::md5sum(all_outputs))
)
setorder(manifest, relative_path)
fwrite(manifest, manifest_file)

writeLines(c(
  paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "No statistical model was rerun by this V3 visualisation script.",
  "Standalone panels have no letter tags; all multi-panel figures have combined exports.",
  "Grey plot subtitles/captions are removed from redrawn charts; their interpretation remains in panel legends and caveat tables.",
  "The two inherited regional plots retain their original layouts; their composed overview embeds PNG and is raster in the PDF container.",
  "PDF files are vector output and PNG files are 600 dpi.",
  paste0("DE display threshold: raw P < ", DE_P_THRESHOLD,
         " and |log2FC| >= ", DE_ABS_LOG2FC_THRESHOLD, "."),
  "All 50 Hallmark pathways are shown in Supplementary Figure 1A/B.",
  "Fig2C uses the exact F-statistic permutation P; absolute-log2FC P is separately reported.",
  paste("HGNC display mapping:", if (nrow(hgnc_mapping))
    "provided; verified ID mapping saved with figures" else
    "NOT PROVIDED; old symbol labels require checking before client delivery"),
  "MHC boundary analysis is post hoc; it filters existing MR tests, reports both retained original and recalculated BH q-values, and does not rerun coloc.",
  "Historical drug strings are review candidates only; MG indications and mechanism directions were not verified.",
  "The 7,852/396 pre-MHC reanalysis baseline is not reconciled to the submitted 7,598/365.",
  "PheWAS is descriptive, not enrichment.",
  "Historical DICE rows do not provide a complete denominator.",
  "Physical 1-Mb groups must not be described as LD-independent loci.",
  "qPCR composition adjustment and strict LD clumping remain blocked by missing inputs."
), file.path(OUTPUT_ROOT, "V4_README.txt"))

old_wd <- getwd()
zip_ok <- tryCatch({
  if (file.exists(zip_file)) unlink(zip_file, force = TRUE)
  setwd(OUTPUT_ROOT)
  utils::zip(
    zipfile = basename(zip_file),
    files = c("01_main_figures", "02_supplementary_figures", "03_source_tables",
              "04_logs", "05_code", "06_individual_panels", "07_analysis_tables",
              "V4_INDIVIDUAL_PANEL_INDEX.csv", "V4_ANALYSIS_TABLE_INDEX.csv",
              "V4_FIGURE_TO_TABLE_MAP.csv", "V4_FIGURE_INDEX.csv",
              "V4_REVIEWER_REQUIREMENT_COVERAGE.csv",
              "V4_VISUALISATION_STATUS.csv", "V4_FILE_MANIFEST.csv",
              "V4_README.txt", "V4_DRAFT_PANEL_LEGENDS.md",
              "V4_REQUIRED_VISUAL_QA.csv")
  )
  file.exists(zip_file)
}, error = function(e) {
  warning("ZIP creation failed, but all V3 outputs were saved: ", e$message)
  FALSE
}, finally = setwd(old_wd))

message("\nV4 figures and tables complete: ", OUTPUT_ROOT)
if (zip_ok) message("V3 package: ", zip_file)
message("Main figure files: ", length(list.files(DIR_MAIN)))
message("Supplementary figure files: ", length(list.files(DIR_SUPP)))
print(status)
