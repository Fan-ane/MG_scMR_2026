suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript make_target_sample_level_summary.R <output_dir>")

output_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
results_dir <- file.path(output_dir, "results")
obj <- readRDS(file.path(results_dir, "seurat_440_validation_refined.rds"))
expr_mat <- GetAssayData(obj, assay = "RNA", slot = "data")

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

per_sample <- list()
summary_rows <- list()
for (i in seq_len(nrow(targets))) {
  gene <- targets$gene[i]
  ct <- targets$target_cell_type[i]
  cells <- colnames(obj)[obj$cell_type_detail_refined == ct]
  for (sample_id in sort(unique(obj$sample))) {
    scells <- cells[obj$sample[cells] == sample_id]
    group <- unique(as.character(obj$group[obj$sample == sample_id]))[1]
    vals <- if (gene %in% rownames(obj) && length(scells) > 0) as.numeric(expr_mat[gene, scells]) else numeric(0)
    per_sample[[length(per_sample) + 1]] <- data.frame(
      gene = gene,
      target_cell_type = ct,
      discovery_direction = targets$discovery_direction[i],
      sample = sample_id,
      group = group,
      n_cells = length(scells),
      mean_expr = ifelse(length(vals) > 0, mean(vals), NA_real_),
      pct_expressing = ifelse(length(vals) > 0, mean(vals > 0) * 100, NA_real_),
      stringsAsFactors = FALSE
    )
  }
}
per_sample <- bind_rows(per_sample)
write.csv(per_sample, file.path(results_dir, "target_gene_sample_level_values_refined.csv"), row.names = FALSE)

for (key in unique(paste(per_sample$gene, per_sample$target_cell_type, sep = "||"))) {
  dat <- per_sample[paste(per_sample$gene, per_sample$target_cell_type, sep = "||") == key, ]
  hc <- dat$mean_expr[dat$group == "HC"]
  mg <- dat$mean_expr[dat$group == "MG"]
  hc <- hc[!is.na(hc)]
  mg <- mg[!is.na(mg)]
  direction <- ifelse(length(hc) > 0 && length(mg) > 0 && mean(mg) > mean(hc), "up",
                      ifelse(length(hc) > 0 && length(mg) > 0 && mean(mg) < mean(hc), "down", "not_tested"))
  p <- ifelse(length(hc) >= 2 && length(mg) >= 2,
              tryCatch(t.test(mg, hc)$p.value, error = function(e) NA_real_),
              NA_real_)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    gene = dat$gene[1],
    target_cell_type = dat$target_cell_type[1],
    discovery_direction = dat$discovery_direction[1],
    expected_MG_direction = ifelse(dat$discovery_direction[1] == "risk", "up", "down"),
    HC_sample_mean = ifelse(length(hc) > 0, mean(hc), NA_real_),
    MG_sample_mean = ifelse(length(mg) > 0, mean(mg), NA_real_),
    observed_MG_direction_sample_mean = direction,
    sample_t_p = p,
    stringsAsFactors = FALSE
  )
}
summary_rows <- bind_rows(summary_rows) %>%
  mutate(direction_consistent_sample_mean = observed_MG_direction_sample_mean == expected_MG_direction)
write.csv(summary_rows, file.path(results_dir, "target_gene_sample_level_summary_refined.csv"), row.names = FALSE)
