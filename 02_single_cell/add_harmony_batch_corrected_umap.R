suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript add_harmony_batch_corrected_umap.R <output_dir>")

output_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
results_dir <- file.path(output_dir, "results")
input_rds <- file.path(results_dir, "seurat_440_validation_refined.rds")
output_rds <- file.path(results_dir, "seurat_440_validation_refined_batch_corrected.rds")

message("Reading refined Seurat object: ", input_rds)
obj <- readRDS(input_rds)

if (!"pca" %in% Reductions(obj)) {
  stop("The input object does not contain PCA reduction.")
}
if (!"sample" %in% colnames(obj[[]])) {
  stop("The input object does not contain sample metadata.")
}

message("Running Harmony batch correction by sample...")
obj <- RunHarmony(
  object = obj,
  group.by.vars = "sample",
  reduction = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  project.dim = FALSE,
  verbose = TRUE
)

message("Computing neighbors, clusters, and UMAP on Harmony embedding...")
obj <- FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = 1:30,
  graph.name = "harmony_snn",
  verbose = FALSE
)
obj <- FindClusters(
  object = obj,
  graph.name = "harmony_snn",
  resolution = 0.6,
  cluster.name = "harmony_clusters",
  verbose = FALSE
)
obj <- RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony",
  reduction.key = "hUMAP_",
  verbose = FALSE
)

saveRDS(obj, output_rds)
message("Batch-corrected object written to: ", output_rds)
