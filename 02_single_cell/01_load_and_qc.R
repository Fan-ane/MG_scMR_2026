# ============================================================
# SESN3 Virtual Knockout — Full Analysis Pipeline
# Part 1: Load data, QC, integration, clustering, annotation
# ============================================================

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)

set.seed(123)

# ---- 1. Load 10X data ----

data_dir <- "c:/Users/zyf/Desktop/课题一生信文章/验证部分工作/scRNA_analysis"

# 3 HC samples
h1 <- Read10X(file.path(data_dir, "H1_matrix"))
h2 <- Read10X(file.path(data_dir, "H2_matrix"))
h3 <- Read10X(file.path(data_dir, "H3_matrix"))

# 4 MG samples
m1 <- Read10X(file.path(data_dir, "M0707-01_matrix_10X"))
m2 <- Read10X(file.path(data_dir, "M0709-02_matrix_10X"))
m3 <- Read10X(file.path(data_dir, "M0709-03_matrix_10X"))
m4 <- Read10X(file.path(data_dir, "M210715-01_matrix_10X"))

# Create Seurat objects
h1_obj <- CreateSeuratObject(h1, project = "HC1", min.cells = 3, min.features = 200)
h2_obj <- CreateSeuratObject(h2, project = "HC2", min.cells = 3, min.features = 200)
h3_obj <- CreateSeuratObject(h3, project = "HC3", min.cells = 3, min.features = 200)
m1_obj <- CreateSeuratObject(m1, project = "MG1", min.cells = 3, min.features = 200)
m2_obj <- CreateSeuratObject(m2, project = "MG2", min.cells = 3, min.features = 200)
m3_obj <- CreateSeuratObject(m3, project = "MG3", min.cells = 3, min.features = 200)
m4_obj <- CreateSeuratObject(m4, project = "MG4", min.cells = 3, min.features = 200)

# Add metadata
h1_obj$condition <- "HC"; h1_obj$sample <- "HC1"
h2_obj$condition <- "HC"; h2_obj$sample <- "HC2"
h3_obj$condition <- "HC"; h3_obj$sample <- "HC3"
m1_obj$condition <- "MG"; m1_obj$sample <- "MG1"
m2_obj$condition <- "MG"; m2_obj$sample <- "MG2"
m3_obj$condition <- "MG"; m3_obj$sample <- "MG3"
m4_obj$condition <- "MG"; m4_obj$sample <- "MG4"

# ---- 2. QC and filtering ----

# Calculate QC metrics
for (obj in list(h1_obj, h2_obj, h3_obj, m1_obj, m2_obj, m3_obj, m4_obj)) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
}

# Filter
filter_cells <- function(obj) {
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
}

h1_obj <- filter_cells(h1_obj)
h2_obj <- filter_cells(h2_obj)
h3_obj <- filter_cells(h3_obj)
m1_obj <- filter_cells(m1_obj)
m2_obj <- filter_cells(m2_obj)
m3_obj <- filter_cells(m3_obj)
m4_obj <- filter_cells(m4_obj)

# ---- 3. Merge and normalize ----

merged <- merge(h1_obj, y = c(h2_obj, h3_obj, m1_obj, m2_obj, m3_obj, m4_obj),
                add.cell.ids = c("HC1", "HC2", "HC3", "MG1", "MG2", "MG3", "MG4"))

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000)
merged <- ScaleData(merged, vars.to.regress = "percent.mt")
merged <- RunPCA(merged, features = VariableFeatures(merged))

# ---- 4. Integration (Harmony) ----

library(harmony)
merged <- RunHarmony(merged, group.by.vars = "sample", dims.use = 1:30)

# ---- 5. Clustering and UMAP ----

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.8)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)

# ---- 6. Cell type annotation ----

# Canonical markers for PBMC immune cells
canonical_markers <- list(
  "CD4 T cells" = c("CD3D", "CD4", "IL7R"),
  "CD8 T cells" = c("CD3D", "CD8A", "CD8B"),
  "NK cells" = c("NKG7", "GNLY", "KLRD1"),
  "B cells" = c("MS4A1", "CD79A", "CD79B"),
  "Monocytes" = c("CD14", "FCGR3A", "LYZ"),
  "Dendritic cells" = c("FCER1A", "CST3", "CLEC10A"),
  "Plasma cells" = c("MZB1", "SDC1", "JCHAIN")
)

# Score each cluster
DotPlot(merged, features = unlist(canonical_markers), group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Assign cell types based on marker expression
# (Manual assignment - adjust based on DotPlot results)
# This will be refined after examining the dot plot

# ---- 7. Save ----

saveRDS(merged, file.path(data_dir, "merged_seurat_annotated.rds"))
cat("Part 1 complete. Merged object saved.\n")
cat("Dimensions:", ncol(merged), "cells,", nrow(merged), "genes\n")
cat("Clusters:", length(unique(merged$seurat_clusters)), "\n")
