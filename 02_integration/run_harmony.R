#!/usr/bin/env Rscript
# Harmony integration of the part 1 atlas object: singlet/condition filtering,
# normalization, PCA, Harmony, UMAP, Leiden clustering, cluster 12 removal.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(reticulate)
})

# ---- CONFIG -----------------------------------------------------------------
in_rds  <- "output/scRNA_consolidation_part1.rds"
out_dir <- "output"
exclude_conditions <- c(
  "empty_well_or_discarded_patchclamp",
  "Other_Age",
  "Other_Parkinson",
  "Other_T2Dreversed"
)
n_harmony_dims <- 20
resolutions    <- c(0.05, 0.1, 0.2)
# -----------------------------------------------------------------------------

reticulate::import("leidenalg")  # fails early if leidenalg is missing (needed for algorithm = 4)
set.seed(123)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(in_rds)
obj <- subset(obj, subset = doublet_final == "Singlet")
obj <- subset(obj,
              subset = Condition %in% exclude_conditions,
              invert = TRUE)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")
)
obj <- RunPCA(obj, features = VariableFeatures(obj))

obj <- RunHarmony(
  obj,
  group.by.vars = c("orig.ident", "Platform"),
  reduction.use = "pca",
  max.iter.harmony = 50,
  dims.use = 1:50,
  lambda = NULL
)

# Use the Harmony dimensions with the largest SD, not simply 1:N
dims <- order(obj[["harmony"]]@stdev, decreasing = TRUE)[seq_len(n_harmony_dims)]

obj <- RunUMAP(obj,
               dims = dims,
               reduction = "harmony",
               reduction.name = "umap_harmony")

obj <- FindNeighbors(obj,
                     dims = dims,
                     reduction = "harmony",
                     annoy.metric = "cosine")

obj <- FindClusters(obj,
                    resolution = resolutions,
                    algorithm = 4,
                    method = "igraph")
saveRDS(obj,
        file.path(out_dir, "seurat.integrated_harmony_20_20.rds"))

# Remove cluster 12 (resolution 0.2) and shorten condition labels
obj <- SetIdent(obj, value = "RNA_snn_res.0.2")
obj <- subset(obj, idents = 12, invert = TRUE)
obj$Condition <- ifelse(
  obj$Condition == "Normal",
  "ND",
  ifelse(obj$Condition == "auto_antibody", "AAB", obj$Condition)
)

saveRDS(obj,
        file.path(
          out_dir,
          "seurat.integrated_harmony_20_20_clean.rds"
        ))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_harmony.txt"))