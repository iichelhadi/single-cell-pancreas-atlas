#!/usr/bin/env Rscript
# Apply the trained scPred model to the integrated atlas.
# The scPred version must match the one used for training.


suppressPackageStartupMessages({
  library(Seurat)
  library(scPred)
})

# ---- CONFIG -----------------------------------------------------------------
in_rds    <- "output/seurat.integrated_harmony_20_20_clean.rds"
model_rds <- "output/scPred_mda.rds"
out_dir   <- "output"
# -----------------------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
x <- readRDS(in_rds)
model <- readRDS(model_rds)
stopifnot("RNA_snn_res.0.2" %in% colnames(x@meta.data))
Idents(x) <- "RNA_snn_res.0.2"

# NOTE: the original run used a manual patch to scPred::project_query for RNA assay access.
x <- scPredict(x, model, max.iter.harmony = 50, threshold = 0.9)

saveRDS(x,
        file.path(
          out_dir,
          "scRNA_integrated_harmony_classified_20_20.rds"
        ))
