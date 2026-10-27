#!/usr/bin/env Rscript
# Train an scPred model on previously annotated singlet cells.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(scPred)
})

# ---- CONFIG -----------------------------------------------------------------
in_rds       <- "output/scRNA_consolidation_part1.rds"  # needs inferred_cell_type, doublet_final
epsilon_file <- "GSE150724_epsilon_cells.txt"       # one cell barcode per line
out_dir      <- "output"
train_studies <- c("EMTAB5061", "GSE81547", "GSE81608", "GSE83139", "GSE150724")
drop_labels   <- c("Unclassified Cell", "unknown", "Islet", "Co-Expression Cell")
# -----------------------------------------------------------------------------

set.seed(123)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(in_rds)
stopifnot(all(
  c(
    "inferred_cell_type",
    "study",
    "doublet_final",
    "Condition",
    "orig.ident",
    "Platform"
  ) %in% colnames(x@meta.data)
))

# Label epsilon cells, then restrict to training studies and singlets
Idents(x) <- "inferred_cell_type"
Idents(x, cells = intersect(readLines(epsilon_file, warn = FALSE), Cells(x))) <- "Epsilon Cell"
x$inferred_cell_type <- as.character(Idents(x))
x <- subset(x, subset = study %in% train_studies &
              doublet_final == "Singlet")

x$inferred_cell_type[x$inferred_cell_type == "Mesenchymal"] <- "Pancreatic Stellate Cell"
x$inferred_cell_type[is.na(x$inferred_cell_type)] <- "unknown"
x <- subset(x,
            subset = !(inferred_cell_type %in% drop_labels) &
              Condition != "Unlabeled")

x <- NormalizeData(x)
x <- FindVariableFeatures(
  x,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

x <- ScaleData(
  x,
  features = VariableFeatures(x),
  vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")
)
x <- RunPCA(x, features = VariableFeatures(x))
x <- RunHarmony(
  x,
  group.by.vars = c("orig.ident", "Platform", "Condition"),
  reduction.use = "pca",
  max.iter.harmony = 50,
  lambda = NULL,
  dims.use = 1:50
)

# Top 20 Harmony dimensions ranked by SD (as in the original)
x <- RunUMAP(x, reduction = "harmony", dims = 1:20)

x <- getFeatureSpace(x, pvar = "inferred_cell_type")
x <- trainModel(
  x,
  seed = 123,
  allowParallel = TRUE,
  model = "mda",
  number = 50,
  returnData = TRUE
)

saveRDS(x, file.path(out_dir, "scPred_harmony.rds"))
saveRDS(get_scpred(x), file.path(out_dir, "scPred_mda.rds"))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))
