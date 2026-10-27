#!/usr/bin/env Rscript
# Apply manual cluster-to-cell-type labels (resolution 0.2).

suppressPackageStartupMessages(library(Seurat))

# ---- CONFIG -----------------------------------------------------------------
in_rds  <- "output/scRNA_integrated_harmony_classified_20_20.rds"
out_rds <- "output/scRNA_integrated_harmony_classified"
out_tsv <- "output/scPred_meta_data_res0.2.tsv"

# Labels in cluster order (clusters 1-11 at RNA_snn_res.0.2)
labels <- c(
  "Alpha",
  "Ductal",
  "Acinar",
  "Beta",
  "Pancreatic Stellate",
  "Endothelial",
  "Delta",
  "Pancreatic Polypeptide",
  "MHC Class II Cell",
  "Bihormonal",
  "Cycling"
)
# -----------------------------------------------------------------------------

x <- readRDS(in_rds)
Idents(x) <- "RNA_snn_res.0.2"
if (nlevels(x) != length(labels))
  stop("Expected ", length(labels), " clusters, found ", nlevels(x))
names(labels) <- levels(x)
x$assigned_label <- unname(labels[as.character(Idents(x))])

# Manual refinement: GHRL-high cells in the pancreatic polypeptide cluster are epsilon cells
epsilon <- Cells(subset(x, subset = GHRL > 4 &
                          assigned_label == "Pancreatic Polypeptide"))
x$assigned_label[Cells(x) %in% epsilon] <- "Epsilon"
Idents(x) <- "assigned_label"

x$Condition[x$Condition == "auto_antibody"] <- "AAB"
x$Condition[x$Condition == "Normal"] <- "ND"

write.table(
  x@meta.data,
  out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(x, out_rds)
