#!/usr/bin/env Rscript

# SoupX ambient RNA correction for 10x Chromium data
# Usage: Rscript run_soupx.R samples.tsv

library(SoupX)
library(DropletUtils)

# Load sample information
samples <- read.delim("samples.tsv")

for (i in seq_len(nrow(samples))) {
  sample <- samples$sample_id[i]
  input <- samples$cellranger_outs[i]
  output <- samples$output_dir[i]
  
  message("Processing: ", sample)
  
  # Ambient RNA correction
  sc <- load10X(input)
  sc <- autoEstCont(sc)
  counts <- adjustCounts(sc)
  
  # Save corrected counts
  write10xCounts(output, counts)
  
  # Preserve original Cell Ranger gene annotations
  features <- file.path(input, "filtered_feature_bc_matrix", "features.tsv.gz")
  
  stopifnot(file.exists(features))
  stopifnot(nrow(counts) == length(readLines(gzfile(features))))
  
  file.copy(features, file.path(output, "features.tsv.gz"), overwrite = TRUE)
  
  # Remove alternative gene annotation files
  unlink(file.path(output, c(
    "genes.tsv", "genes.tsv.gz", "features.tsv"
  )))
  
  message("Completed: ", sample)
}