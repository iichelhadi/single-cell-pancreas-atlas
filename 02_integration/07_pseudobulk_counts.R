################################################################################
# Pseudobulk count matrices for the cluster comparisons tested with DESeq2
# in Figures 3L, 4L, 5A-C and 6I.
#
# Standardises the original per-comparison scripts (C3_vs_C4_DEG_analysis.R
# and the other Lineage 1 C*_vs_C* scripts). Only the pseudobulk construction
# is kept; DESeq2 is run in the figure scripts.
#
# Pseudobulk unit: study x Condition x cluster (monocle_clusters2), i.e. the
# raw RNA counts of all cells from one study, condition and cluster are summed.
# The number of cells in each pseudobulk sample is stored in coldata$counts.
#
# Run after 02-04 from the directory holding the lineage objects; the
# Figure* output folders are created there if absent.
#
# Inputs:
#   Lin1.beta.reintegrated.rds, lin2.beta.reintegrated.rds,
#   lin3.beta.reintegrated.rds
#
# Outputs (cts + coldata TSV per comparison), named as read by the figure
# scripts and written to the folder of the figure script that reads them
# (Lineage 1 for Fig 3L and Lineage 2 for Fig 5A both use beta_C3_v_C4):
#   Figure3_S5/beta_C3_v_C4_pseudobulk_*    Lineage 1 C3 v C4   Figure3_FigureS5_clean.R
#   Figure4_S6/beta_C2_v_C5_pseudobulk_*    Lineage 2 C2 v C5   Figure4_FigureS6_clean.R
#   Figure5_S8/beta_C3_v_C4_pseudobulk_*    Lineage 2 C3 v C4   Figure5_FigureS8_clean.R
#   Figure5_S8/beta_C2_v_C3_pseudobulk_*    Lineage 2 C2 v C3   Figure5_FigureS8_clean.R
#   Figure5_S8/beta_C3_v_C5_pseudobulk_*    Lineage 2 C3 v C5   Figure5_FigureS8_clean.R
#   Figure6_S10/beta_lin3_T2D_pseudobulk_*  Lineage 3 C2 v C3   Figure6_FigureS10_clean.R
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

set.seed(123)

################################################################################
### Comparisons (one row per pseudobulk table)
################################################################################

comparisons <- data.frame(
  object = c('Lin1.beta.reintegrated.rds',
             'lin2.beta.reintegrated.rds',
             'lin2.beta.reintegrated.rds',
             'lin2.beta.reintegrated.rds',
             'lin2.beta.reintegrated.rds',
             'lin3.beta.reintegrated.rds'),
  cluster_a = c('3', '2', '3', '2', '3', '2'),
  cluster_b = c('4', '5', '4', '3', '5', '3'),
  out_name  = c('beta_C3_v_C4',
                'beta_C2_v_C5',
                'beta_C3_v_C4',
                'beta_C2_v_C3',
                'beta_C3_v_C5',
                'beta_lin3_T2D'),
  out_dir   = c('Figure3_S5', 'Figure4_S6', 'Figure5_S8',
                'Figure5_S8', 'Figure5_S8', 'Figure6_S10')
)

current_object <- ''

for (i in seq_len(nrow(comparisons))) {

  comp <- comparisons[i, ]
  print(paste('########################', comp$out_dir, comp$out_name, '########################'))
  dir.create(comp$out_dir, showWarnings = FALSE)

  if (comp$object != current_object) {
    lineage_obj <- readRDS(comp$object)
    current_object <- comp$object
    print(table(lineage_obj$monocle_clusters2))
  }

  data.integrated <-
    subset(lineage_obj,
           subset = monocle_clusters2 %in% c(comp$cluster_a, comp$cluster_b))
  print(dim(data.integrated))

  ##############################################################################
  ### Generate pseudobulk
  ##############################################################################

  data.integrated@meta.data$pseudobulk <-
    paste0(data.integrated@meta.data$study, '_',
           data.integrated@meta.data$Condition, '_',
           data.integrated@meta.data$monocle_clusters2)

  meta_counts <-
    as.data.frame(t(t(table(
      data.integrated@meta.data$pseudobulk
    ))))

  pseudobulk <- AggregateExpression(data.integrated,
                                    return.seurat = F,
                                    slot = 'counts',
                                    assay = 'RNA',
                                    group.by = c('study',
                                                 'Condition',
                                                 'monocle_clusters2'))

  cts <- as.data.frame(pseudobulk$RNA)
  print(dim(cts))

  coldata <-
    data.frame(
      study = str_split_fixed(colnames(cts), '_', n = 3)[, 1],
      Condition = str_split_fixed(colnames(cts), '_', n = 3)[, 2],
      monocle_clusters2 = str_split_fixed(colnames(cts), '_', n = 3)[, 3])

  rownames(coldata) <- colnames(cts)

  #! check that cell counts are attached to the matching pseudobulk sample
  stopifnot(identical(as.character(meta_counts$Var1), colnames(cts)))
  coldata <- cbind(coldata, counts = meta_counts$Freq)

  #! Figure 6I reads a 'monocle_clusters' column with values 'C2'/'C3';
  #! the original Lineage 3 script was not available, verify against the
  #! existing beta_lin3_T2D_pseudobulk_coldata.tsv
  if (comp$out_name == 'beta_lin3_T2D') {
    coldata$monocle_clusters <- paste0('C', coldata$monocle_clusters2)
  }

  write.table(
    cts,
    sep = '\t',
    col.names = NA,
    row.names = T,
    file.path(comp$out_dir, paste0(comp$out_name, '_pseudobulk_cts.tsv'))
  )

  write.table(
    coldata,
    sep = '\t',
    col.names = NA,
    row.names = T,
    file.path(comp$out_dir, paste0(comp$out_name, '_pseudobulk_coldata.tsv'))
  )

  rm(data.integrated, pseudobulk, cts, coldata, meta_counts)
  gc()
}

################################################################################

sessionInfo()
