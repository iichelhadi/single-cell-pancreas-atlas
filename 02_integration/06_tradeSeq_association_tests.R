################################################################################
# tradeSeq fits and start-vs-end association tests for each lineage branch
#
# Consolidates the five original per-branch scripts:
#   20250303_lineage1_1_association_test.R, 20250303_lineage1_2_association_test.R,
#   lineage2_1_association_test.R, lineage2_2_association_test.R,
#   lineage3_association_test.R
# Only the steps that produce saved outputs are kept. Heatmaps, GSEA, volcano
# and gene plots are in the Figure 3, 4 and 6 scripts.
#
# Run after 01-04, from the same working directory (independent of 05).
#
# Inputs:
#   slingshot_sim_lin1.rds, slingshot_sim_lin2.rds, slingshot_sim_lin3.rds
#   Lineage1_sling_perc_expr.tsv, Lineage2_sling_perc_expr.tsv,
#   Lineage3_sling_perc_expr.tsv
#
# Outputs:
#   sce_lin1_1.rds, sce_lin1_2.rds, sce_lin2_1.rds, sce_lin2_2.rds, sce_lin3.rds
#   lin1_1_pseudotime_start_end_association_filt.tsv (and lin1_2, lin2_1,
#   lin2_2, lin3)
################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
})

set.seed(123)

################################################################################
### Branch settings (one row per branch)
################################################################################

branches <- data.frame(
  branch    = c('lin1_1', 'lin1_2', 'lin2_1', 'lin2_2', 'lin3'),
  sim_file  = c('slingshot_sim_lin1.rds', 'slingshot_sim_lin1.rds',
                'slingshot_sim_lin2.rds', 'slingshot_sim_lin2.rds',
                'slingshot_sim_lin3.rds'),
  perc_file = c('Lineage1_sling_perc_expr.tsv', 'Lineage1_sling_perc_expr.tsv',
                'Lineage2_sling_perc_expr.tsv', 'Lineage2_sling_perc_expr.tsv',
                'Lineage3_sling_perc_expr.tsv'),
  lineage   = c('Lineage1', 'Lineage2', 'Lineage1', 'Lineage2', 'Lineage1')
)

## Values recorded in the original scripts, for checking the rerun:
## branch   cells   genes in filtered start-vs-end table
## lin1_1   12283   2006
## lin1_2   13372   6185
## lin2_1   12312   3810
## lin2_2   10981   2657
## lin3      8882   2991

for (i in seq_len(nrow(branches))) {

  b <- branches[i, ]
  print(paste('########################', b$branch, '########################'))

  sim <- readRDS(b$sim_file)

  perc_expr <- read_tsv(b$perc_file)
  colnames(perc_expr) <- c('feature_id', 'perc_expr')

  ##############################################################################
  ### Cells on this branch
  ##############################################################################

  pseudotime  <- slingPseudotime(SlingshotDataSet(sim))[, b$lineage]
  cellWeights <- slingCurveWeights(SlingshotDataSet(sim))[, b$lineage]

  validCells <- !is.na(pseudotime)

  filter_counts <- sim@assays@data@listData$counts[, validCells]
  print(dim(filter_counts))

  U <- model.matrix( ~ 0 + colData(sim)$study[validCells])
  print(dim(U))

  ##############################################################################
  ### fitGAM
  ##############################################################################

  sce <-
    tradeSeq::fitGAM(
      counts = filter_counts,
      pseudotime = pseudotime[validCells],
      cellWeights = cellWeights[validCells],
      U = U,
      nknots = 9,
      verbose = T
    )

  colData(sce) <-
    cbind(barcode = rownames(colData(sce)), colData(sce))
  rownames(colData(sce)) <- colData(sce)$barcode

  ##############################################################################
  ### associationTest
  ##############################################################################

  rowData(sce)$assocRes <-
    associationTest(sce, lineages = FALSE, l2fc = log2(2))

  rowData(sce)$assocRes$fdr <-
    p.adjust(rowData(sce)$assocRes$pvalue, method = "fdr")

  rowData(sce)$assocRes <-
    rowData(sce)$assocRes[order(rowData(sce)$assocRes$pvalue), ]

  rowData(sce)$assocRes$feature_id <-
    rownames(rowData(sce)$assocRes)

  saveRDS(sce, paste0('sce_', b$branch, '.rds'))

  ##############################################################################
  ### startVsEndTest: FDR <= 0.05 and expressed in >= 5% of lineage cells
  ##############################################################################

  pseudotime_start_end_association <- startVsEndTest(sce)
  pseudotime_start_end_association$padj <-
    p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
  pseudotime_start_end_association$feature_id <-
    rownames(pseudotime_start_end_association)

  pseudotime_start_end_association_filt <-
    pseudotime_start_end_association[pseudotime_start_end_association$padj <= 0.05, ]
  print(dim(pseudotime_start_end_association_filt))

  pseudotime_start_end_association_filt <-
    left_join(pseudotime_start_end_association_filt,
              perc_expr, by = 'feature_id')

  pseudotime_start_end_association_filt <-
    pseudotime_start_end_association_filt[pseudotime_start_end_association_filt$perc_expr >= 5, ]
  print(dim(pseudotime_start_end_association_filt))

  write_tsv(
    pseudotime_start_end_association_filt,
    paste0(b$branch, '_pseudotime_start_end_association_filt.tsv')
  )

  rm(sim, sce, filter_counts)
  gc()
}

################################################################################

sessionInfo()
