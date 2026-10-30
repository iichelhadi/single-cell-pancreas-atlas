################################################################################
# Project: scRNA beta cell atlas — Figure 5H pySCENIC inputs
# File: Figure5H_prepare_pySCENIC_inputs.R
# Author: Elhadi Iich
#
# Description:
#   Writes the expression matrix and TF list used by the Figure 5H pySCENIC
#   consensus workflow. Setup lines taken from the original Lineage 2 C3
#   SCENIC script (SCENIC2.R): C3 cells, SCT data, genes detected in >= 1%
#   of cells.
#
# Outputs:
#   - exprMat_filtered.rds
#   - 1.1_exprMatrix_filtered_t.txt   (original single GRNBoost2 run)
#   - 1.1_inputTFs.txt                (TF list, all GRNBoost2 runs)
#   - expression_matrix.csv           (cells x genes, cell IDs in column 1;
#                                      run_grnboost2_multi_v3.py, pyscenic ctx/aucell)
#
# Next steps:
#   run_grnboost2_multi_v3.py --n-runs 50 -> build_consensus_pyscenic_v2.R ->
#   run_pyscenic_ctx_aucell_v2.pbs -> compare_thresholds_post_ctx.py ->
#   Figure5H_define_C3_lineage_groups.R -> Figure5H_SCENIC_consensus_analysis.R ->
#   Figure5H_SCENIC_figures.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(SCENIC)
  library(RcisTarget)
  library(Seurat)
})

set.seed(123)

nCores = 20
org = 'hgnc'
datasetTitle = 'Beta_lin2_C3'

system('mkdir int')

lin2_obj <- readRDS('./revision/Figure4_S6/lin2.beta.reintegrated.rds')
dim(lin2_obj)
# [1] 15360 13853

lin2_obj <- subset(lin2_obj, subset = monocle_clusters2 == '3')
dim(lin2_obj)
# [1] 15360  2762

### Expression matrix (SCT data)
exprMat <- GetAssayData(lin2_obj, assay = 'SCT', slot = 'data')
exprMat <- as.matrix(exprMat)

### SCENIC options (needed by exportsForArboreto for the TF list)
hg19_dbs <-
  list('500bp' = 'hg19-500bp-upstream-10species.mc9nr.feather',
       '5kb' = 'hg19-tss-centered-5kb-10species.mc9nr.feather',
       '10kb' = 'hg19-tss-centered-10kb-10species.mc9nr.feather')

db_path = './Indices/SCENIC/cisTarget_databases/hg19/'   # path to hg19 cisTarget databases

data(list = "motifAnnotations_hgnc", package = "RcisTarget")

scenicOptions <- initializeScenic(
  org = org,
  dbDir = db_path,
  dbs = hg19_dbs,
  datasetTitle = datasetTitle,
  nCores = nCores
)

scenicOptions@settings$seed <- 123

### Filter genes expressed in at least 1% of cells
min_cells <- 0.01 * ncol(exprMat)
genesKept <- rowSums(exprMat > 0) >= round(min_cells)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

saveRDS(exprMat_filtered, 'exprMat_filtered.rds')

### Export for GRNBoost2
exportsForArboreto(exprMat_filtered, scenicOptions, dir = "./")   # 1.1_exprMatrix_filtered_t.txt, 1.1_inputTFs.txt

### pySCENIC input: cells as rows, cell IDs as first column
write.csv(t(exprMat_filtered), 'expression_matrix.csv', row.names = TRUE)

