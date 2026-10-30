################################################################################
### R1-6: HPAP104 exclusion and root-selection sensitivity
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(slingshot)
  library(SingleCellExperiment)
})
set.seed(123)

################################################################################
### 1. Load original beta-cell trajectory
################################################################################

data.integrated <- readRDS('beta.integrated.rds')
cds <- readRDS('cds_all.rds')
df <- data.integrated@meta.data

sim <- SingleCellExperiment(assays = list(counts = GetAssayData(
  data.integrated, assay = 'SCT', slot = 'counts'
)))
reducedDims(sim)$UMAP <- reducedDims(cds)$UMAP
sim$monocle_clusters <- as.character(df$monocle_clusters)

sim <- slingshot(
  sim,
  clusterLabels = 'monocle_clusters',
  reducedDim = 'UMAP',
  start.clus = '3'
)

################################################################################
### 2. Exclude HPAP104 and refit trajectories
################################################################################

keep_cells <- colnames(sim)[df$orig.ident != 'HPAP104']
sim_no_HPAP104 <- sim[, keep_cells]
sim_no_HPAP104 <- slingshot(
  sim_no_HPAP104,
  clusterLabels = 'monocle_clusters',
  reducedDim = 'UMAP',
  start.clus = '3'
)

lineages_HPAP104 <- slingLineages(sim_no_HPAP104)
lineages_HPAP104 <- data.frame(
  lineage = names(lineages_HPAP104),
  topology = sapply(lineages_HPAP104, paste, collapse = ' -> ')
)
write.table(
  lineages_HPAP104,
  'R1-6_HPAP104_exclusion_lineages.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

original_pt <- slingPseudotime(sim)
new_pt <- slingPseudotime(sim_no_HPAP104)
shared <- intersect(rownames(original_pt), rownames(new_pt))

## Compare matching lineage columns, as in the original revision.
pt_cor <- sapply(seq_len(min(ncol(original_pt), ncol(new_pt))), function(i) {
  cor(original_pt[shared, i], new_pt[shared, i], method = 'spearman', use = 'complete.obs')
})
write.table(
  data.frame(
    lineage = paste0('Lineage', seq_along(pt_cor)),
    spearman_rho = pt_cor
  ),
  'R1-6_HPAP104_pseudotime_correlations.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

################################################################################
### 3. Slingshot without a specified root
################################################################################

sim_unrooted <- slingshot(sim, clusterLabels = 'monocle_clusters', reducedDim = 'UMAP')
lineages_unrooted <- slingLineages(sim_unrooted)
write.table(
  data.frame(
    lineage = names(lineages_unrooted),
    topology = sapply(lineages_unrooted, paste, collapse = ' -> ')
  ),
  'R1-6_unrooted_lineages.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
