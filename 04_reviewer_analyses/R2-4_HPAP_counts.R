################################################################################
### R2-4: HPAP donor and study contributions to beta-cell clusters
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

################################################################################
### 1. Load beta-cell metadata
################################################################################

data.integrated <- readRDS('beta.integrated.rds')
df <- data.integrated@meta.data

################################################################################
### 2. HPAP donor counts by cluster
################################################################################

hpap_by_cluster <- df %>%
  filter(study == 'HPAP') %>%
  count(orig.ident, monocle_clusters, name = 'n_cells') %>%
  arrange(orig.ident, monocle_clusters)
write.table(
  hpap_by_cluster,
  'HPAP_donor_by_cluster_counts.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

hpap_wide <- hpap_by_cluster %>%
  pivot_wider(
    names_from = monocle_clusters,
    values_from = n_cells,
    values_fill = 0,
    names_prefix = 'cluster_'
  )
write.table(
  hpap_wide,
  'HPAP_donor_by_cluster_counts_wide.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

################################################################################
### 3. Study contributions by cluster
################################################################################

study_by_cluster <- df %>%
  count(study, monocle_clusters, name = 'n_cells') %>%
  arrange(study, monocle_clusters)
write.table(
  study_by_cluster,
  'study_by_cluster_counts.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
