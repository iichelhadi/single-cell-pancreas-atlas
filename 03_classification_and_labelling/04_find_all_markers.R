suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(ggsci)
  library(scCustomize)
  library(MAST)
})

set.seed(123)

merged_seurat <- readRDS('output/seurat.integrated_harmony_20_20_clean.rds')

merged_seurat <- SetIdent(merged_seurat, value = 'RNA_snn_res.0.2')

DEGs <- FindAllMarkers(object = merged_seurat, verbose = T, test.use = 'MAST',
                       latent.vars = c('orig.ident','Platform'),
                       max.cells.per.ident = 10000,
                       only.pos = F)

DEGs <- DEGs %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write_tsv(DEGs, 'DEGs_allmarkers_assigned_label.tsv')
