################################################################################
### R2-6: Lineage 2 C5 versus atlas acinar cells
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(presto)
})
set.seed(123)

################################################################################
### 1. Identify lineage 2 C5 and atlas acinar cells
################################################################################

## Important: C5 is taken from the lineage-2 object, not global beta cluster 5.
lin2 <- readRDS('lin2.beta.reintegrated.rds')
merged_seurat <- readRDS('scRNA_integrated_harmony_classified.rds')

c5_cells <- colnames(lin2)[lin2$monocle_clusters2 == '5']
acinar_cells <- colnames(merged_seurat)[grepl('Acinar', merged_seurat$assigned_label, ignore.case = TRUE)]

## Retain cell identifiers present in the original atlas.

comparison_obj <- subset(merged_seurat, cells = c(c5_cells, acinar_cells))
comparison_obj$comparison_group <- ifelse(colnames(comparison_obj) %in% c5_cells, 'L2_C5', 'Acinar')
DefaultAssay(comparison_obj) <- 'RNA'

################################################################################
### 2. Marker expression comparison
################################################################################

markers <- c(
  'CPA1',
  'PRSS1',
  'PRSS2',
  'REG1A',
  'CTRB2',
  'INS',
  'IAPP',
  'PCSK1',
  'CPE',
  'PDX1',
  'NKX6-1'
)
markers <- markers[markers %in% rownames(comparison_obj)]

p <- DotPlot(comparison_obj, features = markers, group.by = 'comparison_group') + RotatedAxis()
ggsave('C5_vs_acinar_DotPlot.pdf',
       p,
       width = 8,
       height = 4.5)

################################################################################
### 3. Differential marker comparison
################################################################################

Idents(comparison_obj) <- 'comparison_group'
c5_vs_acinar <- presto::wilcoxauc(
  comparison_obj,
  group_by = 'comparison_group',
  assay = 'data',
  seurat_assay = 'RNA'
)
marker_results <- c5_vs_acinar %>%
  filter(group == 'L2_C5', feature %in% markers) %>%
  arrange(match(feature, markers))

write.table(
  marker_results,
  'C5_vs_acinar_marker_comparison.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
