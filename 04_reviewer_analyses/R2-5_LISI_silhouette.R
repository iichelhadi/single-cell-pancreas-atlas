################################################################################
### R2-5: Integration quality (LISI and silhouette)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(lisi)
  library(cluster)
})
set.seed(123)

################################################################################
### 1. Global LISI
################################################################################

merged_seurat <- readRDS('scRNA_integrated_harmony_classified.rds')
harmony_emb <- Embeddings(merged_seurat, 'harmony')
meta_lisi <- merged_seurat@meta.data %>%
  select(study, Platform, assigned_label)

lisi_res <- compute_lisi(harmony_emb,
                         meta_lisi,
                         label_colnames = c('study', 'Platform', 'assigned_label'))
global_summary <- data.frame(
  metric = c('Study', 'Platform', 'Cell_type'),
  median_LISI = c(
    median(lisi_res$study),
    median(lisi_res$Platform),
    median(lisi_res$assigned_label)
  ),
  mean_LISI = c(
    mean(lisi_res$study),
    mean(lisi_res$Platform),
    mean(lisi_res$assigned_label)
  )
)
write.table(
  global_summary,
  'global_LISI_summary.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

################################################################################
### 2. LISI for each major cell type
################################################################################

major_labels <- names(which(table(merged_seurat$assigned_label) >= 1000))

study_lisi_by_celltype <- map_dfr(major_labels, function(ct) {
  cells_ct <- colnames(merged_seurat)[merged_seurat$assigned_label == ct]
  meta_ct <- merged_seurat@meta.data[cells_ct, , drop = FALSE]
  eligible_studies <- names(which(table(meta_ct$study) >= 50))
  if (length(eligible_studies) < 2)
    return(NULL)
  
  cells_balanced <- unlist(lapply(eligible_studies, function(s) {
    cells_s <- rownames(meta_ct)[meta_ct$study == s]
    sample(cells_s, min(length(cells_s), 500))
  }))
  
  emb <- harmony_emb[cells_balanced, , drop = FALSE]
  res <- compute_lisi(emb, merged_seurat@meta.data[cells_balanced, 'study', drop = FALSE], label_colnames = 'study')
  data.frame(
    cell_type = ct,
    n_cells = length(cells_balanced),
    n_studies = length(eligible_studies),
    median_study_LISI = median(res$study),
    mean_study_LISI = mean(res$study)
  )
})
write.table(
  study_lisi_by_celltype,
  'study_LISI_by_celltype.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

################################################################################
### 3. Study silhouette for each major cell type
################################################################################

silhouette_by_celltype <- map_dfr(major_labels, function(ct) {
  cells_ct <- colnames(merged_seurat)[merged_seurat$assigned_label == ct]
  meta_ct <- merged_seurat@meta.data[cells_ct, , drop = FALSE]
  eligible_studies <- names(which(table(meta_ct$study) >= 50))
  if (length(eligible_studies) < 2)
    return(NULL)
  
  cells_balanced <- unlist(lapply(eligible_studies, function(s) {
    cells_s <- rownames(meta_ct)[meta_ct$study == s]
    sample(cells_s, min(length(cells_s), 500))
  }))
  
  emb <- harmony_emb[cells_balanced, 1:20, drop = FALSE]
  study_factor <- factor(merged_seurat$study[cells_balanced])
  sil <- silhouette(as.integer(study_factor), dist(emb))
  
  data.frame(
    cell_type = ct,
    n_cells = length(cells_balanced),
    n_studies = length(eligible_studies),
    mean_study_silhouette = mean(sil[, 'sil_width']),
    median_study_silhouette = median(sil[, 'sil_width'])
  )
})
write.table(
  silhouette_by_celltype,
  'study_silhouette_by_celltype.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
