################################################################################
### R1-3: Reintegrate alpha cells with atlas cluster 11
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(monocle3)
  library(SingleCellExperiment)
})
set.seed(123)

################################################################################
### 1. Alpha + cluster 11 integration
################################################################################

atlas <- readRDS('scRNA_integrated_harmony_classified.rds')

cells_use <- colnames(atlas)[atlas$Condition %in% c('ND', 'T2D') &
                               (atlas$assigned_label == 'Alpha' |
                                  atlas$RNA_snn_res.0.2 == '11')]
data <- subset(atlas, cells = cells_use)
rm(atlas)

data.list <- SplitObject(data, split.by = 'study')
data.list <- data.list[sapply(data.list, ncol) > 100]
rm(data)

data.list <- lapply(
  data.list,
  SCTransform,
  vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),
  return.only.var.genes = FALSE
)
features <- SelectIntegrationFeatures(data.list)
data.list <- PrepSCTIntegration(data.list, anchor.features = features)
anchors <- FindIntegrationAnchors(
  object.list = data.list,
  normalization.method = 'SCT',
  k.filter = NA,
  anchor.features = features,
  dims = 1:30
)

data.integrated <- IntegrateData(anchors, dims = 1:30, normalization.method = 'SCT')
data.integrated <- RunPCA(data.integrated, verbose = FALSE)
saveRDS(data.integrated, 'alpha_C11.integrated.rds')

################################################################################
### 2. Monocle3 clustering
################################################################################

expr <- GetAssayData(data.integrated, assay = 'SCT', slot = 'counts')
gene_metadata <- data.frame(gene_short_name = rownames(expr),
                            row.names = rownames(expr))
cds <- new_cell_data_set(expr,
                         cell_metadata = data.integrated@meta.data,
                         gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim = 10)
reducedDims(cds)$PCA <- Embeddings(data.integrated, 'pca')[, 1:10]
cds <- reduce_dimension(cds, reduction_method = 'UMAP')
cds <- cluster_cells(cds, cluster_method = 'leiden', resolution = 0.00001)

data.integrated$monocle_clusters <- as.character(clusters(cds))

################################################################################
### 3. Reanalysis of atlas cluster 11
################################################################################

df <- data.integrated@meta.data
cluster11 <- df %>%
  group_by(monocle_clusters) %>%
  summarise(
    n_total = n(),
    n_atlas_C11 = sum(RNA_snn_res.0.2 == '11'),
    percent_atlas_C11 = 100 * n_atlas_C11 / n_total,
    .groups = 'drop'
  )

write.table(
  cluster11,
  'R1-3_cluster11_reintegration.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
saveRDS(data.integrated, 'alpha_C11.clustered.rds')
saveRDS(cds, 'alpha_C11.cds.rds')
