###############################################################################
# Project: Pancreas scRNA atlas
# Analysis: B2 beta-cell trajectory sensitivity including atlas cluster 10
#
# Description:
#   Repeats beta-cell integration, clustering and Slingshot trajectory inference
#   after inclusion of atlas cluster 10 (bihormonal population).
#
# Generates:
#   Supplementary Fig S3G:
#     Reintegrated beta-cell UMAP with Slingshot trajectories
#
#   Supplementary Fig S3H:
#     Same trajectory with atlas cluster 10 cells highlighted
#
#   Source-data outputs:
#     FigS3G_S3H_cell_coordinates.tsv
#     FigS3G_S3H_trajectory_curves.tsv
#     FigS3G_S3H_cluster10_distribution.tsv
#
# Input:
#   scRNA_integrated_harmony_classified.rds
#
# Output:
#   beta_c10.integrated.rds
###############################################################################


###############################################################################
# Packages
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(slingshot)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(tidyverse)
  library(patchwork)
  library(reticulate)
})

set.seed(123)

cores <- 1
resolution <- 0.00004


###############################################################################
# Input / output
###############################################################################

input_file <- "scRNA_integrated_harmony_classified.rds"

integrated_file <- "beta_c10.integrated.rds"

outdir <- "Fig2_S3"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. Subset beta cells together with atlas cluster 10
###############################################################################

atlas <- readRDS(input_file)

beta_c10 <- subset(
  atlas,
  cells = colnames(atlas)[
    atlas$Condition %in% c("ND", "T2D") &
      (
        atlas$assigned_label == "Beta" |
          as.character(atlas$RNA_snn_res.0.2) == "10"
      )
  ]
)


###############################################################################
# 2. SCTransform integration
###############################################################################

data.list <- SplitObject(
  beta_c10,
  split.by = "study"
)

rm(beta_c10, atlas)
gc()

# Match the primary beta-cell analysis:
# remove datasets contributing <=100 cells
data.list <- data.list[
  vapply(
    data.list,
    ncol,
    integer(1)
  ) > 100
]


length(data.list) # 12


data.list <- lapply(
  data.list,
  SCTransform,
  vars.to.regress = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  ),
  return.only.var.genes = FALSE,
  verbose = FALSE
)

features <- SelectIntegrationFeatures(
  object.list = data.list
)

data.list <- PrepSCTIntegration(
  object.list = data.list,
  anchor.features = features
)

anchors <- FindIntegrationAnchors(
  object.list = data.list,
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:30,
  k.filter = NA
)

data.integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30
)

rm(data.list, anchors)
gc()

data.integrated <- RunPCA(
  data.integrated,
  verbose = FALSE
)

saveRDS(
  data.integrated,
  integrated_file
)


###############################################################################
# 3. Monocle3 UMAP and clustering
###############################################################################

expr <- GetAssayData(
  data.integrated,
  assay = "SCT",
  slot = "counts"
)

gene_metadata <- data.frame(
  gene_short_name = rownames(expr),
  row.names = rownames(expr)
)

cell_metadata <- data.integrated@meta.data

cds <- new_cell_data_set(
  expression_data = expr,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

rm(expr)
gc()


# Use the same PCA representation as the integrated Seurat object
reducedDims(cds)$PCA <-
  Embeddings(data.integrated, "pca")[, 1:20, drop = FALSE]


# Treat all cells as one partition
partition <- factor(
  rep(1, ncol(cds))
)

names(partition) <- colnames(cds)

cds@clusters$UMAP$partitions <- partition


set.seed(123)

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  cores = cores,
  umap.fast_sgd = FALSE
)

cds <- cluster_cells(
  cds,
  cluster_method = "leiden",
  resolution = resolution,
  verbose = TRUE
)


colData(cds)$monocle_clusters <-
  as.character(
    cds@clusters$UMAP$clusters
  )

data.integrated$monocle_clusters <-
  colData(cds)$monocle_clusters


###############################################################################
# 4. Slingshot trajectory inference
###############################################################################

sce <- SingleCellExperiment(
  assays = list(
    counts = GetAssayData(
      data.integrated,
      assay = "SCT",
      slot = "counts"
    )
  )
)

reducedDims(sce)$UMAP <-
  reducedDims(cds)$UMAP

colData(sce)$monocle_clusters <-
  colData(cds)$monocle_clusters


# In the expanded clustering, cluster 4 corresponds to the original
# primary-analysis root cluster 3.
sce <- slingshot(
  sce,
  clusterLabels = "monocle_clusters",
  reducedDim = "UMAP",
  start.clus = "4"
)


print(
  slingLineages(
    SlingshotDataSet(sce)
  )
)


###############################################################################
# 5. Extract trajectory curves
###############################################################################

curves <- slingCurves(
  SlingshotDataSet(sce)
)

trajectory_data <- lapply(
  names(curves),
  function(lineage_name) {
    
    x <- as.data.frame(
      curves[[lineage_name]]$s
    )
    
    colnames(x) <- c(
      "UMAP_1",
      "UMAP_2"
    )
    
    x$lineage <- lineage_name
    x$curve_order <- seq_len(nrow(x))
    
    x
  }
) %>%
  bind_rows()


###############################################################################
# 6. Cell-level plotting data
###############################################################################

cell_umap <- as.data.frame(
  reducedDims(cds)$UMAP
)

colnames(cell_umap) <- c(
  "UMAP_1",
  "UMAP_2"
)

cell_umap$cell <- rownames(cell_umap)

cell_umap$cluster <-
  as.character(
    colData(cds)$monocle_clusters
  )

cell_umap$atlas_cluster10 <-
  ifelse(
    as.character(colData(cds)$RNA_snn_res.0.2) == "10",
    "Atlas cluster 10",
    "Other beta cells"
  )


###############################################################################
# 7. Cluster-label coordinates
###############################################################################

cluster_centers <- cell_umap %>%
  group_by(cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )


###############################################################################
# 8. Supplementary Fig S3G
#    Reintegrated beta-cell topology
###############################################################################

cluster_palette <- c(
  "1" = "#F8766D",
  "2" = "#C49A00",
  "3" = "#53B400",
  "4" = "#00C094",
  "5" = "#00B6EB",
  "6" = "#A58AFF"
)

lineage_palette <- c(
  "Lineage1" = "blue",
  "Lineage2" = "red",
  "Lineage3" = "#2B2B2B"
)


panel_A <- ggplot(
  cell_umap,
  aes(
    x = UMAP_1,
    y = UMAP_2
  )
) +
  geom_point(
    aes(color = cluster),
    size = 0.35,
    alpha = 0.8
  ) +
  geom_path(
    data = trajectory_data,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = lineage,
      group = lineage
    ),
    linewidth = 1.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = cluster_centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = cluster
    ),
    fontface = "bold",
    size = 5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      cluster_palette,
      lineage_palette
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1,
        linewidth = 1.5
      )
    )
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  )


###############################################################################
# 9. Supplementary Fig S3H
#    Location of cells originating from atlas cluster 10
###############################################################################

panel_B <- ggplot(
  cell_umap,
  aes(
    x = UMAP_1,
    y = UMAP_2
  )
) +
  geom_point(
    aes(color = atlas_cluster10),
    size = 0.35,
    alpha = 0.8
  ) +
  geom_path(
    data = trajectory_data,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = lineage,
      group = lineage
    ),
    linewidth = 1.3,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Other beta cells" = "grey80",
      "Atlas cluster 10" = "black",
      lineage_palette
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1,
        linewidth = 1.5
      )
    )
  ) + 
  coord_equal() +
  theme_classic() +
  theme(
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )+
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  )


###############################################################################
# 10. Combine and save Supplementary Fig S3G_S3H
###############################################################################

FigS3G_S3H <- panel_A | panel_B

FigS3G_S3H


###############################################################################
# 11. Quantify where atlas cluster 10 maps after reintegration
###############################################################################

cluster10_distribution <- cell_umap %>%
  filter(
    atlas_cluster10 == "Atlas cluster 10"
  ) %>%
  dplyr::count(
    cluster,
    name = "n"
  ) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  arrange(
    desc(n)
  )

print(
  cluster10_distribution
)


###############################################################################
# 12. Source data
###############################################################################

write_tsv(
  cell_umap,
  file.path(outdir, "FigS3G_S3H_cell_coordinates.tsv")
)

write_tsv(
  trajectory_data,
  file.path(outdir, "FigS3G_S3H_trajectory_curves.tsv")
)

write_tsv(
  cluster10_distribution,
  file.path(outdir, "FigS3G_S3H_cluster10_distribution.tsv")
)


###############################################################################
# sessionInfo()
###############################################################################



