###############################################################################
# Project: Pancreas scRNA atlas — Figure 2 (alpha cells)
# File: Figure2_alpha_monocle.R
# Author: Elhadi Iich
# Date: 2026-09-07
#
# Description:
#   Generates alpha-cell supplemental panels S3A-S3E:
#   S3A: UMAP by cluster
#   S3B: Age and BMI distributions by cluster
#   S3C: Study composition
#   S3D: Condition composition
#   S3E: GCG and PPY expression
#
# Inputs:
#   - alpha.integrated.rds
#   - cds_alpha.rds
#
# Dependencies:
#   Seurat, patchwork, reticulate, monocle3, slingshot, dplyr, tidyr, ggplot2,
#   ggpubr (used via ::), Matrix (implicitly via Seurat), grid
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(reticulate)
  library(monocle3)
  library(slingshot)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(123)

# Ensure python leidenalg is available
reticulate::py_module_available(module = "leidenalg")
reticulate::import("leidenalg")

palette <- c("#5A5156FF", "#E4E1E3FF", "#F6222EFF",
             "#FE00FAFF", "#16FF32FF", "#3283FEFF",
             "#FEAF16FF", "#B00068FF", "#1CFFCEFF",
             "#90AD1CFF", "#2ED9FFFF", "#DEA0DFFF",
             "#AA0EDFFF", "#F8A19FFF", "#325A9BFF",
             "#C4451CFF", "#1C8356FF", "#85660DFF",
             "#B10DA1FF", "#FBE426FF", "#1CBE4FFF")

cores <- 2
resolution <- 0.00001

# ----------------------------------------------------------------------------- 
# Load integrated object
# -----------------------------------------------------------------------------
data.integrated <- readRDS("alpha.integrated.rds")
dim(data.integrated)

cds <- readRDS('cds_alpha.rds')

message("######################## extract data for cds ######################################")

# expr <- as(as.matrix(data.integrated@assays$SCT@counts), "sparseMatrix")
# gene_metadata <- data.frame(gene_short_name = rownames(expr),
#                             row.names = rownames(expr))
# cell_metadata <- as.data.frame(data.integrated@meta.data)
# 
# cds <- new_cell_data_set(expr,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_metadata)
# 
# colnames(pData(cds))[10] <- "sample_name2"
# fData(cds)$gene_short_name <- rownames(fData(cds))
# 
# message("################################ assign partitions #########################################")
# reacreate.partition <- as.factor(rep(1, length(cds@colData@rownames)))
# names(reacreate.partition) <- cds@colData@rownames
# cds@clusters$UMAP$partitions <- reacreate.partition
# 
# meta <- pData(cds)
# pData(cds) <- meta
# 
# rm(expr); gc()
# 
# # Use Seurat PCA embeddings (first 10 PCs)
# reducedDims(cds)$PCA <- data.integrated@reductions$pca@cell.embeddings
# reducedDims(cds)$PCA <- reducedDims(cds)$PCA[, 1:10]
# 
# # Ensure character columns are factors
# meta <- as.data.frame(colData(cds))
# meta <- dplyr::mutate_if(meta, is.character, as.factor)
# 
# # UMAP + Leiden
# cds <- reduce_dimension(cds, cores = cores, umap.fast_sgd = FALSE,
#                         reduction_method = "UMAP")
# cds <- cluster_cells(cds, verbose = TRUE, resolution = resolution,
#                      cluster_method = "leiden")
# 
# colData(cds)$monocle_clusters <- cds@clusters@listData$UMAP$clusters


data.integrated@meta.data$monocle_clusters <- pData(cds)$monocle_clusters

df <- data.integrated@meta.data %>% as.data.frame()
df$monocle_clusters <- as.character(df$monocle_clusters)

rownames(df) <- rownames(data.integrated@meta.data)

# Order clusters by increasing %T2D
desired_order <- df %>%
  group_by(monocle_clusters, Condition) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  mutate(per = round(100 * n / sum(n), 1)) %>%
  filter(Condition == "T2D") %>%
  arrange(per) %>%
  pull(monocle_clusters)


# Gene subsets for violins
cds_subset_GCG  <- cds[rownames(cds) == "GCG", ]
cds_subset_PPY  <- cds[rownames(cds) == "PPY", ]

# saveRDS(cds, 'cds_alpha.rds')

################################################################################
## Figure S3A — UMAP by cluster
################################################################################
A3 <- plot_cells(cds, show_trajectory_graph = FALSE,
                 color_cells_by = "cluster",
                 label_cell_groups = TRUE, cell_size = 1,
                 group_cells_by = "cluster",
                 group_label_size = 10) +
  ggpubr::clean_theme()
A3
################################################################################
## Figure S3B-1 — Age distribution by cluster
################################################################################

B3_1 <- df %>%
  tidyr::drop_na(Age) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  ggplot(
    aes(
      x = monocle_clusters,
      y = Age,
      fill = monocle_clusters
    )
  ) +
  geom_violin(
    scale = "width",
    trim = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  ggpubr::theme_pubr(
    border = FALSE,
    legend = "none",
    base_size = 10
  ) +
  scale_x_discrete(limits = desired_order) +
  labs(
    x = NULL,
    y = "Age (years)"
  ) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  )

################################################################################
## Figure S3B-2 — BMI distribution by cluster
################################################################################

B3_2 <- df %>%
  tidyr::drop_na(BMI) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  ggplot(
    aes(
      x = monocle_clusters,
      y = BMI,
      fill = monocle_clusters
    )
  ) +
  geom_violin(
    scale = "width",
    trim = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  ggpubr::theme_pubr(
    border = FALSE,
    legend = "none",
    base_size = 10
  ) +
  scale_x_discrete(limits = desired_order) +
  labs(
    x = "Cluster",
    y = expression(bold("BMI (kg/m"^2*")"))
  ) +
  theme(
    axis.title.y = element_text(size = 16),
    axis.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  )
################################################################################
## Figures S3A–S3B combined
################################################################################
A3 / (B3_1 / B3_2)

################################################################################
## Figure S3C — Study composition per cluster (proportion)
################################################################################
C3 <- df %>%
  group_by(monocle_clusters, study) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  mutate(per = paste0(round(100 * n / sum(n), 1), "%")) %>%
  ggplot(aes(x = monocle_clusters, y = n, fill = study)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = "right", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  ylab("Proportion") + xlab(NULL) +
  scale_fill_manual(values = palette[-c(12, 13, 14)]) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text   = element_text(face = "bold", size = 12),
    legend.title  = element_text(face = "bold", size = 14)
  )

################################################################################
## Figure S3D — Condition composition per cluster (with labels)
################################################################################
D3 <- df %>%
  group_by(monocle_clusters, Condition) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  mutate(per = paste0(round(100 * n / sum(n), 1), "%")) %>%
  ggplot(aes(x = monocle_clusters, y = n, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = "right", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  ylab("Proportion") + xlab("Cluster") +
  geom_text(aes(label = per), size = 5, fontface = "bold",
            position = position_fill(vjust = 0.5)) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text   = element_text(face = "bold", size = 12),
    legend.title  = element_text(face = "bold", size = 14)
  )

################################################################################
## Figure S3E — Violin plots (GCG, PPY)
################################################################################
E3_1 <- plot_genes_violin(cds_subset_GCG, group_cells_by = "monocle_clusters") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_discrete(limits = desired_order) + ggpubr::theme_pubr() +
  theme(strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(face = "bold", size = 14)) +
  xlab("Cluster") + ylab("Log10(expression)")

E3_2 <- plot_genes_violin(cds_subset_PPY, group_cells_by = "monocle_clusters") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_discrete(limits = desired_order) + ggpubr::theme_pubr() +
  theme(strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(face = "bold", size = 14)) +
  xlab("Cluster") + ylab("Log10(expression)")

# Show E panels
E3_1 / E3_2

###############################################################################
# SOURCE DATA EXPORT — Figure S3A–S3E
###############################################################################

source_dir <- "../SourceData/Figure2_S3"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Figures S3A and S3B
# Alpha-cell UMAP + age/BMI distributions
###############################################################################

s3_umap <- as.data.frame(
  SingleCellExperiment::reducedDims(cds)$UMAP
)

colnames(s3_umap)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

s3_umap <- s3_umap %>%
  tibble::rownames_to_column("cell_id")


s3_meta <- df %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    monocle_clusters,
    study,
    Condition,
    Age,
    BMI
  )


s3ab_source <- s3_umap %>%
  dplyr::left_join(
    s3_meta,
    by = "cell_id"
  )

readr::write_tsv(
  s3ab_source,
  file.path(
    source_dir,
    "FigS3A_B_alpha_cell_data.tsv.gz"
  )
)


###############################################################################
# Figures S3C and S3D
# Study and condition composition
###############################################################################

s3c_source <- s3_meta %>%
  dplyr::count(
    monocle_clusters,
    study,
    name = "cell_count"
  ) %>%
  dplyr::group_by(monocle_clusters) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count),
    variable = "Study",
    category = as.character(study)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = monocle_clusters,
    variable,
    category,
    cell_count,
    proportion
  )


s3d_source <- s3_meta %>%
  dplyr::count(
    monocle_clusters,
    Condition,
    name = "cell_count"
  ) %>%
  dplyr::group_by(monocle_clusters) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count),
    variable = "Condition",
    category = as.character(Condition)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = monocle_clusters,
    variable,
    category,
    cell_count,
    proportion
  )


readr::write_tsv(
  dplyr::bind_rows(
    s3c_source,
    s3d_source
  ),
  file.path(
    source_dir,
    "FigS3C_D_alpha_composition.tsv"
  )
)


###############################################################################
# Figure S3E
# Exact data used by GCG and PPY violin plots
###############################################################################

s3e_source <- dplyr::bind_rows(
  as.data.frame(E3_1$data) %>%
    dplyr::mutate(panel_gene = "GCG"),
  
  as.data.frame(E3_2$data) %>%
    dplyr::mutate(panel_gene = "PPY")
)

readr::write_tsv(
  s3e_source,
  file.path(
    source_dir,
    "FigS3E_GCG_PPY_violin.tsv.gz"
  )
)

################################################################################
# sessionInfo()
################################################################################
