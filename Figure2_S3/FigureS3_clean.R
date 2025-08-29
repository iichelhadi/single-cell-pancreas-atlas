###############################################################################
# Project: Pancreas scRNA atlas — Figure 2 (alpha cells)
# File: Figure2_alpha_monocle.R
# Author: Elhadi Iich
# Date: 2025-08-29
#
# Description:
#   Reproducible script to build a Monocle3 CellDataSet from an integrated
#   Seurat object, run UMAP + Leiden clustering (via python/reticulate),
#   add clinical metadata (Age, BMI, HbA1c), and generate panels A–F:
#     A) UMAP by cluster
#     B–C) Median Age / BMI per cluster (+MAD)
#     D) Study composition per cluster (stacked proportion)
#     D (alt)) Condition composition per cluster (stacked proportion with labels)
#     E) Violin plots for GCG and PPY
#     F) Median HbA1c per cluster with clinical class
#
# Inputs:
#   - alpha.integrated.rds
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
data.integrated <- readRDS("/hpctmp/eiich01/V2/part2/20241030/part4/alpha/alpha.integrated.rds")
dim(data.integrated)

message("######################## extract data for cds ######################################")

expr <- as(as.matrix(data.integrated@assays$SCT@counts), "sparseMatrix")
gene_metadata <- data.frame(gene_short_name = rownames(expr),
                            row.names = rownames(expr))
cell_metadata <- as.data.frame(data.integrated@meta.data)

cds <- new_cell_data_set(expr,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

colnames(pData(cds))[10] <- "sample_name2"
fData(cds)$gene_short_name <- rownames(fData(cds))

message("################################ assign partitions #########################################")
reacreate.partition <- as.factor(rep(1, length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
cds@clusters$UMAP$partitions <- reacreate.partition

meta <- pData(cds)
pData(cds) <- meta

rm(expr); gc()

# Use Seurat PCA embeddings (first 10 PCs)
reducedDims(cds)$PCA <- data.integrated@reductions$pca@cell.embeddings
reducedDims(cds)$PCA <- reducedDims(cds)$PCA[, 1:10]

# Ensure character columns are factors
meta <- as.data.frame(colData(cds))
meta <- dplyr::mutate_if(meta, is.character, as.factor)

# UMAP + Leiden
cds <- reduce_dimension(cds, cores = cores, umap.fast_sgd = FALSE,
                        reduction_method = "UMAP")
cds <- cluster_cells(cds, verbose = TRUE, resolution = resolution,
                     cluster_method = "leiden")

colData(cds)$monocle_clusters <- cds@clusters@listData$UMAP$clusters
data.integrated@meta.data$monocle_clusters <- pData(cds)$monocle_clusters

# ----------------------------------------------------------------------------- 
# Add HbA1c-derived classes to metadata
# -----------------------------------------------------------------------------
df <- data.integrated@meta.data %>% as.data.frame()
df2 <- df %>%
  group_by(monocle_clusters) %>%
  summarise(
    mean_value   = round(mean(Hba1c, na.rm = TRUE), 2),
    median_value = median(Hba1c, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Hba1c_class = case_when(
      mean_value < 5.7 ~ "ND",
      mean_value >= 5.7 & mean_value < 6.5 ~ "pre-T2D",
      mean_value >= 6.5 ~ "T2D"
    ),
    Hba1c_class2 = case_when(
      median_value < 5.7 ~ "ND",
      median_value >= 5.7 & median_value < 6.5 ~ "pre-T2D",
      median_value >= 6.5 ~ "T2D"
    )
  )

df <- dplyr::left_join(df, df2, by = "monocle_clusters")

## Add to cds & Seurat meta
colData(cds)$mean_value   <- df$mean_value
colData(cds)$median_value <- df$median_value
colData(cds)$Hba1c_class  <- df$Hba1c_class
colData(cds)$Hba1c_class2 <- df$Hba1c_class2

rownames(df) <- rownames(data.integrated@meta.data)
data.integrated@meta.data <- df

df$monocle_clusters <- as.character(df$monocle_clusters)

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

# Count unique samples per cluster
numbers <- lapply(unique(df$monocle_clusters), function(i) {
  length(unique(df$orig.ident[df$monocle_clusters == i]))
})
numbers_df <- data.frame(Clusters = unique(df$monocle_clusters),
                         n_Samples = unlist(numbers))

# Gene subsets for violins
cds_subset_INS  <- cds[rownames(cds) == "INS", ]
cds_subset_GCG  <- cds[rownames(cds) == "GCG", ]
cds_subset_PPY  <- cds[rownames(cds) == "PPY", ]
cds_subset_SST  <- cds[rownames(cds) == "SST", ]
cds_subset_GHRL <- cds[rownames(cds) == "GHRL", ]

################################################################################
## Figure 2A — UMAP by cluster
################################################################################
A2 <- plot_cells(cds, show_trajectory_graph = FALSE,
                 color_cells_by = "cluster",
                 label_cell_groups = TRUE, cell_size = 1,
                 group_cells_by = "cluster",
                 group_label_size = 10) +
  ggpubr::clean_theme()

################################################################################
## Figure 2B — Age by cluster (median + MAD)
################################################################################
B2_1 <- df %>%
  tidyr::drop_na(Age) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  summarise(
    median_value = median(Age, na.rm = TRUE),
    mad_value    = mad(Age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = monocle_clusters, y = median_value, fill = monocle_clusters)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = median_value, ymax = median_value + mad_value),
                width = 0.2) +
  ggpubr::theme_pubr(border = FALSE, legend = "none", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median Age") +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  ) +
  xlab(NULL)

################################################################################
## Figure 2C — BMI by cluster (median + MAD)
################################################################################
B2_2 <- df %>%
  tidyr::drop_na(BMI) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  summarise(
    median_value = median(BMI, na.rm = TRUE),
    mad_value    = mad(BMI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = monocle_clusters, y = median_value, fill = monocle_clusters)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = median_value, ymax = median_value + mad_value),
                width = 0.2) +
  ggpubr::theme_pubr(border = FALSE, legend = "none", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median BMI") +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  )

################################################################################
## Figures 2A–2C combined
################################################################################
A2 / (B2_1 / B2_2)

################################################################################
## Figure 2D — Study composition per cluster (proportion)
################################################################################
C2 <- df %>%
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
## Figure 2D (alt) — Condition composition per cluster (with labels)
################################################################################
D2 <- df %>%
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
## Figure 2E — Violin plots (GCG, PPY)
################################################################################
E2_1 <- plot_genes_violin(cds_subset_GCG, group_cells_by = "monocle_clusters") +
  scale_x_discrete(limits = desired_order) + ggpubr::theme_pubr() +
  theme(strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(face = "bold", size = 14)) +
  xlab("Cluster") + ylab("Log10(expression)")

E2_2 <- plot_genes_violin(cds_subset_PPY, group_cells_by = "monocle_clusters") +
  scale_x_discrete(limits = desired_order) + ggpubr::theme_pubr() +
  theme(strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(face = "bold", size = 14)) +
  xlab("Cluster") + ylab("Log10(expression)")

# Show E panels
E2_1 / E2_2

################################################################################
## Figure 2F — Median HbA1c per cluster (+MAD) with class
################################################################################
F2 <- df %>%
  tidyr::drop_na(Hba1c) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  summarise(
    median_value = median(Hba1c, na.rm = TRUE),
    mad_value    = mad(Hba1c, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Hba1c_class = case_when(
      median_value < 5.7 ~ "ND",
      median_value >= 5.7 & median_value < 6.5 ~ "pre-T2D",
      median_value >= 6.5 ~ "T2D"
    )
  ) %>%
  ggplot(aes(x = monocle_clusters, y = median_value, fill = Hba1c_class)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = median_value, ymax = median_value + mad_value),
                width = 0.2) +
  ylim(0, 9) +
  geom_text(aes(label = paste0(round(median_value, 1), "%")),
            size = 7, vjust = 2, fontface = "bold",
            position = position_dodge(0.8)) +
  ggpubr::theme_pubr(border = FALSE, legend = "top", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median HbA1c(%)") +
  guides(fill = guide_legend(title = "HbA1c class")) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text   = element_text(face = "bold", size = 10),
    legend.title  = element_text(face = "bold", size = 14)
  )

# Print F2 if running interactively
# print(F2)

################################################################################
# sessionInfo()
################################################################################
