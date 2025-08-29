###############################################################################
# Project: scRNA beta cell atlas — Lineage 1 reintegration and Figure 3 suite
# File: Figure3_FigureS5.R
# Author: Elhadi Iich
# Date: 2025-08-27
#
# Description:
#   End-to-end workflow for lineage 1 reintegration and visualization:
#   reclustering, supplementary figures, Slingshot overlays, INS trends,
#   branched heatmaps, GO enrichment, module detection and trajectories,
#   pseudobulk DESeq2, module scores along pseudotime, dot plots, and
#   metadata associations.
#
# Inputs:
#   - Lin1.beta.reintegrated.rds
#   - cds_lin1.rds
#   - slingshot_sim_lin1.rds
#   - Beta_T2D_lineage1_1.rds, Beta_T2D_lineage1_2.rds
#   - Lineage1_sling_perc_expr.tsv
#   - ck1/ck2 RDS objects for branched heatmaps
#   - pseudobulk cts/coldata TSVs for DESeq2
#
###############################################################################

################################################################################
### 1. Environment Setup & Library Loading
################################################################################


suppressPackageStartupMessages({
  # Core tidy data manipulation and visualization
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(scales)

  # String and text processing
  library(stringi)
  library(stringr)

  # Single-cell analysis
  library(Seurat)
  library(patchwork)
  library(reticulate)
  library(monocle3)
  library(slingshot)

  # Enrichment analysis
  library(clusterProfiler)
  library(org.Hs.eg.db)

  # Heatmaps and visualization
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(viridis)

  # Data handling and parallel computing
  library(data.table)
  library(parallel)

  # Utilities
  library(rlang)
})

set.seed(123)
reticulate::py_module_available(module = "leidenalg")
reticulate::import("leidenalg")

palette <- c(
  "#5A5156FF", "#E4E1E3FF", "#F6222EFF", "#FE00FAFF",
  "#16FF32FF", "#3283FEFF", "#FEAF16FF", "#B00068FF",
  "#1CFFCEFF", "#90AD1CFF", "#2ED9FFFF", "#DEA0DFFF",
  "#AA0EDFFF", "#F8A19FFF", "#325A9BFF", "#C4451CFF",
  "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF",
  "#1CBE4FFF"
)

################################################################################
### 2. Data Loading
################################################################################

data.integrated <- readRDS("Lin1.beta.reintegrated.rds")
df              <- data.integrated@meta.data %>% as.data.frame()
cds             <- readRDS("cds_lin1.rds")

desired_order <- c("2", "1", "3", "4")

################################################################################
### 3. Figure 3B–3D
################################################################################

## Figure 3B
B3 <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = "monocle_clusters2",
  label_cell_groups = TRUE,
  cell_size = 1,
  group_cells_by = "cluster",
  group_label_size = 10
) + ggpubr::clean_theme()
print(B3)

## Figure 3C
C3 <- df %>%
  drop_na(Age) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(
    median_value = median(Age, na.rm = TRUE),
    mad_value    = mad(Age, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = monocle_clusters2)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(),
    color    = "black",
    width    = 0.8
  ) +
  geom_errorbar(
    aes(ymin = median_value, ymax = median_value + mad_value),
    width    = 0.2,
    position = position_dodge(0.9)
  ) +
  ggpubr::theme_pubr(border = FALSE, legend = "none", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median Age") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 14)
  )
print(C3)

## Figure 3D
D3 <- df %>%
  drop_na(BMI) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(
    median_value = median(BMI, na.rm = TRUE),
    mad_value    = mad(BMI, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = monocle_clusters2)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(),
    color    = "black",
    width    = 0.8
  ) +
  geom_errorbar(
    aes(ymin = median_value, ymax = median_value + mad_value),
    width    = 0.2,
    position = position_dodge(0.9)
  ) +
  ggpubr::theme_pubr(border = FALSE, legend = "none", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median BMI") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(face = "bold", size = 14)
  )
print(D3)

## Combine 3C + 3D
print(C3 / D3)

################################################################################
### 4. Supplementary Figure S5A–S5F
################################################################################

## S5A
S5A <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = "Study",
  label_cell_groups = FALSE,
  cell_size = 1,
  group_cells_by = "cluster",
  group_label_size = 4
) +
  ggpubr::clean_theme() + NoLegend() +
  scale_color_manual(values = palette[-c(12, 13, 14)])

## S5B
S5B <- df %>%
  group_by(monocle_clusters2, Study) %>%
  dplyr::count() %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  mutate(per = paste0(round(100 * n / sum(n), 1), "%")) %>%
  ggplot(aes(x = monocle_clusters2, y = n, fill = Study)) +
  geom_bar(
    stat = "identity",
    position = "fill",
    color = "black",
    width = 0.8
  ) +
  ggpubr::theme_pubr(border = FALSE, legend = "top", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  ylab("Proportion") + xlab("Cluster") +
  scale_fill_manual(values = palette[-c(12, 13, 14)]) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 14)
  )

## S5C
S5C <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = "Condition",
  label_cell_groups = FALSE,
  cell_size = 1,
  group_cells_by = "cluster",
  group_label_size = 4
) + ggpubr::clean_theme() + NoLegend()

## S5D
S5D <- df %>%
  group_by(monocle_clusters2, Condition) %>%
  dplyr::count() %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  mutate(per = paste0(round(100 * n / sum(n), 1), "%")) %>%
  ggplot(aes(x = monocle_clusters2, y = n, fill = Condition)) +
  geom_bar(
    stat = "identity",
    position = "fill",
    color = "black",
    width = 0.8
  ) +
  ggpubr::theme_pubr(border = FALSE, legend = "top", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  ylab("Proportion") + xlab("Cluster") +
  geom_text(
    aes(label = per),
    size = 8,
    position = position_fill(vjust = 0.5)
  ) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14)
  )

print((S5A / S5B) | (S5C / S5D))

## S5E
S5E <- df %>%
  drop_na(Hba1c) %>%
  ggplot(aes(monocle_clusters2, Hba1c, fill = monocle_clusters2)) +
  geom_boxplot() +
  ggpubr::theme_pubr(border = FALSE, legend = "none") +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "HbA1c (%)") +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14)
  )
print(S5E)

## S5F
S5F <- df %>%
  drop_na(Hba1c) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(
    median_value = median(Hba1c, na.rm = TRUE),
    mad_value    = mad(Hba1c, na.rm = TRUE)
  ) %>%
  mutate(
    Hba1c_class = case_when(
      median_value < 5.7                        ~ "ND",
      median_value >= 5.7 & median_value < 6.5 ~ "pre-T2D",
      median_value >= 6.5                       ~ "T2D"
    )
  ) %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = Hba1c_class)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(),
    color    = "black",
    width    = 0.8
  ) +
  geom_errorbar(aes(ymin = median_value, ymax = median_value + mad_value), width = 0.2) +
  ylim(0, 7.5) +
  scale_fill_manual(values = c("#F8766D", "#FEAF16FF", "#619CFF")) +
  geom_text(
    aes(label = paste0(round(median_value, 1), "%")),
    size = 8,
    vjust = 2,
    position = position_dodge(0.8)
  ) +
  ggpubr::theme_pubr(border = FALSE, legend = "top", base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median HbA1c (%)") +
  guides(fill = guide_legend(title = "HbA1c class")) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14)
  )
print(S5F)

################################################################################
### 5. Slingshot Curves & UMAP Labels
################################################################################

sim <- readRDS("slingshot_sim_lin1.rds")

coldata <- data.frame(
  study        = sim$study,
  Cluster      = sim$monocle_clusters2,
  Hba1c_class2 = sim$Hba1c_class_median,
  Lineage1_1   = sim$slingPseudotime_1,
  Lineage1_2   = sim$slingPseudotime_2
)
rownames(coldata) <- colnames(sim)

lineages <- slingshot::slingCurves(SlingshotDataSet(sim))
lineage_data <- lapply(names(lineages), function(lineage_name) {
  curve_df <- as.data.frame(lineages[[lineage_name]]$s)
  colnames(curve_df) <- c("Dim_1", "Dim_2")
  curve_df$lineage <- lineage_name
  curve_df$path_id <- paste(lineage_name, seq_len(nrow(curve_df)), sep = "_")
  curve_df
})
combined_lineages <- do.call(rbind, lineage_data)
combined_lineages$lineage <- ifelse(
  combined_lineages$lineage == "Lineage1",
  "Lineage1_1",
  "Lineage1_2"
)

umap_coords <- reducedDims(cds)$UMAP
umap_df <- data.frame(
  Dim_1   = umap_coords[, 1],
  Dim_2   = umap_coords[, 2],
  cluster = clusters(cds, reduction_method = "UMAP")
)
label_coords <- aggregate(cbind(Dim_1, Dim_2) ~ cluster, data = umap_df, FUN = mean)

B3_slingshot <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by          = "cluster",
  label_groups_by_cluster = FALSE,
  label_cell_groups       = FALSE,
  cell_size               = 1,
  group_label_size        = 10
) +
  geom_path(
    data = combined_lineages,
    aes(x = Dim_1, y = Dim_2, color = lineage, group = lineage),
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "1"          = "#F8766D",
      "2"          = "#53B400",
      "3"          = "#00bec4",
      "4"          = "#c67bff",
      "Lineage1_1" = "#325A9BFF",
      "Lineage1_2" = "#2ED9FFFF"
    )
  ) +
  geom_text(data = label_coords, aes(x = Dim_1, y = Dim_2, label = cluster), size = 10) +
  ggpubr::clean_theme() +
  theme(
    legend.text  = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14)
  )
print(B3_slingshot)

################################################################################
### 6. Slingshot Branch Plots & INS Expression (Figures 3G–3H)
################################################################################

lin1_1 <- readRDS(
  "Beta_T2D_lineage1_1.rds"
)
dim(lin1_1)

lin1_2 <- readRDS(
  "Beta_T2D_lineage1_2.rds"
)
dim(lin1_2)

###### Figure 3G
exp1_1 <- lin1_1@assays$RNA@data %>% t()

meta <- lin1_1@meta.data[, c("Lineage1_1", "monocle_clusters2")] %>% as.data.frame()
colnames(meta)[1:2] <- c("Lineage1_1", "Clusters")
all(rownames(meta) %in% rownames(exp1_1))

plt_data <- cbind(meta, exp1_1)

mycolors <- c("#F8766D", "#53B400", "#00bec4")

G3 <- ggscatter(
  data = plt_data,
  x = "Lineage1_1",
  y = "INS",
  color = "Clusters",
  size = 0.6
) +
  geom_smooth(se = FALSE, color = "black", method = "lm",
              formula = y ~ splines::ns(x, df = 3)) +
  ggpubr::theme_pubr() +
  scale_color_manual(values = mycolors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  xlab("Pseudotime") +
  theme(
    legend.position = "right",
    axis.title.y  = element_text(face = "bold", size = 16),
    axis.title.x  = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 16),
    legend.title  = element_text(face = "bold", size = 14),
    legend.text   = element_text(face = "bold", size = 14, hjust = 1, vjust = 1)
  )
print(G3)

###### Figure 3H
exp1_2 <- lin1_2@assays$RNA@data %>% t()

meta <- lin1_2@meta.data[, c("Lineage1_2", "monocle_clusters2")] %>% as.data.frame()
colnames(meta)[1:2] <- c("Lineage1_2", "Clusters")
all(rownames(meta) %in% rownames(exp1_2))

plt_data <- cbind(meta, exp1_2)

mycolors <- c("#F8766D", "#53B400", "#A58AFF")

H3 <- ggscatter(
  data = plt_data,
  x = "Lineage1_2",
  y = "INS",
  color = "Clusters",
  size = 0.6
) +
  geom_smooth(se = FALSE, color = "black", method = "lm",
              formula = y ~ splines::ns(x, df = 3)) +
  ggpubr::theme_pubr() +
  scale_color_manual(values = mycolors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  xlab("Pseudotime") +
  theme(
    legend.position = "right",
    axis.title.y  = element_text(face = "bold", size = 16),
    axis.title.x  = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 16),
    legend.title  = element_text(face = "bold", size = 14),
    legend.text   = element_text(face = "bold", size = 14, hjust = 1, vjust = 1)
  )
print(H3)

################################################################################
### 7. Heatmap Data Prep, Drawing & Legends (Figure 3I)
################################################################################

# 1) Source custom helper functions
source("custom_ClusterData.R")
source("custom_visCluster.R")

# 2) Define genes of interest
genes_of_interest <- c(
  "INS", "MAFA", "ABCC8", "IAPP", "G6PC2", "PAX6", "NKX6-1", "PDX1", "SLC30A8",
  "BHLHE41", "HERPUD1", "DDIT3", "HSP90B1", "TNFRSF12A", "VGF",
  "HSPA1A", "HSPA1B", "HSPA5", "GCG", "FXYD2"
)

coldata <- data.frame(
  study        = sim$study,
  Cluster      = sim$monocle_clusters2,
  Hba1c_class2 = sim$Hba1c_class_median,
  Lineage1_1   = sim$slingPseudotime_1,
  Lineage1_2   = sim$slingPseudotime_2
)
rownames(coldata) <- colnames(sim)

mat1_raw <- lin1_1@assays$SCT@data

# 4) Branch 1: order → smooth (parallel) → z-score ---------------------
cells1     <- intersect(colnames(mat1_raw), rownames(coldata))
ord1       <- order(coldata[cells1, "Lineage1_1"])
cells1_o   <- cells1[ord1]
cells1_rev <- rev(cells1_o)
pt1        <- coldata[cells1_rev, "Lineage1_1"]
mat1_ord   <- mat1_raw[, cells1_rev]

# 1) Define a “never-fail” smooth helper
safe_smooth <- function(x, df = 3) {
  tryCatch(
    stats::smooth.spline(x, df = df)$y,
    error = function(e) rep(NA_real_, length(x))
  )
}

ncores <- min(detectCores(), 3)

# Branch 1
sm1_list <- mclapply(
  seq_len(nrow(mat1_ord)),
  function(i) safe_smooth(mat1_ord[i, ], df = 3),
  mc.cores = ncores
)
gc()
sm1 <- do.call(rbind, sm1_list)
rownames(sm1) <- rownames(mat1_ord)
mat1_z <- t(scale(t(sm1)))
colnames(mat1_z) <- colnames(mat1_ord)
# saveRDS(mat1_z, "mat1_z.rds") ### mat1_z generated from RNA assay for module detection and SCT for heatmap

# 5) Branch 2: same steps -----------------------------------------------
mat2_raw <- lin1_2@assays$SCT@data
# mat2_raw <- mat2_raw[debug_genes, ]

cells2   <- intersect(colnames(mat2_raw), rownames(coldata))
ord2     <- order(coldata[cells2, "Lineage1_2"])
cells2_o <- cells2[ord2]
pt2      <- coldata[cells2_o, "Lineage1_2"]
mat2_ord <- mat2_raw[, cells2_o]

sm2_list <- mclapply(
  seq_len(nrow(mat2_ord)),
  function(i) safe_smooth(mat2_ord[i, ], df = 3),
  mc.cores = ncores
)
gc()
sm2 <- do.call(rbind, sm2_list)
rownames(sm2) <- rownames(mat2_ord)
mat2_z <- t(scale(t(sm2)))
colnames(mat2_z) <- colnames(mat2_ord)
# saveRDS(mat2_z, "mat2_z.rds") ### mat2_z generated from RNA assay for module detection and SCT for heatmap

# 6) Cluster rows -------------------------------------------------------
# ck1 <- clusterData_custom(exp = mat1_z, scaleData = TRUE, min.std = 0)
# ck2 <- clusterData_custom(exp = mat2_z, scaleData = TRUE, min.std = 0)


# 3) Load precomputed clustering results
ck1 <- readRDS("ck1_branched_heatmap.rds")
ck2 <- readRDS("ck2_branched_heatmap.rds")

# 4) Harmonize rows & build common matrices
shared     <- intersect(ck1$wide.res$gene, ck2$wide.res$gene)
common_all <- unique(c(shared, genes_of_interest))
row_order  <- ck1$wide.res$gene[ck1$wide.res$gene %in% common_all]

mat1_common <- ck1$wide.res %>%
  filter(gene %in% row_order) %>%
  arrange(match(gene, row_order)) %>%
  dplyr::select(-cluster) %>%
  column_to_rownames("gene") %>%
  as.matrix()

mat2_common <- ck2$wide.res %>%
  filter(gene %in% row_order) %>%
  arrange(match(gene, row_order)) %>%
  dplyr::select(-cluster) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# 5) Build top annotations
library(circlize)
pt1_fun <- colorRamp2(seq(min(pt1), max(pt1), length.out = 100), plasma(100))
pt2_fun <- colorRamp2(seq(min(pt2), max(pt2), length.out = 100), plasma(100))

col_annot1 <- HeatmapAnnotation(
  HbA1c = anno_simple(
    coldata[cells1_rev, "Hba1c_class2"],
    col = setNames(c("#F8766D", "#FEAF16FF", "#619CFF"),
                   sort(unique(coldata[cells1_rev, "Hba1c_class2"])))
  ),
  Cluster = anno_simple(
    coldata[cells1_rev, "Cluster"],
    col = setNames(c("#F8766D", "#53B400", "#00B6EB"),
                   sort(unique(coldata[cells1_rev, "Cluster"])))
  ),
  Pseudotime = anno_simple(pt1, col = pt1_fun),
  Lineage    = anno_simple(rep("Lineage1_1", length(pt1)),
                           col = setNames("#325A9BFF", "Lineage1_1")),
  annotation_height  = unit(c(4, 4, 4, 4), "mm"),
  annotation_label   = "",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

col_annot2 <- HeatmapAnnotation(
  HbA1c = anno_simple(
    coldata[cells2_o, "Hba1c_class2"],
    col = setNames(c("#F8766D", "#FEAF16FF", "#619CFF"),
                   sort(unique(coldata[cells2_o, "Hba1c_class2"])))
  ),
  Cluster = anno_simple(
    coldata[cells2_o, "Cluster"],
    col = setNames(c("#F8766D", "#53B400", "#A58AFF"),
                   sort(unique(coldata[cells2_o, "Cluster"])))
  ),
  Pseudotime = anno_simple(pt2, col = pt2_fun),
  Lineage    = anno_simple(rep("Lineage1_2", length(pt2)),
                           col = setNames("#2ED9FFFF", "Lineage1_2")),
  annotation_height  = unit(c(4, 4, 4, 4), "mm"),
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

# 6) Row markers on the right heatmap
pos_interest <- which(rownames(mat2_common) %in% genes_of_interest)
row_mark <- rowAnnotation(
  mark = anno_mark(
    at        = pos_interest,
    labels    = rownames(mat2_common)[pos_interest],
    labels_gp = gpar(fontsize = 14)
  )
)

# 7) Draw heatmaps side by side
col_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "#00FFFF", "#FF0000"))

ht1 <- Heatmap(
  mat1_common, name = NULL, col = col_fun,
  top_annotation      = col_annot1,
  cluster_rows        = FALSE,
  cluster_columns     = FALSE,
  show_row_names      = FALSE,
  show_column_names   = FALSE,
  show_row_dend       = FALSE,
  show_heatmap_legend = FALSE,
  border              = TRUE
)

ht2 <- Heatmap(
  mat2_common, name = "Expression (z)", col = col_fun,
  top_annotation      = col_annot2,
  right_annotation    = row_mark,
  cluster_rows        = FALSE,
  cluster_columns     = FALSE,
  show_row_names      = FALSE,
  show_column_names   = FALSE,
  show_row_dend       = FALSE,
  show_heatmap_legend = TRUE,
  border              = TRUE
)

# Cairo::CairoPNG("20250807_combined_heatmaps4.png",
#   width = 10 * 300, height = 8 * 300, res = 300)
draw(ht1 + ht2, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()

# 8) Build legends
library(ComplexHeatmap)

# HbA1c legend
HbA1c_classes <- sort(unique(c(coldata[cells1_rev, "Hba1c_class2"], coldata[cells2_o, "Hba1c_class2"])))
HbA1c_colors  <- setNames(c("#F8766D", "#FEAF16FF", "#619CFF")[seq_along(HbA1c_classes)], HbA1c_classes)
HbA1c_legend  <- Legend(labels = names(HbA1c_colors), title = "HbA1c Class", legend_gp = gpar(fill = HbA1c_colors))

# Cluster legend
clusters_all   <- sort(unique(c(coldata[cells1_rev, "Cluster"], coldata[cells2_o, "Cluster"])))
cluster_colors <- setNames(c("#F8766D", "#53B400", "#00B6EB", "#A58AFF")[seq_along(clusters_all)], clusters_all)
Cluster_legend <- Legend(labels = names(cluster_colors), title = "Cluster", legend_gp = gpar(fill = cluster_colors))

# Lineage legend
lineage_colors <- c("Lineage1_1" = "#325A9BFF", "Lineage1_2" = "#2ED9FFFF")
Lineage_legend <- Legend(labels = names(lineage_colors), title = "Lineage", legend_gp = gpar(fill = lineage_colors))

# Pseudotime legend
pseudotime_fun <- colorRamp2(seq(min(pt1, pt2), max(pt1, pt2), length.out = 100), plasma(100))
Pseudotime_legend <- Legend(
  col_fun = pseudotime_fun, title = "Pseudotime",
  legend_height = unit(3, "cm"),
  at = round(c(0, 5, 10, 15, 20), 1)
)

# Expression legend
expr_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "#00FFFF", "#FF0000"))
Expression_legend <- Legend(
  col_fun = expr_fun, title = "Expression (z)",
  legend_height = unit(2, "cm")
)

# 9) Pack and draw all legends
all_legends <- packLegend(
  HbA1c_legend,
  Cluster_legend,
  Lineage_legend,
  Pseudotime_legend,
  Expression_legend,
  direction = "vertical"
)
draw(all_legends, just = "left")
dev.off()

################################################################################
### 8. Pseudotime Association Volcano Plots (Figures 3J–3K)
################################################################################

sce <- readRDS("sce_lin1_1.rds")
perc_expr <- read_tsv("Lineage1_sling_perc_expr.tsv")
colnames(perc_expr) <- c("feature_id", "perc_expr")

pseudotime_association         <- rowData(sce)$assocRes
pseudotime_start_end_association <- tradeSeq::startVsEndTest(sce)
pseudotime_start_end_association$padj <- p.adjust(
  pseudotime_start_end_association$pvalue, method = "fdr"
)
pseudotime_start_end_association$feature_id <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association1_1 <- left_join(
  pseudotime_start_end_association, perc_expr, by = "feature_id"
)

# Genes to highlight for lineage1_1
genes <- c(
  "DBP", "FXYD2", "HSPA5", "HSP90B1", "TRIB3",
  "SLC30A8", "NKX6-1", "MAFA", "INS", "IAPP", "PDX1", "G6PC2", "DDIT3", "CEBPD"
)

create_volcano <- function(data, logFC_col) {
  data <- data %>%
    dplyr::select(padj, !!sym(logFC_col), feature_id, perc_expr) %>%
    mutate(
      significance = case_when(
        padj <= 0.05 & !!sym(logFC_col) >  1 ~ "UP",
        padj <= 0.05 & !!sym(logFC_col) < -1 ~ "DOWN",
        TRUE                                 ~ "NS"
      ),
      neg_log10_padj = ifelse(padj == 0, 20, -log10(padj))
    )

  gene_data  <- filter(data, feature_id %in% genes)
  other_data <- filter(data, !(feature_id %in% genes))

  ggplot() +
    geom_jitter(
      data = other_data,
      aes(x = !!sym(logFC_col), y = neg_log10_padj, size = perc_expr),
      color = "lightgray", alpha = 0.8, width = 0.2, height = 0.2, show.legend = FALSE
    ) +
    geom_jitter(
      data = gene_data,
      aes(x = !!sym(logFC_col), y = neg_log10_padj, color = significance, size = perc_expr),
      alpha = 0.8, width = 0.2, height = 0.2
    ) +
    scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "darkgrey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    ggrepel::geom_text_repel(
      data = gene_data,
      aes(x = !!sym(logFC_col), y = neg_log10_padj, label = feature_id),
      max.overlaps = 40, fontface = "bold", size = 4, box.padding = 0.5, force = 5, nudge_y = 2
    ) +
    xlim(-10, 10) +
    ylim(0, max(data$neg_log10_padj) + 5) +
    labs(x = "Log2 Fold Change", y = "-Log10 FDR", size = "Percent Expression") +
    theme_pubr(legend = "right") +
    theme(
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text    = element_text(face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 14),
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(face = "bold", size = 14)
    )
}

### Figure 3J
volcano_1_1 <- create_volcano(pseudotime_start_end_association1_1, "logFClineage1")
print(volcano_1_1)

### Figure 3K
sce <- readRDS("/hpctmp/eiich01/V2/part2/20241030/part4/beta/lineage1_reintegration_test/association_test/sce_lin1_2.rds")
pseudotime_association         <- rowData(sce)$assocRes
pseudotime_start_end_association <- tradeSeq::startVsEndTest(sce)
pseudotime_start_end_association$padj <- p.adjust(
  pseudotime_start_end_association$pvalue, method = "fdr"
)
pseudotime_start_end_association$feature_id <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association1_2 <- left_join(
  pseudotime_start_end_association, perc_expr, by = "feature_id"
)

# Genes to highlight for lineage1_2
genes <- c("FXYD2", "SLC30A8", "TRIB3", "TNFRSF12A", "NKX6-1", "MAFA",
           "INS", "IAPP", "PDX1", "G6PC2", "CEBPB", "VGF")

volcano_1_2 <- create_volcano(pseudotime_start_end_association1_2, "logFClineage1")
print(volcano_1_2)

################################################################################
### 9. Pseudobulk DESeq2 Analysis (Figures 3L–3M)
################################################################################

library(DESeq2)

cts <- read.table(
  "beta_C3_v_C4_pseudobulk_cts.tsv",
  header = TRUE, sep = "\t", row.names = 1
)
colnames(cts) <- gsub(pattern = "[.]", replacement = "-", x = colnames(cts))

coldata <- read.table(
  "beta_C3_v_C4_pseudobulk_coldata.tsv",
  header = TRUE, sep = "\t", row.names = 1
)
coldata$monocle_clusters2 <- as.character(coldata$monocle_clusters2)

perc_expr <- read_tsv("Lineage1_sling_perc_expr.tsv")
colnames(perc_expr) <- c("feature_id", "perc_expr")

cts <- cts[rownames(cts) %in% perc_expr$feature_id[perc_expr$perc_expr >= 5], ]

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = coldata,
  design    = ~ study + monocle_clusters2
)
size_factors <- coldata$counts / median(coldata$counts)
names(size_factors) <- colnames(dds)
sizeFactors(dds) <- size_factors

dds <- DESeq(dds, parallel = FALSE)
dds <- dds[which(mcols(dds)$betaConv), ]

res <- results(dds, contrast = c("monocle_clusters2", "3", "4"), alpha = 0.05)
res_filt <- as.data.frame(res) %>% filter(complete.cases(.))
res_filt$genes <- rownames(res_filt)
res_filt <- left_join(res_filt, perc_expr, by = c("genes" = "feature_id"))

res_filt_2 <- res_filt
res_filt_2$change <- ifelse(
  (res_filt_2$log2FoldChange > 1 & round(res_filt_2$padj, 2) <= 0.05), "UP",
  ifelse((res_filt_2$log2FoldChange < -1 & round(res_filt_2$padj, 2) <= 0.05), "DOWN", "NS")
)

genes_l <- c(
  "HSPA5", "HERPUD1", "HSP90B1", "BHLHE41", "DDIT3", "SLC30A8",
  "NKX6-1", "PAX6", "MAFA", "INS", "IAPP", "PDX1", "G6PC2"
)

############### Figure 3L
res_filt_2 %>%
  mutate(
    highlight = case_when(
      genes %in% genes_l & change == "UP"   ~ "UP",
      genes %in% genes_l & change == "DOWN" ~ "DOWN",
      genes %in% genes_l & change == "NS"   ~ "NS",
      TRUE ~ NA_character_
    )
  ) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = res_filt_2, aes(size = perc_expr), color = "lightgrey", alpha = 0.8) +
  geom_point(aes(color = highlight, size = perc_expr), alpha = 0.8, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) +
  scale_color_manual(
    values = c("UP" = "red", "DOWN" = "blue", "NS" = "darkgrey"),
    na.translate = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = res_filt_2 %>% filter(genes %in% genes_l),
    aes(label = genes),
    max.overlaps = 40,
    fontface = "bold",
    size = 4,
    box.padding = 0.5
  ) +
  xlim(c(-5, 5)) + ylim(0, 5) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Padj",
    size = "Percent Expression"
  ) +
  guides(color = "none") +
  ggpubr::theme_pubr(legend = "right") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14)
  )

################################################################################
### 10. GO Enrichment Dotplots (Figures 3M & S5G–S5H)
################################################################################

res_filt_2_filt <- res_filt_2[res_filt_2$padj < 0.05, ]
genelist <- setNames(res_filt_2_filt$log2FoldChange, res_filt_2_filt$genes)

set.seed(123)
ego <- enrichGO(
  gene          = names(genelist)[abs(genelist) > 1],
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500
)
ego <- simplify(ego)
ids <- c("GO:0031331", "GO:0048193", "GO:0016236", "GO:0008654", "GO:0007034", "GO:0042594")

dotplot(ego %>% filter(ID %in% ids), showCategory = 10) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "\n", " "))) +
  ylab("Term") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(size = 18),
    axis.text.y  = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14)
  )

## S5G – Lineage1_1 GSEA
pseudotime_start_end_association_filt <- read_tsv(
  "lin1_1_pseudotime_start_end_association_filt.tsv"
)
genelist <- setNames(
  pseudotime_start_end_association_filt$logFClineage1,
  pseudotime_start_end_association_filt$feature_id
)
genelist <- sort(genelist, decreasing = TRUE)
genelist <- genelist[abs(genelist) > 1]

set.seed(123)
gse <- gseGO(
  geneList      = genelist,
  ont           = "BP",
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  eps           = 0,
  minGSSize     = 25,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)
gse <- simplify(gse)
go_ids <- c("GO:0034976", "GO:0097190", "GO:0033365", "GO:0042391", "GO:0099536")

gse@result %>%
  filter(ID %in% go_ids) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
  ylab("Term") +
  ggpubr::theme_pubr() +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(size = 16),
    axis.text.y  = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14, angle = 45, hjust = 1, vjust = 1)
  )

## S5H – Lineage1_2 GSEA
pseudotime_start_end_association_filt <- read_tsv(
  "lin1_2_pseudotime_start_end_association_filt.tsv"
)
genelist <- setNames(
  pseudotime_start_end_association_filt$logFClineage1,
  pseudotime_start_end_association_filt$feature_id
)
genelist <- sort(genelist, decreasing = TRUE)
genelist <- genelist[abs(genelist) > 1]

set.seed(123)
gse1_2 <- gseGO(
  geneList      = genelist,
  ont           = "BP",
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  eps           = 0,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)
gse1_2 <- simplify(gse1_2)
go_ids <- c("GO:0019932", "GO:0030001", "GO:0001578", "GO:0006882")

gse1_2@result %>%
  filter(ID %in% go_ids) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
  ylab("Term") +
  ggpubr::theme_pubr() +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(size = 16),
    axis.text.y  = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14, angle = 45, hjust = 1, vjust = 1)
  )

################################################################################
### 11. SCpubr DotPlots (Figure S5I)
################################################################################

genes <- c("INS", "CACNA1C", "CACNA1D", "ABCC8", "PCSK1", "PDX1", "GNAS")

p1 <- SCpubr::do_DotPlot(
  lin1_1,
  features          = genes,
  viridis.palette   = "F",
  viridis.direction = 1,
  use_viridis       = TRUE,
  flip              = FALSE,
  legend.position   = "right",
  legend.length     = 5,
  plot.title        = "Lineage1_1",
  group.by          = "monocle_clusters2",
  scale             = TRUE
) +
  scale_y_discrete(limits = c("2", "1", "3")) +
  ylab("Cluster") +
  NoLegend() +
  theme(
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text    = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5)
  )

p2 <- SCpubr::do_DotPlot(
  lin1_2,
  features          = genes,
  viridis.palette   = "F",
  viridis.direction = 1,
  use_viridis       = TRUE,
  flip              = FALSE,
  legend.position   = "bottom",
  legend.width      = 1,
  legend.length     = 7,
  plot.title        = "Lineage1_2",
  group.by          = "monocle_clusters2",
  scale             = TRUE
) +
  scale_y_discrete(limits = c("2", "1", "4")) +
  ylab("Cluster") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    plot.title   = element_text(hjust = 0.5),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text    = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14, hjust = 1, vjust = 1)
  )

print(p1 / p2)

################################################################################
### 12. Unified Pseudotime Trajectory Plots (Figures S5J–S5L)
################################################################################

# Extract expression data and pseudotime
exp <- data.integrated@assays$RNA@data %>% t()
gc()

sce <- readRDS("sce_lin1.rds")

pseudotime <- as.data.frame(colData(sce)$crv)
pseudotime$UnifiedPseudotime <- NA_real_

for (jj in 1:2) {
  lineage_col <- paste0("cellWeights.Lineage", jj)
  pseud_col   <- paste0("pseudotime.Lineage", jj)
  mask        <- pseudotime[[lineage_col]] == 1
  pseudotime$UnifiedPseudotime[mask] <- pseudotime[[pseud_col]][mask]
}

pseudotime$AdjustedPseudotime <- pseudotime$UnifiedPseudotime

for (lineage in 1:2) {
  mask   <- pseudotime[[paste0("cellWeights.Lineage", lineage)]] == 1
  min_pt <- min(pseudotime$UnifiedPseudotime[mask], na.rm = TRUE)
  max_pt <- max(pseudotime$UnifiedPseudotime[mask], na.rm = TRUE)
  pseudotime$AdjustedPseudotime[mask] <-
    (pseudotime$UnifiedPseudotime[mask] - min_pt) / (max_pt - min_pt)
}

pseudotime$Membership <- apply(
  pseudotime[, c("cellWeights.Lineage1", "cellWeights.Lineage2")],
  1,
  function(w) {
    lineages <- c("Lineage1_1", "Lineage1_2")
    paste0(lineages[w > 0], collapse = ", ")
  }
)

pseudotime$DominantLineage <- apply(
  pseudotime[, c("cellWeights.Lineage1", "cellWeights.Lineage2")],
  1,
  function(w) {
    lineages <- c("Lineage1_1", "Lineage1_2")
    if (all(w == 0)) return(NA_character_)
    lineages[which.max(w)]
  }
)

point_colors    <- c("Lineage1_1" = "#325A9BFF", "Lineage1_2" = "#2ED9FFFF")
smoother_colors <- point_colors

################## Figure S5J
# genes <- c("NDUFA1","UQCRC1","ATP5F1C","MAFA","ABCC8","PCSK1","CACNA1C","CACNA1D")
# l     <- c(3,3,3,3,3,3,1.5,1.5)

################## Figure S5K
# genes <- c("HSPA5","HSP90B1","DDIT3","SAR1B","CDKN1A","GADD45A","ATM","ATR")
# l     <- c(3,3,2,2,3,3,2,2)

################## Figure S5L
genes <- c("RPL3", "RPL13", "EEF1A1", "EIF3E", "EIF4A2", "RPLP2", "EEF2")
l     <- c(6, 5, 6, 3, 4, 4, 4)

plt_list <- list()
i <- 1
for (gene in genes) {
  plt_data <- data.frame(
    AdjustedPseudotime = pseudotime$AdjustedPseudotime,
    Membership         = pseudotime$Membership,
    DominantLineage    = pseudotime$DominantLineage,
    expr               = exp[, gene]
  )

  p <- ggplot(plt_data, aes(x = AdjustedPseudotime, y = expr)) +
    geom_point(aes(color = DominantLineage), size = 0.6, alpha = 0.1) +
    geom_smooth(
      data = subset(plt_data, grepl("Lineage1_1", Membership)),
      aes(group = "Lineage1_1"),
      size = 2, color = smoother_colors["Lineage1_1"],
      se = FALSE, method = "lm", formula = y ~ splines::ns(x, df = 2)
    ) +
    geom_smooth(
      data = subset(plt_data, grepl("Lineage1_2", Membership)),
      aes(group = "Lineage1_2"),
      size = 2, color = smoother_colors["Lineage1_2"],
      se = FALSE, method = "lm", formula = y ~ splines::ns(x, df = 2)
    ) +
    scale_color_manual(name = "Lineages", values = point_colors) +
    labs(x = "Pseudotime", y = gene) +
    theme_pubr(legend = "right") +
    theme(
      axis.title.y = element_text(face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text    = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(face = "bold", size = 14)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.5))) +
    ylim(0, l[i])

  plt_list[[gene]] <- p
  i <- i + 1
}
print(wrap_plots(plt_list, guides = "collect", nrow = 2))

################################################################################
### 13. Module Detection & GO Functions (Figures 3N–3Q)
################################################################################

mat1_z <- readRDS("mat1_z.rds")
mat2_z <- readRDS("mat2_z.rds")

pt1 <- coldata[colnames(mat1_z), "Lineage1_1"]
pt2 <- coldata[colnames(mat2_z), "Lineage1_2"]

# Get Slingshot curve weights
wts <- slingCurveWeights(sim)  # sds is your SlingshotDataSet object

# Cells that are shared between both lineages (i.e., near bifurcation)
shared_cells <- which(wts[, "Lineage1"] > 0.1 & wts[, "Lineage2"] > 0.1)

# Estimate bifurcation pseudotime as the max pseudotime in shared cells
pt_shared      <- slingPseudotime(sim)[shared_cells, ]
bifurcation_pt1 <- max(pt_shared[, "Lineage1"], na.rm = TRUE)
bifurcation_pt2 <- max(pt_shared[, "Lineage2"], na.rm = TRUE)

shared_coords <- reducedDim(sim)[shared_cells, ]

pre_branch_cells <- which(pt1 <= bifurcation_pt1 & pt2 <= bifurcation_pt2)
post_branch_l1   <- which(pt1 > bifurcation_pt1)
post_branch_l2   <- which(pt2 > bifurcation_pt2)

w1 <- pseudotime[, "cellWeights.Lineage1"]
w2 <- pseudotime[, "cellWeights.Lineage2"]

# Get pseudotime
pt1 <- slingPseudotime(sim)[, "Lineage1"]
pt2 <- slingPseudotime(sim)[, "Lineage2"]

# Estimate bifurcation point: cells with similar weights
shared_cells  <- which(abs(w1 - w2) < 0.2 & !is.na(pt1) & !is.na(pt2))
bifurcation_pt <- max(pmin(pt1[shared_cells], pt2[shared_cells]))

cat("Estimated bifurcation pseudotime:", round(bifurcation_pt, 2), "\n")
### Estimated bifurcation pseudotime: 12.1

# Minimum pseudotime across the two lineages
pt_min <- pmin(pt1, pt2, na.rm = TRUE)

# Label as Pre or Post based on bifurcation point
bifurcation_pt <- round(bifurcation_pt, 2)
branch_state <- ifelse(pt_min <= bifurcation_pt, "Pre-branch", "Post-branch")

# Add to colData
colData(cds)$Branch_state <- branch_state

#### Figure 3N
plot_cells(
  cds, show_trajectory_graph = FALSE,
  color_cells_by = "Branch_state",
  label_cell_groups = FALSE, cell_size = 1,
  group_cells_by = "cluster",
  group_label_size = 10
) +
  ggpubr::clean_theme() +
  theme(
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14, angle = 0, hjust = 1, vjust = 1)
  )

################################################################################
################################################################################

# Quick function to make elbow plot for choosing optimal PCs
make_elbow_plot <- function(expression_matrix, pseudotime,
                            sd_threshold = 0.90, nbins = 100, max_pcs = 15) {

  # Prepare data (same as your main function)
  row_sds <- matrixStats::rowSds(expression_matrix, na.rm = TRUE)
  names(row_sds) <- rownames(expression_matrix)
  sd_thresh <- quantile(row_sds, sd_threshold, na.rm = TRUE)
  variable_genes <- names(row_sds)[row_sds > sd_thresh]

  # FIX: Only keep genes that actually exist in the matrix
  variable_genes <- intersect(variable_genes, rownames(expression_matrix))
  cat("Using", length(variable_genes), "variable genes (top", (1 - sd_threshold) * 100, "%)\n")

  mat_sub <- expression_matrix[variable_genes, ]

  # Bin cells by pseudotime
  bin_id   <- cut(pseudotime, breaks = nbins, labels = FALSE)
  bin_list <- split(seq_along(bin_id), bin_id)
  bin_list <- bin_list[lengths(bin_list) > 0]

  binned_mat <- sapply(bin_list, function(ix) {
    rowMeans(mat_sub[, ix, drop = FALSE], na.rm = TRUE)
  })

  # Run PCA
  pca_result <- prcomp(binned_mat, center = TRUE, scale. = TRUE)

  # Get variance explained
  variance_explained <- summary(pca_result)$importance[2, ]
  cumulative_var     <- summary(pca_result)$importance[3, ]

  # Create data frame
  pc_df <- data.frame(
    PC                 = 1:min(max_pcs, length(variance_explained)),
    VarianceExplained  = variance_explained[1:min(max_pcs, length(variance_explained))],
    CumulativeVariance = cumulative_var[1:min(max_pcs, length(variance_explained))]
  )

  # Create elbow plot
  p1 <- ggplot(pc_df, aes(x = PC, y = VarianceExplained)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    labs(
      title = "Scree Plot - Variance Explained by Each PC",
      x = "Principal Component", y = "Proportion of Variance Explained"
    ) +
    theme_minimal() +
    scale_x_continuous(breaks = pc_df$PC) +
    theme(plot.title = element_text(hjust = 0.5))

  p2 <- ggplot(pc_df, aes(x = PC, y = CumulativeVariance)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed", alpha = 0.6, color = "red") +
    labs(
      title = "Cumulative Variance Explained",
      x = "Principal Component", y = "Cumulative Proportion"
    ) +
    theme_minimal() +
    scale_x_continuous(breaks = pc_df$PC) +
    annotate("text", x = max(pc_df$PC) * 0.7, y = c(0.8, 0.9),
             label = c("80%", "90%"), vjust = -0.5, color = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  # Print recommendations
  pcs_80 <- which(cumulative_var >= 0.80)[1]
  pcs_90 <- which(cumulative_var >= 0.90)[1]

  cat("PCA Variance Summary:\n")
  cat("- For 80% variance:", pcs_80, "PCs\n")
  cat("- For 90% variance:", pcs_90, "PCs\n")

  # Fix the elbow calculation to avoid potential errors
  if (length(variance_explained) >= 3) {
    elbow_pc <- which.max(diff(diff(variance_explained[1:min(10, length(variance_explained))]))) + 1
    cat("- Elbow appears around PC", elbow_pc, "\n")
  }

  # Combine plots
  combined_plot <- p1 / p2

  list(
    variance_df         = pc_df,
    elbow_plot          = combined_plot,
    recommended_pcs_80  = pcs_80,
    recommended_pcs_90  = pcs_90
  )
}

# Main function to identify pseudotime modules using PCA approach
identify_pseudotime_modules_pca <- function(expression_matrix, pseudotime,
                                            lineage_name = "Lineage",
                                            nbins = 100,
                                            min_genes = 10,
                                            sd_threshold = 0.90,
                                            perc_expr = 5,
                                            n_pcs = 5) {

  cat("Processing", lineage_name, "with", ncol(expression_matrix), "cells and", nrow(expression_matrix), "genes\n")

  row_sds <- matrixStats::rowSds(expression_matrix, na.rm = TRUE)
  names(row_sds) <- rownames(expression_matrix)
  sd_thresh <- quantile(row_sds, sd_threshold, na.rm = TRUE)
  variable_genes <- names(row_sds)[row_sds > sd_thresh]
  variable_genes <- intersect(variable_genes, rownames(expression_matrix))

  cat("Using", length(variable_genes), "variable genes (top", (1 - sd_threshold) * 100, "%)\n")

  mat_sub <- expression_matrix[variable_genes, ]

  bin_id   <- cut(pseudotime, breaks = nbins, labels = FALSE)
  bin_list <- split(seq_along(bin_id), bin_id)
  bin_list <- bin_list[lengths(bin_list) > 0]
  cat("Created", length(bin_list), "non-empty pseudotime bins\n")

  binned_mat <- sapply(bin_list, function(ix) {
    rowMeans(mat_sub[, ix, drop = FALSE], na.rm = TRUE)
  })
  rownames(binned_mat) <- rownames(mat_sub)

  pca_result <- prcomp(binned_mat, center = TRUE, scale. = TRUE)
  coords     <- pca_result$x[, 1:min(n_pcs, ncol(pca_result$x))]

  cat(
    "PCA variance explained by first", ncol(coords), "PCs:",
    round(sum(summary(pca_result)$importance[2, 1:ncol(coords)]) * 100, 1), "%\n"
  )

  dist_mat <- dist(coords)
  hc <- hclust(dist_mat, method = "ward.D2")

  modules_raw <- dynamicTreeCut::cutreeDynamic(
    dendro = hc,
    distM  = as.matrix(dist_mat),
    deepSplit = 2,
    minClusterSize = min_genes
  )
  names(modules_raw) <- rownames(coords)

  module_list <- split(names(modules_raw), modules_raw)
  module_list <- module_list[names(module_list) != "0"]

  # Optional filtering by % expression
  if (file.exists("Lineage1_sling_perc_expr.tsv")) {
    perc_df <- readr::read_tsv("Lineage1_sling_perc_expr.tsv", show_col_types = FALSE)
    names(perc_df)[1:2] <- c("gene", "pct_expr")
    perc_cutoff <- perc_expr
    expr_pct_vec <- setNames(perc_df$pct_expr, perc_df$gene)
    module_list <- lapply(module_list, function(glist) {
      glist[glist %in% names(expr_pct_vec) & expr_pct_vec[glist] >= perc_cutoff]
    })
  }

  module_list <- Filter(function(g) length(g) >= min_genes, module_list)
  cat("Found", length(module_list), "modules with >=", min_genes, "genes\n")
  cat("Module sizes:", sapply(module_list, length), "\n")

  module_df <- data.frame(
    gene    = names(modules_raw),
    module  = modules_raw,
    lineage = lineage_name,
    stringsAsFactors = FALSE
  )

  mean_expr_modules <- sapply(module_list, function(genes) {
    genes <- intersect(genes, rownames(mat_sub))
    if (length(genes) > 0) {
      colMeans(mat_sub[genes, , drop = FALSE], na.rm = TRUE)
    } else {
      rep(NA, ncol(mat_sub))
    }
  })

  df_plot <- data.frame(
    pseudotime = pseudotime,
    mean_expr_modules
  )

  df_long <- tidyr::pivot_longer(
    df_plot, cols = -pseudotime, names_to = "Module", values_to = "Expression"
  )
  df_long$Module <- paste0(lineage_name, "_M", df_long$Module)

  list(
    modules           = module_list,
    module_df         = module_df,
    binned_matrix     = binned_mat,
    pca_result        = pca_result,
    plot_data         = df_long,
    mean_expr_modules = mean_expr_modules
  )
}

# Function for GO enrichment analysis
perform_go_enrichment <- function(module_list, lineage_name = "Lineage",
                                  org_db = org.Hs.eg.db, ont = "BP") {

  cat("Performing GO enrichment for", length(module_list), "modules in", lineage_name, "\n")

  ego_list <- list()

  for (i in seq_along(module_list)) {
    module_name <- names(module_list)[i]
    if (is.null(module_name)) module_name <- i

    cat("Processing module", module_name, "with", length(module_list[[i]]), "genes\n")

    ego <- enrichGO(
      gene          = module_list[[i]],
      OrgDb         = org_db,
      keyType       = "SYMBOL",
      ont           = ont,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500
    )

    if (!is.null(ego) && nrow(ego@result) > 0) {
      # ego_list[[module_name]] <- simplify(ego)
      ego_list[[module_name]] <- ego  # Removed simplify()
      cat("Found", nrow(ego@result), "enriched terms\n")
    } else {
      ego_list[[module_name]] <- NULL
      cat("No enriched terms found\n")
    }
  }

  ego_list
}

# Assuming `ego_pre` is your list of enrichGO results per module
get_module_names_from_go <- function(ego_list) {
  module_names <- sapply(names(ego_list), function(mod_id) {
    ego <- ego_list[[mod_id]]
    if (!is.null(ego) && nrow(ego@result) > 0) {
      top_term <- ego@result %>%
        dplyr::filter(p.adjust < 0.05) %>%
        dplyr::slice_min(p.adjust, n = 1) %>%
        dplyr::pull(Description)
      if (length(top_term) > 0) return(top_term)
    }
    return(NA)
  })
  names(module_names) <- paste0("Module_", names(module_names))
  module_names
}

# Function to plot module trajectories
plot_module_trajectories <- function(results, title = NULL, smooth_method = "loess") {
  if (is.null(title)) title <- "Module Expression Trajectories"

  ggplot(results$plot_data, aes(x = pseudotime, y = Expression, color = Module)) +
    geom_smooth(method = smooth_method, se = FALSE, size = 1.2) +
    labs(title = title, x = "Pseudotime", y = "Average Expression (log-normalized)") +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title   = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))
}

###################################################

# Combine scaled expression matrices used in the module detection
shared_cells <- intersect(colnames(mat1_z), colnames(mat2_z))
mat1_unique  <- mat1_z[, setdiff(colnames(mat1_z), shared_cells)]
mat2_unique  <- mat2_z[, setdiff(colnames(mat2_z), shared_cells)]
shared_expr  <- mat1_z[, shared_cells]
combined_expr_matrix <- cbind(mat1_unique, mat2_unique, shared_expr)
dim(combined_expr_matrix)

# Shared aging trajectory
cells_pre <- which(pt1 <= bifurcation_pt & pt2 <= bifurcation_pt)

# Sublineage-specific post-bifurcation
cells_post_l1 <- which(pt1 > bifurcation_pt)
cells_post_l2 <- which(pt2 > bifurcation_pt)

# Filter expression matrices
mat_pre <- combined_expr_matrix[, cells_pre]
mat_l1  <- combined_expr_matrix[, cells_post_l1]
mat_l2  <- combined_expr_matrix[, cells_post_l2]

pt_pre <- pt1[cells_pre]
pt_l1  <- pt1[cells_post_l1]
pt_l2  <- pt2[cells_post_l2]

elbow_results_pt_pre <- make_elbow_plot(mat_pre, pt_pre)
elbow_results_pt_pre$elbow_plot

results_pre <- identify_pseudotime_modules_pca(
  expression_matrix = mat_pre,
  pseudotime = pt_pre,
  lineage_name = "PreBifurcation",
  nbins = 25,
  min_genes = 100,
  perc_expr = 10,
  sd_threshold = 0.75,
  n_pcs = 3
)

ego_pre <- perform_go_enrichment(results_pre$modules, "pre")
module_labels <- get_module_names_from_go(ego_pre)
print(module_labels)

################################################################################

manuscript_colors <- list(
  pre_colors = c(
    "Aerobic Respiration (1)"                 = "#E31A1C",
    "Aerobic Respiration (2)"                 = "#FF7F00",
    "Translation + Vesicle Trafficking"       = "#1F78B4",
    "Insulin/Glucose Secretion"               = "#33A02C",
    "Translation + Hormone & Stress Response" = "#6A3D9A"
  ),
  l1_colors = c(
    "Ubiquitin–Proteasome Pathway"       = "#E31A1C",
    "RNA splicing + mRNA processing"     = "#FF7F00",
    "Unannotated"                        = "#999999",
    "Cytoplasmic Translation"            = "#6A3D9A"
  ),
  l2_colors = c(
    "Ubiquitin–Proteasome Pathway"       = "#E31A1C",
    "RNA splicing"                       = "#FF7F00",
    "Protein Degradation (ERAD)"         = "#1F78B4",
    "ER Stress & Vesicle Transport"      = "#B15928"
  )
)

################################################################################

library(geomtextpath)

plot_pre <- plot_module_trajectories(results_pre, "Pre-branching")
p_all <- plot_pre +
  ggpubr::theme_pubr(legend = "right") +
  scale_color_manual(
    labels = c(
      "PreBifurcation_MX1" = "Aerobic Respiration (1)",
      "PreBifurcation_MX2" = "Aerobic Respiration (2)",
      "PreBifurcation_MX3" = "Translation + Vesicle Trafficking",
      "PreBifurcation_MX4" = "Insulin/Glucose Secretion",
      "PreBifurcation_MX5" = "Translation + Hormone & Stress \nResponse"
    ),
    values = c(
      "PreBifurcation_MX1" = "#E31A1C",
      "PreBifurcation_MX2" = "#FF7F00",
      "PreBifurcation_MX3" = "#1F78B4",
      "PreBifurcation_MX4" = "#33A02C",
      "PreBifurcation_MX5" = "#6A3D9A"
    )
  ) +
  ylab("Average Expression") +
  geom_vline(xintercept = 12.1, linetype = "dotted", color = "blue", size = 1.5) +
  geom_textvline(label = "Pre-branching", xintercept = 12.1, vjust = -0.2, size = 5, fontface = "bold") +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14)
  )

p_all

#################################################################################

results_l1 <- identify_pseudotime_modules_pca(
  expression_matrix = mat_l1,
  pseudotime = pt_l1,
  lineage_name = "Lineage1_Post",
  nbins = 25,
  min_genes = 100,
  perc_expr = 10,
  sd_threshold = 0.75,
  n_pcs = 3
)

ego_lin1_1 <- perform_go_enrichment(results_l1$modules, "post_l1_1")
plot_go_enrichment(ego_lin1_1, "3")

############
module_labels <- get_module_names_from_go(ego_lin1_1)
print(module_labels)

##########################

plot_lin1_1 <- plot_module_trajectories(results_l1, "lineage1_1")
p1_1 <- plot_lin1_1 +
  ylim(-1, 1) +
  ggpubr::theme_pubr(legend = "right") +
  scale_color_manual(
    labels = c(
      "Lineage1_Post_MX1" = "Ubiquitin–Proteasome Pathway",
      "Lineage1_Post_MX2" = "RNA splicing + mRNA processing",
      "Lineage1_Post_MX3" = "Epigenetic–proteostatic control*",
      "Lineage1_Post_MX4" = "Cytoplasmic Translation"
    ),
    values = c(
      "Lineage1_Post_MX1" = "#E31A1C",
      "Lineage1_Post_MX2" = "#FF7F00",
      "Lineage1_Post_MX3" = "#999999",
      "Lineage1_Post_MX4" = "#6A3D9A"
    )
  ) +
  xlim(11.7, 18) +
  geom_vline(xintercept = 12.1, linetype = "dotted", color = "blue", size = 1.5) +
  geom_textvline(label = "Post-branching", xintercept = 12.1, vjust = -0.3, size = 5, fontface = "bold") +
  ylab(NULL) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14)
  )
p1_1

################################################################################
################################################################################

results_l2 <- identify_pseudotime_modules_pca(
  expression_matrix = mat_l2,
  pseudotime = pt_l2,
  lineage_name = "Lineage2_Post",
  nbins = 25,
  min_genes = 100,
  perc_expr = 10,
  sd_threshold = 0.75,
  n_pcs = 3
)

ego_lin1_2 <- perform_go_enrichment(results_l2$modules, "post_l1_2")
# plot_go_enrichment(ego_lin1_2, "3")

################################################################################
module_labels <- get_module_names_from_go(ego_lin1_2)
print(module_labels)

################################################################################

plot_lin1_2 <- plot_module_trajectories(results_l2, "lineage1_2")
p1_2 <- plot_lin1_2 +
  ylim(-1, 1) +
  ggpubr::theme_pubr(legend = "right") +
  scale_color_manual(
    labels = c(
      "Lineage2_Post_MX1" = "Ubiquitin–Proteasome Pathway",
      "Lineage2_Post_MX2" = "RNA splicing",
      "Lineage2_Post_MX3" = "Protein Degradation (ERAD)",
      "Lineage2_Post_MX4" = "ER Stress & Vesicle Transport"
    ),
    values = c(
      "Lineage2_Post_MX1" = "#E31A1C",
      "Lineage2_Post_MX2" = "#FF7F00",
      "Lineage2_Post_MX3" = "#1F78B4",
      "Lineage2_Post_MX4" = "#B15928"
    )
  ) +
  xlim(11.8, 16) +
  geom_vline(xintercept = 12.1, linetype = "dotted", color = "blue", size = 1.5) +
  geom_textvline(label = "Post-branching", xintercept = 12.1, vjust = -0.3, size = 5, fontface = "bold") +
  ylab(NULL) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 14)
  )

p1_2

################################################################################

# Cairo::CairoPNG("module_GO_terms.png", width = 23 * 300, height = 3.5 * 300, res = 300)
p_all | p1_1 | p1_2
dev.off()

################################################################################
# sessionInfo()
################################################################################

