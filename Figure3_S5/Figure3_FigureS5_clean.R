###############################################################################
# Project: scRNA beta cell atlas — Lineage 1 reintegration and Fig 3 suite
# File: Fig3_FigS5.R
# Author: Elhadi Iich
# Date: 2026-09-07
#
# Description:
#   End-to-end workflow for lineage 1 reintegration and visualization:
#   reclustering, supplementary Figs, Slingshot overlays, INS trends,
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
### 3. Fig 3B–3D
################################################################################

## Fig 3B
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

## Fig 3C
C3 <- df %>%
  drop_na(Age) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  ggplot(aes(x = monocle_clusters2, y = Age, fill = monocle_clusters2)) +
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
print(C3)

## Fig 3D
D3 <- df %>%
  drop_na(BMI) %>%
  mutate(monocle_clusters = as.character(monocle_clusters2)) %>%
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
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  )
print(D3)

## Combine 3C + 3D
print(C3 / D3)

################################################################################
### 4. Supplementary Fig S5A–S5D
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

################################################################################
################################################################################
# Reconstruct harmonized donor IDs for donor-level HbA1c analysis
################################################################################

annotation_dir <- "../metadata"

ann81608 <- readr::read_csv(
  file.path(annotation_dir, "GSE81608.csv"),
  show_col_types = FALSE
)

ann153855 <- readr::read_csv(
  file.path(annotation_dir, "GSE153855.csv"),
  show_col_types = FALSE
)

ann81608 <- ann81608 %>%
  transmute(
    sample_name = as.character(Run),
    annotated_donor = as.character(donor_id)
  ) %>%
  filter(
    !is.na(sample_name),
    sample_name != "",
    !is.na(annotated_donor),
    annotated_donor != ""
  ) %>%
  distinct()

ann153855 <- ann153855 %>%
  transmute(
    sample_name = as.character(Run),
    annotated_biosample = as.character(BioSample)
  ) %>%
  filter(
    !is.na(sample_name),
    sample_name != "",
    !is.na(annotated_biosample),
    annotated_biosample != ""
  ) %>%
  distinct()

fallback_studies <- c(
  "GSE98887",
  "GSE124742",
  "EMTAB5061",
  "GSE86469",
  "GSE81547"
)

signature_value <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x
}

df <- df %>%
  left_join(ann81608, by = "sample_name") %>%
  left_join(ann153855, by = "sample_name") %>%
  mutate(
    donor_signature = paste(
      signature_value(Age),
      signature_value(BMI),
      signature_value(Hba1c),
      signature_value(Gender),
      signature_value(Ethnicity),
      signature_value(Condition),
      sep = "|"
    ),
    donor_id = case_when(
      orig.ident == "GSE81608" ~
        paste0("GSE81608::", annotated_donor),
      
      orig.ident == "GSE153855" ~
        paste0("GSE153855::", annotated_biosample),
      
      Study %in% fallback_studies ~
        paste0(Study, "::metadata::", donor_signature),
      
      TRUE ~
        paste0(Study, "::", orig.ident)
    )
  )

## S5E
hba1c_donor_cluster <- df %>%
  filter(
    !is.na(Hba1c),
    !is.na(monocle_clusters2)
  ) %>%
  group_by(donor_id, monocle_clusters2) %>%
  summarise(
    Hba1c = first(Hba1c),
    .groups = "drop"
  )

S5E <- ggplot(
  hba1c_donor_cluster,
  aes(
    x = factor(monocle_clusters2, levels = desired_order),
    y = Hba1c,
    fill = factor(monocle_clusters2, levels = desired_order)
  )
) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(width = 0.15, size = 1.2) +
  geom_hline(yintercept = 5.7, linetype = "dashed") +
  ggpubr::theme_pubr(legend = "none") +
  labs(
    x = "Cluster",
    y = "Donor HbA1c (%)"
  )+
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text  = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14)
  )
print(S5E)


################################################################################
### 5. Slingshot Curves & UMAP Labels
################################################################################

sim <- readRDS("slingshot_sim_lin1.rds")

coldata <- data.frame(
  study        = sim$study,
  Cluster      = sim$monocle_clusters2,
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
### 6. Slingshot Branch Plots & INS Expression (Figs 3G–3H)
################################################################################

lin1_1 <- readRDS(
  "Beta_T2D_lineage1_1.rds"
)
dim(lin1_1)

lin1_2 <- readRDS(
  "Beta_T2D_lineage1_2.rds"
)
dim(lin1_2)

###### Fig 3G
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

###### Fig 3H
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
### 7. Heatmap Data Prep, Drawing & Legends (Fig 3I)
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
  Cluster_legend,
  Lineage_legend,
  Pseudotime_legend,
  Expression_legend,
  direction = "vertical"
)
draw(all_legends, just = "left")
dev.off()

################################################################################
### 8. Pseudotime Association Volcano Plots (Figs 3J–3K)
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

### Fig 3J
volcano_1_1 <- create_volcano(pseudotime_start_end_association1_1, "logFClineage1")
print(volcano_1_1)

### Fig 3K
sce <- readRDS("sce_lin1_2.rds")
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
### 9. Pseudobulk DESeq2 Analysis (Figs 3L–3M)
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
  (res_filt_2$log2FoldChange > 1 & res_filt_2$padj <= 0.05), "UP",
  ifelse((res_filt_2$log2FoldChange < -1 & res_filt_2$padj <= 0.05), "DOWN", "NS")
)

genes_l <- c(
  "HSPA5", "HERPUD1", "HSP90B1", "BHLHE41", "DDIT3", "SLC30A8",
  "NKX6-1", "PAX6", "MAFA", "INS", "IAPP", "PDX1", "G6PC2"
)

############### Fig 3L
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
### 10. GO Enrichment Dotplots (Figs 3M & S5F–S5G)
################################################################################

genes_C3 <- res_filt_2 %>%
  filter( padj < 0.05, log2FoldChange > 1 ) %>%
  pull(genes)

set.seed(123)

ego <- enrichGO(
  gene          = genes_C3,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500
) 


dotplot(ego) +
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

## S5F – Lineage1_1 GSEA
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

## S5G – Lineage1_2 GSEA
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
### 11. SCpubr DotPlots (Fig S5H)
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
### 12. Unified Pseudotime Trajectory Plots (Figs S5I–S5J)
################################################################################

# Extract expression data and pseudotime
exp <- data.integrated@assays$RNA@data %>% t()
gc()

sce <- readRDS("sce_lin1.rds")

crv <- as.data.frame(colData(sce)$crv)

pt <- as.matrix(
  crv[, c(
    "pseudotime.Lineage1",
    "pseudotime.Lineage2"
  )]
)

wt <- as.matrix(
  crv[, c(
    "cellWeights.Lineage1",
    "cellWeights.Lineage2"
  )]
)

rownames(pt) <- colnames(sce)
rownames(wt) <- colnames(sce)

colnames(pt) <- c("Lineage1_1", "Lineage1_2")
colnames(wt) <- c("Lineage1_1", "Lineage1_2")

point_colors <- c(
  "Lineage1_1" = "#325A9BFF",
  "Lineage1_2" = "#2ED9FFFF"
)

smoother_colors <- point_colors


################## Fig S5I
# genes <- c("NDUFA1","UQCRC1","ATP5F1C","MAFA","ABCC8","PCSK1","CACNA1C","CACNA1D")
# l     <- c(3,3,3,3,3,3,1.5,1.5)

################## Fig S5J
################## Fig S5J — genes increased along Lineage1_2

################## Fig S5J — preferentially increased along Lineage1_2

genes <- c(
  "EIF4A2",   # translation initiation
  "SEC61B",   # ER protein translocation
  "HSPA8",    # protein folding / proteostasis
  "SLC30A8",  # beta-cell insulin granule function
  "G6PC2",    # beta-cell glucose sensing
  "NEUROD1"   # beta-cell/endocrine identity
)

# l     <- c(3,3,2,2,3,3,2,2)

################## Fig S5J — declining genes along Lineage1_2



plt_list <- list()

for (i in seq_along(genes)) {
  
  gene <- genes[i]
  
  gene_expr <- exp[, gene]
  
  plot_data <- lapply(
    seq_len(ncol(pt)),
    function(j) {
      
      keep <- !is.na(pt[, j]) & wt[, j] > 0
      
      x <- pt[keep, j]
      
      tibble(
        cell = rownames(pt)[keep],
        lineage = colnames(pt)[j],
        pseudotime = x,
        weight = wt[keep, j],
        expr = as.numeric(gene_expr[rownames(pt)[keep]])
      )
    }
  ) %>%
    bind_rows()
  
  p <- ggplot(
    plot_data,
    aes(x = pseudotime, y = expr)
  ) +
    geom_point(
      aes(color = lineage),
      size = 0.6,
      alpha = 0.05
    ) +
    geom_smooth(
      data = subset(plot_data, lineage == "Lineage1_1"),
      aes(
        group = lineage,
        weight = weight
      ),
      size = 2,
      color = smoother_colors["Lineage1_1"],
      se = FALSE,
      method = "lm",
      formula = y ~ splines::ns(x, df = 2)
    ) +
    geom_smooth(
      data = subset(plot_data, lineage == "Lineage1_2"),
      aes(
        group = lineage,
        weight = weight
      ),
      size = 2,
      color = smoother_colors["Lineage1_2"],
      se = FALSE,
      method = "lm",
      formula = y ~ splines::ns(x, df = 2)
    ) +
    scale_color_manual(
      name = "Lineages",
      values = point_colors
    ) +
    labs(
      x = "Pseudotime",
      y = gene
    ) +
    theme_pubr(legend = "right") +
    theme(
      axis.title.y = element_text(face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 14)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 3, alpha = 0.5)
      )
    ) #+
    # coord_cartesian(
    #   ylim = c(0, l[i])
    # )
  
  plt_list[[gene]] <- p
}

print(
  wrap_plots(
    plt_list,
    guides = "collect",
    nrow = 2
  )
)


################################################################################
### 13. Module Detection & GO Functions (Figs 3N–3Q)
################################################################################

mat1_z <- readRDS("mat1_z.rds")
mat2_z <- readRDS("mat2_z.rds")

################################################################################
# Build pre- and post-branch matrices using CELL NAMES
################################################################################

shared_expr_cells <- intersect(
  colnames(mat1_z),
  colnames(mat2_z)
)

mat1_unique <- mat1_z[
  ,
  setdiff(colnames(mat1_z), shared_expr_cells),
  drop = FALSE
]

mat2_unique <- mat2_z[
  ,
  setdiff(colnames(mat2_z), shared_expr_cells),
  drop = FALSE
]

shared_expr <- mat1_z[
  ,
  shared_expr_cells,
  drop = FALSE
]

combined_expr_matrix <- cbind(
  mat1_unique,
  mat2_unique,
  shared_expr
)


pt_sim <- slingPseudotime(sim)
wt_sim <- slingCurveWeights(sim)

pt1_all <- pt_sim[, "Lineage1"]
pt2_all <- pt_sim[, "Lineage2"]

w1 <- wt_sim[, "Lineage1"]
w2 <- wt_sim[, "Lineage2"]

names(pt1_all) <- rownames(pt_sim)
names(pt2_all) <- rownames(pt_sim)
names(w1) <- rownames(wt_sim)
names(w2) <- rownames(wt_sim)

shared_idx <- which(
  abs(w1 - w2) < 0.2 &
    !is.na(pt1_all) &
    !is.na(pt2_all)
)

bifurcation_pt <- max(
  pmin(
    pt1_all[shared_idx],
    pt2_all[shared_idx]
  ),
  na.rm = TRUE
)

cat("Bifurcation:", bifurcation_pt, "\n")

# Pre-branch cells
pre_cells <- intersect(
  colnames(combined_expr_matrix),
  names(pt1_all)[
    !is.na(pt1_all) &
      !is.na(pt2_all) &
      pt1_all <= bifurcation_pt &
      pt2_all <= bifurcation_pt
  ]
)

# Lineage1_1 post-branch
post_l1_cells <- intersect(
  colnames(mat1_z),
  names(pt1_all)[
    !is.na(pt1_all) &
      pt1_all > bifurcation_pt
  ]
)

# Lineage1_2 post-branch
post_l2_cells <- intersect(
  colnames(mat2_z),
  names(pt2_all)[
    !is.na(pt2_all) &
      pt2_all > bifurcation_pt
  ]
)

mat_pre <- combined_expr_matrix[
  ,
  pre_cells,
  drop = FALSE
]

mat_l1 <- mat1_z[
  ,
  post_l1_cells,
  drop = FALSE
]

mat_l2 <- mat2_z[
  ,
  post_l2_cells,
  drop = FALSE
]

pt_pre <- pmin(
  pt1_all[pre_cells],
  pt2_all[pre_cells]
)

pt_l1 <- pt1_all[post_l1_cells]
pt_l2 <- pt2_all[post_l2_cells]

cat(
  "Pre-branch:", ncol(mat_pre), "\n",
  "Lineage1_1 post-branch:", ncol(mat_l1), "\n",
  "Lineage1_2 post-branch:", ncol(mat_l2), "\n"
)

################################################################################
# Branch-state annotation for Fig 3N
################################################################################

pt_min <- pmin(
  pt1_all,
  pt2_all,
  na.rm = TRUE
)

branch_state <- ifelse(
  pt_min <= bifurcation_pt,
  "Pre-branch",
  "Post-branch"
)

names(branch_state) <- names(pt1_all)

colData(cds)$Branch_state <- branch_state[colnames(cds)]


#### Fig 3N
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

###############################################################################
# SOURCE DATA EXPORT — Fig 3 and Fig S5
###############################################################################

source_dir <- "../SourceData/Fig3_S5"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# Fig 3B, Fig S5A, S5C and Fig 3N
# UMAP coordinates + plotted annotations
###############################################################################

umap_coords_source <- as.data.frame(
  reducedDims(cds)$UMAP
)

colnames(umap_coords_source)[1:2] <- c("UMAP_1", "UMAP_2")

umap_coords_source <- umap_coords_source %>%
  tibble::rownames_to_column("cell_id")

cds_meta_source <- as.data.frame(
  SummarizedExperiment::colData(cds)
) %>%
  tibble::rownames_to_column("cell_id")

# Monocle's own cluster assignment used in the Slingshot overlay
monocle3_cluster_source <- data.frame(
  cell_id = colnames(cds),
  monocle3_cluster = as.character(
    monocle3::clusters(cds, reduction_method = "UMAP")
  )
)

fig3_umap_source <- umap_coords_source %>%
  dplyr::left_join(
    cds_meta_source %>%
      dplyr::select(
        cell_id,
        dplyr::any_of(c(
          "monocle_clusters2",
          "Study",
          "Condition",
          "Branch_state"
        ))
      ),
    by = "cell_id"
  ) %>%
  dplyr::left_join(
    monocle3_cluster_source,
    by = "cell_id"
  )

readr::write_tsv(
  fig3_umap_source,
  file.path(
    source_dir,
    "Fig3B_S5A_S5C_3N_UMAP.tsv.gz"
  )
)


###############################################################################
# Slingshot curves plotted over Fig 3 UMAP
###############################################################################

fig3_slingshot_source <- combined_lineages %>%
  dplyr::select(
    lineage,
    path_id,
    Dim_1,
    Dim_2
  )

readr::write_tsv(
  fig3_slingshot_source,
  file.path(
    source_dir,
    "Fig3_slingshot_curves.tsv"
  )
)


###############################################################################
# Figs 3C and 3D
# Cell-weighted Age and BMI distributions
###############################################################################

fig3cd_source <- df %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    monocle_clusters2,
    Age,
    BMI
  )

readr::write_tsv(
  fig3cd_source,
  file.path(
    source_dir,
    "Fig3C_D_age_BMI.tsv.gz"
  )
)


###############################################################################
# Figs S5B and S5D
# Study and disease-condition composition
###############################################################################

s5b_source <- df %>%
  dplyr::count(
    monocle_clusters2,
    Study,
    name = "cell_count"
  ) %>%
  dplyr::group_by(monocle_clusters2) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count),
    variable = "Study",
    category = as.character(Study)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = monocle_clusters2,
    variable,
    category,
    cell_count,
    proportion
  )

s5d_source <- df %>%
  dplyr::count(
    monocle_clusters2,
    Condition,
    name = "cell_count"
  ) %>%
  dplyr::group_by(monocle_clusters2) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count),
    variable = "Condition",
    category = as.character(Condition)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = monocle_clusters2,
    variable,
    category,
    cell_count,
    proportion
  )

readr::write_tsv(
  dplyr::bind_rows(s5b_source, s5d_source),
  file.path(
    source_dir,
    "FigS5B_D_composition.tsv"
  )
)


###############################################################################
# Fig S5E
# Donor-level HbA1c values
###############################################################################

readr::write_tsv(
  hba1c_donor_cluster,
  file.path(
    source_dir,
    "FigS5E_HbA1c.tsv"
  )
)


###############################################################################
# Figs 3G and 3H
# INS expression versus Slingshot pseudotime
###############################################################################

fig3g_source <- data.frame(
  cell_id = colnames(lin1_1),
  lineage = "Lineage1_1",
  pseudotime = lin1_1@meta.data[
    colnames(lin1_1),
    "Lineage1_1"
  ],
  cluster = lin1_1@meta.data[
    colnames(lin1_1),
    "monocle_clusters2"
  ],
  INS = as.numeric(
    lin1_1@assays$RNA@data[
      "INS",
      colnames(lin1_1)
    ]
  )
)

fig3h_source <- data.frame(
  cell_id = colnames(lin1_2),
  lineage = "Lineage1_2",
  pseudotime = lin1_2@meta.data[
    colnames(lin1_2),
    "Lineage1_2"
  ],
  cluster = lin1_2@meta.data[
    colnames(lin1_2),
    "monocle_clusters2"
  ],
  INS = as.numeric(
    lin1_2@assays$RNA@data[
      "INS",
      colnames(lin1_2)
    ]
  )
)

readr::write_tsv(
  dplyr::bind_rows(
    fig3g_source,
    fig3h_source
  ),
  file.path(
    source_dir,
    "Fig3G_H_INS_pseudotime.tsv.gz"
  )
)




###############################################################################
# Figs 3J and 3K
# tradeSeq start-vs-end statistics
###############################################################################

readr::write_tsv(
  pseudotime_start_end_association1_1,
  file.path(
    source_dir,
    "Fig3J_pseudotime_DE.tsv"
  )
)

readr::write_tsv(
  pseudotime_start_end_association1_2,
  file.path(
    source_dir,
    "Fig3K_pseudotime_DE.tsv"
  )
)


###############################################################################
# Fig 3L
# Pseudobulk DESeq2 results
###############################################################################

readr::write_tsv(
  res_filt_2,
  file.path(
    source_dir,
    "Fig3L_pseudobulk_DE.tsv"
  )
)


###############################################################################
# Fig 3M
# GO enrichment underlying dot plot
###############################################################################

fig3m_source <- as.data.frame(
  ego@result
)

readr::write_tsv(
  fig3m_source,
  file.path(
    source_dir,
    "Fig3M_GO_enrichment.tsv"
  )
)


###############################################################################
# Figs S5F and S5G
# GSEA results
###############################################################################

s5f_ids <- c(
  "GO:0034976",
  "GO:0097190",
  "GO:0033365",
  "GO:0042391",
  "GO:0099536"
)

s5g_ids <- c(
  "GO:0019932",
  "GO:0030001",
  "GO:0001578",
  "GO:0006882"
)

s5f_source <- as.data.frame(gse@result) %>%
  dplyr::filter(ID %in% s5f_ids) %>%
  dplyr::mutate(panel = "S5F_Lineage1_1")

s5g_source <- as.data.frame(gse1_2@result) %>%
  dplyr::filter(ID %in% s5g_ids) %>%
  dplyr::mutate(panel = "S5G_Lineage1_2")

readr::write_tsv(
  dplyr::bind_rows(
    s5f_source,
    s5g_source
  ),
  file.path(
    source_dir,
    "FigS5F_G_GSEA.tsv"
  )
)


###############################################################################
# Fig S5H
# Dot-plot values for selected beta-cell genes
###############################################################################

s5h_genes <- c(
  "INS",
  "CACNA1C",
  "CACNA1D",
  "ABCC8",
  "PCSK1",
  "PDX1",
  "GNAS"
)

s5h_l1 <- Seurat::DotPlot(
  lin1_1,
  features = s5h_genes,
  group.by = "monocle_clusters2"
)$data %>%
  as.data.frame() %>%
  dplyr::mutate(lineage = "Lineage1_1")

s5h_l2 <- Seurat::DotPlot(
  lin1_2,
  features = s5h_genes,
  group.by = "monocle_clusters2"
)$data %>%
  as.data.frame() %>%
  dplyr::mutate(lineage = "Lineage1_2")

readr::write_tsv(
  dplyr::bind_rows(
    s5h_l1,
    s5h_l2
  ),
  file.path(
    source_dir,
    "FigS5H_dotplot.tsv"
  )
)


###############################################################################
# Fig S5J — cell-level expression and pseudotime
###############################################################################

s5j_genes <- c(
  "EIF4A2",
  "SEC61B",
  "HSPA8",
  "SLC30A8",
  "G6PC2",
  "NEUROD1"
)

# Re-load trajectory object so this export does not depend on earlier objects
sce_s5j <- readRDS("sce_lin1.rds")

crv_s5j <- as.data.frame(
  SummarizedExperiment::colData(sce_s5j)$crv
)

pt_s5j <- as.matrix(
  crv_s5j[, c(
    "pseudotime.Lineage1",
    "pseudotime.Lineage2"
  ), drop = FALSE]
)

wt_s5j <- as.matrix(
  crv_s5j[, c(
    "cellWeights.Lineage1",
    "cellWeights.Lineage2"
  ), drop = FALSE]
)

rownames(pt_s5j) <- colnames(sce_s5j)
rownames(wt_s5j) <- colnames(sce_s5j)

colnames(pt_s5j) <- c(
  "Lineage1_1",
  "Lineage1_2"
)

colnames(wt_s5j) <- c(
  "Lineage1_1",
  "Lineage1_2"
)

# Pull only the six genes plotted in S5J
s5j_expr <- Seurat::GetAssayData(
  data.integrated,
  assay = "RNA",
  layer = "data"
)[s5j_genes, , drop = FALSE]

s5j_source <- lapply(
  s5j_genes,
  function(gene) {
    
    lapply(
      seq_len(ncol(pt_s5j)),
      function(j) {
        
        keep <- !is.na(pt_s5j[, j]) & wt_s5j[, j] > 0
        
        cells <- intersect(
          rownames(pt_s5j)[keep],
          colnames(s5j_expr)
        )
        
        data.frame(
          cell_id = cells,
          gene = gene,
          lineage = colnames(pt_s5j)[j],
          pseudotime = pt_s5j[cells, j],
          weight = wt_s5j[cells, j],
          expression = as.numeric(
            s5j_expr[gene, cells]
          )
        )
      }
    ) %>%
      dplyr::bind_rows()
  }
) %>%
  dplyr::bind_rows()

readr::write_tsv(
  s5j_source,
  file.path(
    source_dir,
    "FigS5J_gene_pseudotime.tsv.gz"
  )
)

###############################################################################
# Figs 3O–3Q
#
# These were already generated by Fig3_OPQ_tradeSeq_modules_clean.R.
# Copy the existing source-data files rather than regenerating them.
###############################################################################

opq_dir <- "Fig3_OPQ_tradeSeq_outputs"

opq_files <- c(
  "Fig3O_module_trajectory_source_data.tsv",
  "Fig3P_module_trajectory_source_data.tsv",
  "Fig3Q_module_trajectory_source_data.tsv",
  "Fig3O_tradeSeq_prebranch_modules.tsv",
  "Fig3P_tradeSeq_Lineage1_1_modules.tsv",
  "Fig3Q_tradeSeq_Lineage1_2_modules.tsv",
  "Fig3O-Q_module_annotations.tsv"
)

for (f in opq_files) {
  
  src <- file.path(opq_dir, f)
  
  if (file.exists(src)) {
    file.copy(
      src,
      file.path(source_dir, f),
      overwrite = TRUE
    )
  } else {
    warning("Missing existing source-data file: ", src)
  }
}


################################################################################
# sessionInfo()
################################################################################


