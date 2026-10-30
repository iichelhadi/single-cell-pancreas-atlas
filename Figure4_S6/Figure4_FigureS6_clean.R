###############################################################################
# Project: scRNA beta cell atlas — Lineage 2 reintegration and Figure 4 suite
# File: Figure4_FigureS6.R
# Author: Elhadi Iich
# Date: 2026-09-07
#
# Description:
#   End-to-end workflow for lineage 2 reintegration and visualization:
#   reclustering, supplementary figures, Slingshot overlays, INS trends,
#   branched heatmaps, GO enrichment, module detection and trajectories,
#   pseudobulk DESeq2, module scores along pseudotime, dot plots, and
#   metadata associations.
#
# Inputs:
#   - lin2.beta.reintegrated.rds
#   - cds_lin2.rds
#   - slingshot_sim_lin2.rds
#   - Beta_T2D_lineage2_1.rds, Beta_T2D_lineage2_2.rds
#   - Lineage2_sling_perc_expr.tsv
#   - ck1/ck2 RDS objects for branched heatmaps
#   - pseudobulk cts/coldata TSVs for DESeq2
#
###############################################################################

################################################################################
### 1. Environment Setup & Library Loading
################################################################################

suppressPackageStartupMessages({
  # Data manipulation and I/O
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(data.table)
  
  # Visualization
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(flextable)
  library(officer)
  
  # Parallel computing
  library(parallel)
  
  # Single-cell and trajectory analysis
  library(Seurat)
  library(monocle3)
  library(slingshot)
  library(reticulate)
  
  # Clustering and statistics
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(matrixStats)
  library(fastcluster)
  library(dynamicTreeCut)
})

set.seed(123)

reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')

palette <- c("#5A5156FF", "#E4E1E3FF", "#F6222EFF",
             "#FE00FAFF", "#16FF32FF" ,"#3283FEFF",
             "#FEAF16FF", "#B00068FF", "#1CFFCEFF",
             "#90AD1CFF", "#ffb300" ,"#DEA0DFFF",
             "#AA0EDFFF", "#F8A19FFF", "#F62424" ,
             "#C4451CFF", "#1C8356FF" ,"#85660DFF",
             "#B10DA1FF", "#FBE426FF", "#1CBE4FFF")

################################################################################
### 2. Data Loading
################################################################################
data.integrated <-
  readRDS('lin2.beta.reintegrated.rds')
df <- data.integrated@meta.data %>% as.data.frame()
cds <- readRDS('cds_lin2.rds')

perc_expr <- read_tsv('Lineage2_sling_perc_expr.tsv')

colnames(perc_expr) <- c('feature_id', 'perc_expr')

desired_order = c('1','4','3','5','2')

################################################################################
### 3. Reclustering (Figure 4B–4D)
################################################################################
### Figure 4B
B4 <- plot_cells(cds, show_trajectory_graph = F,
                 color_cells_by = 'monocle_clusters2',
                 label_cell_groups=T, cell_size = 1,
                 group_label_size = 10)+
  ggpubr::clean_theme()+ 
  scale_color_manual(values = c('1'='#F8766D',
                               '2'='#C49A00',
                               '3'='#53B400',
                               '4'='#00C094',
                               '5'='#00B6EB'))
B4

##Figure 4C ####
C4 <- df %>%
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
    x = "Cluster",
    y = "Age (years)"
  ) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  ) + 
  scale_fill_manual(values = c('1'='#F8766D',
                               '2'='#C49A00',
                               '3'='#53B400',
                               '4'='#00C094',
                               '5'='#00B6EB'))

### Figure 4D
D4 <- df %>%
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
  ) + 
  scale_fill_manual(values = c('1'='#F8766D',
                               '2'='#C49A00',
                               '3'='#53B400',
                               '4'='#00C094',
                               '5'='#00B6EB'))


(C4/D4)

################################################################################
### 4. Supplementary Figure S6A–S6F
################################################################################
##Figure S6A ####
S6A <- plot_cells(cds, show_trajectory_graph = F,
                  color_cells_by = 'Study',
                  label_cell_groups=FALSE, cell_size = 1,
                  # group_cells_by = 'cluster',
                  group_label_size = 4)+
  ggpubr::clean_theme()+NoLegend()+
  scale_color_manual(values = palette[-c(12,13,14)])

S6B <- df %>%
  group_by(monocle_clusters2, Study) %>%
  dplyr::count() %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>% 
  group_by(monocle_clusters2) %>% 
  mutate(per = paste0(round(100 *n/sum(n),1),'%')) %>%
  ggplot(aes(x = monocle_clusters2, y = n, fill = Study)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
  scale_x_discrete(limits = desired_order) + 
  # geom_text(aes( label = per),position = position_fill(vjust = 0.5)) +
  ylab('Proportion') +xlab("Cluster")+
  scale_fill_manual(values = palette[-c(12,13,14)])+
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text =  element_text(face = 'bold',size = 18), 
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 14))

##Figure S6C ####
S6C <- plot_cells(cds, show_trajectory_graph = F, 
                  color_cells_by = 'Condition',
                  label_cell_groups=FALSE, cell_size = 1,
                  # group_cells_by = 'cluster', 
                  group_label_size = 4)+
  ggpubr::clean_theme()+NoLegend()

##Figure S6D ####
S6D <- df %>%
  group_by(monocle_clusters2, Condition) %>%
  dplyr::count() %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>% 
  group_by(monocle_clusters2) %>% 
  mutate(per = paste0(round(100 *n/sum(n),1),'%')) %>%
  ggplot(aes(x = monocle_clusters2, y = n, fill = Condition)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
  scale_x_discrete(limits = desired_order)+
  ylab('Proportion') +xlab("Cluster") + 
  geom_text(aes( label = per),size = 8,position = position_fill(vjust = 0.5))+
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text =  element_text(face = 'bold',size = 18), 
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

## Figures S6A-D ####
(S6A/S6B)|(S6C/S6D)

##Figure S6E
# S6E <- df %>% drop_na(Hba1c) %>%
#   mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
#   group_by(monocle_clusters2) %>%
#   summarise(
#     median_value = median(Hba1c, na.rm = TRUE),
#     mad_value = mad(Hba1c, na.rm = TRUE)
#   ) %>%
#   mutate(
#     Hba1c_class = case_when(
#       median_value < 5.7 ~ "ND",
#       median_value >= 5.7 & median_value < 6.5 ~ "pre-T2D",
#       median_value >= 6.5 ~ "T2D"
#     )
#   ) %>%
#   ggplot(aes(x = monocle_clusters2, 
#              y = median_value, 
#              fill = Hba1c_class)) +
#   geom_bar(
#     stat = "identity",
#     position = position_dodge(),
#     color = 'black',
#     width = 0.8
#   ) +
#   geom_errorbar(
#     aes(ymin = median_value, #- mad_value, 
#         ymax = median_value + mad_value),
#     width = 0.2,
#     position = position_dodge(0.9)
#   ) + ylim(0, 9) +
#   scale_fill_manual(values = c('#F8766D','#FEAF16FF','#619CFF'))+
#   geom_text(
#     aes(label = paste0(round(median_value, 1), '%')),
#     size = 8, vjust = 2,
#     position = position_dodge(0.8)) +
#   ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
#   scale_x_discrete(limits = desired_order) +
#   labs(x = "Cluster", y = "Median HbA1c (%)")+ 
#   guides(fill = guide_legend(title = "HbA1c class"))+
#   theme(axis.title.y = element_text(face = 'bold', size = 16),
#         axis.text =  element_text(face = 'bold', size = 18), 
#         axis.title.x = element_text(face = 'bold', size = 16),
#         legend.text = element_text(face = 'bold', size = 14),
#         legend.title = element_text(face = 'bold', size = 14))



################################################################################
### Slingshot Curves & UMAP Labels
################################################################################

sim <- readRDS(
  "slingshot_sim_lin2.rds"
)

lineages <- slingshot::slingCurves(SlingshotDataSet(sim))

# Assuming 'lineages' is already obtained from slingshot::slingCurves
lineage_data <- lapply(names(lineages), function(lineage_name) {
  curve_df <- as.data.frame(lineages[[lineage_name]]$s)
  colnames(curve_df) <- c("Dim-1", "Dim-2")  # Rename columns based on your structure
  curve_df$lineage <- lineage_name
  curve_df$path_id <- paste(lineage_name, seq_len(nrow(curve_df)), sep="_")
  return(curve_df)
})

# Combine all lineage data frames into one
combined_lineages <- do.call(rbind, lineage_data)
combined_lineages$lineage <- ifelse(combined_lineages$lineage == 'Lineage1',
                                    'Lineage2_1','Lineage2_2')


# Step 1: Get UMAP coordinates and clusters
umap_coords <- reducedDims(cds)$UMAP

umap_df <- data.frame(
  Dim_1 = umap_coords[, 1],
  Dim_2 = umap_coords[, 2],
  cluster = clusters(cds, reduction_method = "UMAP")
)

# Step 2: Calculate cluster centers for label positions
label_coords <- aggregate(
  cbind(Dim_1, Dim_2) ~ cluster, 
  data = umap_df, 
  FUN = mean
)

# Step 3: Rename columns in combined_lineages for consistency
colnames(combined_lineages)[colnames(combined_lineages) == "Dim-1"] <- "Dim_1"
colnames(combined_lineages)[colnames(combined_lineages) == "Dim-2"] <- "Dim_2"

# Step 4: Build the plot
B4 <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = 'cluster',
  label_groups_by_cluster = FALSE,
  label_cell_groups = FALSE, # turn off default labels so we add our own
  cell_size = 1,
  group_label_size = 10
) + 
  geom_path(
    data = combined_lineages,
    aes(
      x = Dim_1,
      y = Dim_2,
      color = lineage,
      group = lineage
    ),
    size = 2
  ) + 
  scale_color_manual(
    values = c(
      '1' = '#F8766D',
      '2' = '#C49A00',
      '3' = '#53B400',
      '4' = '#00C094',
      '5' = '#00B6EB',
      'Lineage2_1' = '#F62424',
      'Lineage2_2' = '#ffb300'
    )
  ) + 
  geom_text(
    data = label_coords,
    aes(x = Dim_1, y = Dim_2, label = cluster),
    size = 10
  )

# Step 5: Show the plot
B4+ggpubr::clean_theme()+ 
  theme(legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

################################################################################
### 5. Slingshot Branch Plots & INS Expression (Figures 4E–4H)
################################################################################

E4 <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = 'Lineage2_1',
  cell_size = 1,
  group_label_size = 10
)+ggpubr::clean_theme()+
  theme(legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

F4 <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = 'Lineage2_2',
  cell_size = 1,
  group_label_size = 10
)+ggpubr::clean_theme()+
  theme(legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

E4/F4


lin2_1 <-
  readRDS(
    'Beta_T2D_lineage2_1.rds'
  )


lin2_2 <-
  readRDS(
    'Beta_T2D_lineage2_2.rds'
  )


### Figure 4G
exp2_1 <- lin2_1@assays$RNA@data
exp2_1 <- exp2_1 %>% t()

meta <- lin2_1@meta.data[, c('Lineage2_1', 'monocle_clusters2')]
meta <- as.data.frame(meta)
colnames(meta)[1:2] <- c('Lineage2_1', 'Clusters')
all(rownames(meta) %in% rownames(exp2_1))

plt_data <- cbind(meta, exp2_1)

mycolors <- c('#F8766D','#C49A00', '#53B400','#00C094')

G4 <- ggscatter(
  data = plt_data,
  x = 'Lineage2_1',
  y = 'INS',
  color = 'Clusters',
  size = 0.6
) +
  geom_smooth(
    se = FALSE,
    color = 'black',
    method = 'lm',
    formula = y ~ splines::ns(x, df = 3)
  ) +
  ggpubr::theme_pubr() +
  scale_color_manual(values = mycolors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  xlab('Pseudotime') +
  theme(
    legend.position = 'right',
    axis.title.y  = element_text(face = 'bold', size = 16),
    axis.title.x  = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 16),
    legend.title  = element_text(face = 'bold', size = 14),
    legend.text   = element_text(
      face = 'bold',
      size = 14,
      hjust = 1,
      vjust = 1
    )
  )
print(G4)

### Figure 4H
exp2_2 <- lin2_2@assays$RNA@data
exp2_2 <- exp2_2 %>% t()

meta <- lin2_2@meta.data[, c('Lineage2_2', 'monocle_clusters2')]
meta <- as.data.frame(meta)
colnames(meta)[1:2] <- c('Lineage2_2', 'Clusters')
all(rownames(meta) %in% rownames(exp2_2))

plt_data <- cbind(meta, exp2_2)

mycolors <- c('#F8766D', '#53B400', '#00C094','#00B6EB')

H4 <- ggscatter(
  data = plt_data,
  x = 'Lineage2_2',
  y = 'INS',
  color = 'Clusters',
  size = 0.6
) +
  geom_smooth(
    se = FALSE,
    color = 'black',
    method = 'lm',
    formula = y ~ splines::ns(x, df = 3)
  ) +
  ggpubr::theme_pubr() +
  scale_color_manual(values = mycolors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  xlab('Pseudotime') +
  theme(
    legend.position = 'right',
    axis.title.y  = element_text(face = 'bold', size = 16),
    axis.title.x  = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 16),
    legend.title  = element_text(face = 'bold', size = 14),
    legend.text   = element_text(
      face = 'bold',
      size = 14,
      hjust = 1,
      vjust = 1
    )
  )
print(H4)

################################################################################
### 6. Heatmap Data Prep, Drawing & Legends (Figure S6E)
################################################################################

### Figure S6E

coldata <- data.frame(
  study        = sim$study,
  Cluster      = sim$monocle_clusters2,
  Lineage2_1   = sim$slingPseudotime_1,
  Lineage2_2   = sim$slingPseudotime_2
)
rownames(coldata) <- colnames(sim)

# 2) Source helper functions -----------------------------------------
source(
  "custom_ClusterData.R"
)
source(
  "custom_visCluster.R"
)

# 3) Genes of interest + random subset -------------------------------
genes_of_interest <- c('INS','MAFA','ABCC8','IAPP', 'G6PC2', 'PAX6',
                       'HSPA1A','HSPA1B','HSPA5','MT1E','JUND',
                       'PRSS1','CTRB2','REG1A','CPA1','GCG','CLU')


lin2_1_obj <- readRDS(
  "Beta_T2D_lineage2_1.rds"
)
mat1_raw <- lin2_1_obj@assays$SCT@data

# 4) Branch 1: order -> smooth (parallel) -> z-score ---------------------
cells1   <- intersect(colnames(mat1_raw), rownames(coldata))
ord1     <- order(coldata[cells1, "Lineage2_1"])
cells1_o <- cells1[ord1]
cells1_rev <- rev(cells1_o)
pt1        <- coldata[cells1_rev, "Lineage2_1"]
mat1_ord   <- mat1_raw[, cells1_rev]


# Define a “never-fail” smooth helper
safe_smooth <- function(x, df = 3) {
  tryCatch(
    stats::smooth.spline(x, df = df)$y,
    error = function(e) rep(NA_real_, length(x))
  )
}

ncores   <- min(detectCores(), 4)
# Branch 1
sm1_list <- mclapply(
  seq_len(nrow(mat1_ord)),
  function(i) safe_smooth(mat1_ord[i, ], df = 3),
  mc.cores = ncores
)

sm1 <- do.call(rbind, sm1_list)
rownames(sm1) <- rownames(mat1_ord)
mat1_z <- t(scale(t(sm1)))
colnames(mat1_z) <- colnames(mat1_ord)
# saveRDS(mat1_z, 'mat1_z.rds') ### mat1_z generated from RNA assay for module detection and SCT for heatmap


# 5) Branch 2: same steps -----------------------------------------------
lin2_2_obj <- readRDS(
  "Beta_T2D_lineage2_2.rds"
)
mat2_raw <- lin2_2_obj@assays$SCT@data

cells2   <- intersect(colnames(mat2_raw), rownames(coldata))
ord2     <- order(coldata[cells2, "Lineage2_2"])
cells2_o <- cells2[ord2]
pt2      <- coldata[cells2_o, "Lineage2_2"]
mat2_ord <- mat2_raw[, cells2_o]

sm2_list <- mclapply(
  seq_len(nrow(mat2_ord)),
  function(i) safe_smooth(mat2_ord[i, ], df = 3),
  mc.cores = 3
)

sm2 <- do.call(rbind, sm2_list)
rownames(sm2) <- rownames(mat2_ord)
mat2_z <- t(scale(t(sm2)))
colnames(mat2_z) <- colnames(mat2_ord)
# saveRDS(mat2_z, 'mat2_z.rds') ### mat2_z generated from RNA assay for module detection and SCT for heatmap

# 6) Cluster rows -------------------------------------------------------
#! enabled: generate the clustering once if the saved objects are absent
#! (save names now match the readRDS() calls below; previously '20250624_' prefix)
if (!file.exists('ck1_branched_heatmap.rds')) {
  ck1 <- clusterData_custom(exp = mat1_z, scaleData = TRUE, min.std = 0)
  ck2 <- clusterData_custom(exp = mat2_z, scaleData = TRUE, min.std = 0)
  saveRDS(ck1, 'ck1_branched_heatmap.rds')
  saveRDS(ck2, 'ck2_branched_heatmap.rds')
}

ck1 <- readRDS('ck1_branched_heatmap.rds')
ck2 <- readRDS('ck2_branched_heatmap.rds')

# 7) Force shared row set & order --------------------------------------
shared     <- intersect(ck1$wide.res$gene, ck2$wide.res$gene)
common_all <- unique(c(shared, genes_of_interest))
row_order  <- ck1$wide.res$gene[ck1$wide.res$gene %in% common_all]

mat1_common <- ck1$wide.res %>%
  dplyr::filter(gene %in% row_order) %>%
  dplyr::arrange(match(gene, row_order)) %>%
  dplyr::select(-cluster) %>%
  column_to_rownames("gene") %>%
  as.matrix()

mat2_common <- ck2$wide.res %>%
  dplyr::filter(gene %in% row_order) %>%
  dplyr::arrange(match(gene, row_order)) %>%
  dplyr::select(-cluster) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# 8) Build top annotations ---------------------------------------------
pt1_fun <- colorRamp2(
  seq(min(pt1), max(pt1), length.out = 100),
  plasma(100)
)

col_annot1 <- HeatmapAnnotation(
  Cluster  = anno_simple(as.character(coldata[cells1_rev,"Cluster"]),
                         col = setNames(c("#F8766D", "#C49A00", "#53B400","#00C094"),
                                        sort(unique(coldata[cells1_rev,"Cluster"])))),
  Pseudotime = anno_simple(pt1, col = pt1_fun),
  Lineage = anno_simple(rep(x = 'Lineage2_1',dim(coldata[cells1_rev,])[1]),
                        col = setNames("#F62424",nm = 'Lineage2_1')),
  annotation_height = unit(c(4,4,4,4), "mm"),
  annotation_label = "",
  annotation_name_gp = gpar(fontsize=12,fontface = 'bold')
)

pt2_fun <- colorRamp2(
  seq(min(pt2), max(pt2), length.out = 100),
  plasma(100)
)
col_annot2 <- HeatmapAnnotation(
  Cluster  = anno_simple(as.character(coldata[cells2_o,"Cluster"]),
                         col = setNames(c("#F8766D", "#53B400", "#00C094", "#00B6EB"),
                                        sort(unique(coldata[cells2_o,"Cluster"])))),
  Pseudotime = anno_simple(pt2, col = pt2_fun),
  Lineage = anno_simple(rep(x = 'Lineage2_2',dim(coldata[cells2_o,])[1]),
                        col = setNames("#ffb300",nm = 'Lineage2_2')),
  annotation_height = unit(c(4,4,4,4), "mm"),
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

# 9) Row markers on right panel only -------------------------------------
pos_interest <- which(rownames(mat2_common) %in% genes_of_interest)
row_mark     <- rowAnnotation(
  mark = anno_mark(
    at     = pos_interest,
    labels = rownames(mat2_common)[pos_interest],
    labels_gp = gpar(fontsize=14)
  )
)

# 10) Build & draw side‐by‐side -----------------------------------------
col_fun <- colorRamp2(c(-2,0,2), c("#0000FF","#00FFFF","#FF0000"))

ht1 <- Heatmap(
  mat1_common, name                = NULL,
  col                 = col_fun,
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
  mat2_common, name                = "Expression (z)",
  col                 = col_fun,
  top_annotation      = col_annot2,
  right_annotation    = row_mark,
  cluster_rows        = FALSE,
  cluster_columns     = FALSE,
  show_row_names      = FALSE,
  show_column_names   = FALSE,
  show_row_dend       = FALSE,
  show_heatmap_legend = TRUE,
  # column_title        = "Branched Trajectory",
  border              = TRUE
)

# Cairo::CairoPNG("20250625_combined_heatmaps4.png",
#     width  = 10*300, height = 8*300, res = 300)
draw(ht1 + ht2,
     merge_legend        = TRUE,
     heatmap_legend_side = "right")
dev.off()


library(ComplexHeatmap)
library(circlize)
library(grid)

# pdf("combined_heatmap_legend.pdf", width = 5, height = 6)


# Cluster legend (may have different clusters in each lineage)
clusters <- sort(unique(c(coldata[cells1_rev, "Cluster"], coldata[cells2_o, "Cluster"])))
cluster_colors <- setNames(c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB")[seq_along(clusters)], clusters)

Cluster_legend <- Legend(
  labels     = names(cluster_colors),
  title      = "Cluster",
  legend_gp  = gpar(fill = cluster_colors)
)

# Lineage legend
lineage_colors <- c("Lineage2_1" = "#F62424", "Lineage2_2" = "#ffb300")
Lineage_legend <- Legend(
  labels     = names(lineage_colors),
  title      = "Lineage",
  legend_gp  = gpar(fill = lineage_colors)
)

# Determine actual pseudotime range from both branches
pt_min <- min(c(pt1, pt2), na.rm = TRUE)
pt_max <- max(c(pt1, pt2), na.rm = TRUE)

# Use a plasma gradient over the actual range
pseudotime_fun <- colorRamp2(
  seq(pt_min, pt_max, length.out = 100),
  plasma(100)
)

Pseudotime_legend <- Legend(
  col_fun        = pseudotime_fun,
  title          = "Pseudotime",
  legend_height  = unit(3, "cm"),
  at             = round(c(0,5,10,15,20), 1)  # nicely spaced ticks
)


# Expression legend (z-score)
expr_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "#00FFFF", "#FF0000"))
Expression_legend <- Legend(
  col_fun    = expr_fun,
  title      = "Expression (z)",
  legend_height = unit(2, "cm")
)

# Combine and draw all
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
### Figure S6F
pseudotime_start_end_association_filt <-
  read_tsv(
    'lin2_1_pseudotime_start_end_association_filt.tsv'
  )

genelist <- pseudotime_start_end_association_filt$logFClineage1
names(genelist) <- pseudotime_start_end_association_filt$feature_id
genelist <- genelist[abs(genelist) > 1]
genelist <- genelist[order(genelist,decreasing = T)]

set.seed(123)
gse <- clusterProfiler::gseGO(geneList = genelist,
                              ont = 'BP',
                              OrgDb = org.Hs.eg.db,
                              keyType = 'SYMBOL',
                              eps = 0,
                              minGSSize = 25,
                              maxGSSize = 500,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = 'BH')

gse <- clusterProfiler::simplify(gse)

gse@result %>% as.data.frame() %>%
  filter(
    ID %in% c(
      'GO:1901700',
      'GO:0045595',
      'GO:0002250',
      'GO:0042127',
      'GO:0010038',
      'GO:0006873'
    )
  ) %>% 
  ggplot(., aes( x = NES, y = reorder(Description, NES),
                 fill = p.adjust)) + ylab('Term') +
  geom_bar(stat = 'identity') + ggpubr::theme_pubr() +
  scale_fill_continuous(
    low = "red",high = "blue",
    guide = guide_colorbar(reverse = TRUE))+
  theme(
    axis.title.y = element_text(face = 'bold', size = 14),
    axis.text =  element_text(size = 16),
    axis.text.y = element_text(face = 'bold', size = 14),
    axis.title.x = element_text(face = 'bold', size = 15),
    legend.title = element_text(face = 'bold', size = 14),
    legend.text = element_text(face = 'bold', size = 14,angle = 45,
                               hjust = 1,vjust = 1))


### Figure S6G
pseudotime_start_end_association_filt <-
  read_tsv(
    'lin2_2_pseudotime_start_end_association_filt.tsv'
  )

genelist <- pseudotime_start_end_association_filt$logFClineage1
names(genelist) <- pseudotime_start_end_association_filt$feature_id
genelist <- genelist[abs(genelist) > 1]
genelist <- genelist[order(genelist,decreasing = T)]


set.seed(123)
gse2_2 <- clusterProfiler::gseGO(geneList = genelist,
                                 ont = 'BP',
                                 OrgDb = org.Hs.eg.db,
                                 keyType = 'SYMBOL',
                                 eps = 0,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 pvalueCutoff = 0.1,
                                 pAdjustMethod = 'BH')

gse2_2 <- clusterProfiler::simplify(gse2_2)

gse2_2@result %>% as.data.frame() %>%
  filter(ID %in% c('GO:0032940',
                   'GO:0009725',
                   'GO:1901700',
                   'GO:0010628',
                   'GO:0042127',
                   'GO:0016042'
  )) %>%
  ggplot(., aes( x = NES, y = reorder(Description, NES),
                 fill = p.adjust)) + ylab('Term') +
  geom_bar(stat = 'identity') + ggpubr::theme_pubr() +
  scale_fill_continuous(
    low = "red",high = "blue",
    guide = guide_colorbar(reverse = TRUE))+
  theme(
    axis.title.y = element_text(face = 'bold', size = 14),
    axis.text =  element_text(size = 16),
    axis.text.y = element_text(face = 'bold', size = 14),
    axis.title.x = element_text(face = 'bold', size = 15),
    legend.title = element_text(face = 'bold', size = 14),
    legend.text = element_text(face = 'bold', size = 14,angle = 45,
                               hjust = 1,vjust = 1))

################################################################################
### 7. Module trajectories and GO enrichment — Figures 4I–4K
################################################################################

# Revised analysis is implemented in:
# Figure4_IJK_tradeSeq_modules_clean.R
#
# The revised workflow uses tradeSeq-fitted GAM trajectories, separates the
# shared pre-bifurcation and lineage-specific post-bifurcation regions, selects
# module number by silhouette analysis, and performs GO enrichment using the
# pseudotime-associated gene set as the background.

################################################################################
### 8. Pseudobulk DESeq2 Analysis (Figures 4L)
################################################################################

library(DESeq2)

cts <-
  read.table(
    'beta_C2_v_C5_pseudobulk_cts.tsv',
    header = T,
    sep = '\t',
    row.names = 1
  )
colnames(cts) <- gsub(colnames(cts),pattern = '[.]',replacement = '-')

coldata <- read.table(
  'beta_C2_v_C5_pseudobulk_coldata.tsv',
  header = T,
  sep = '\t',
  row.names = 1
)

coldata$monocle_clusters2 <- as.character(coldata$monocle_clusters2)

perc_expr <- read_tsv('Lineage2_sling_perc_expr.tsv')
dim(perc_expr)
colnames(perc_expr) <- c('feature_id','perc_expr')

cts <- cts[rownames(cts) %in% perc_expr$feature_id[perc_expr$perc_expr >= 5],]

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ study + monocle_clusters2)


# Calculate size factors proportional to cell counts
size_factors <- coldata$counts / median(coldata$counts)

# Ensure the size factors have names matching the sample names in 'dds'
names(size_factors) <- colnames(dds)

# Assign the size factors to 'dds'
sizeFactors(dds) <- size_factors

###

dds <- DESeq(dds,parallel = F)

which(mcols(dds)$betaConv) %>% length()
dds <- dds[which(mcols(dds)$betaConv),]

resultsNames(dds)
res <-
  results(dds,
          contrast = c('monocle_clusters2', '2', '5'),
          alpha = 0.05)
summary(res)

res_filt <- res %>% as.data.frame() %>% filter(complete.cases(.))
res_filt$genes <- rownames(res_filt)

res_filt <- left_join(res_filt, perc_expr, 
                      by = c('genes'='feature_id'))

res_filt_2 <- res_filt
res_filt_2$change <- ifelse((res_filt_2$log2FoldChange > 1 &
                               res_filt_2$padj <= 0.05 ),
                            'UP',
                            ifelse((
                              res_filt_2$log2FoldChange < -1 &
                                res_filt_2$padj <= 0.05
                            ), 'DOWN', 'NS'))

### Figure 4L
gene_sel <- c('REG1A','REG1B', 'TM4SF4', 'S100A10', 'PRSS1', 'CPA1', 'SPINK1', 
              'CTRB2', 'NKX6-1', 'HSPA5', 'DDIT3', 
              'PDK4', 'INS', 'IAPP', 'G6PC2', 'MAFA','PAX6',
              "EDEM1", "HSPA5", "DDIT3", "EIF2AK3",  # significant
              "ATF4", "XBP1", "ATF6", "ERN1")        # master regulators)


# Create the plot
ggplot(res_filt_2, aes(x = log2FoldChange, y = -log10(padj))) +
  # Plot non-highlighted points
  geom_point(data = subset(res_filt_2, !genes %in% gene_sel), 
             aes(size = perc_expr), color = "grey") +
  # Plot highlighted points
  geom_point(data = subset(res_filt_2, genes %in% gene_sel), 
             aes(size = perc_expr, color = change)) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue")) +  # Manually set colors
  # Add vertical and horizontal lines
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) +
  # Add labels for highlighted points
  ggrepel::geom_text_repel(
    data = subset(res_filt_2, genes %in% gene_sel),
    aes(label = genes),
    max.overlaps = 60,
    box.padding = 0.5, 
    size = 5,
    fontface = "bold"
  ) +
  guides(color = 'none')+
  # Set y-axis limits
  # ylim(0, 6) +
  # Apply publication-ready theme
  ggpubr::theme_pubr(legend = 'right') +
  # Override legend aesthetics to ensure black color
  # guides(size = guide_legend(override.aes = list(color = "black")))+
  scale_size_continuous(name = "Percent Expression")+ 
  xlab('Log2 Fold Change') + ylab('Log10 padj')+
  guides(color = 'none') +
  ggpubr::theme_pubr(legend = 'right')+
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text =  element_text(face = 'bold', size = 15),
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14))


################################################################################
### 10.             Cluster violin plots for Figures 4M–4N
################################################################################

### Figure 4M
# Define gene sets for module scores
lin2_obj <-
  readRDS(
    "lin2.beta.reintegrated.rds"
  )

er_stress_genes <- c("HSPA5","DDIT3","EIF2AK3","XBP1","ATF4","HERPUD1",
                     "HSP90B1","DNAJC3","EDEM1")

exocrine_genes <- c("REG1B", "REG1A", "PRSS1", "SPINK1","CPA1","PRSS2","CPB1",
                    'AMY2A', 'CEL', 'PNLIP','CTRB2')

core_genes     <- c("INS","PDX1","NKX6-1","PCSK1")

maturity_genes <- c("MAFA","G6PC2","IAPP","SLC30A8","ABCC8",'UCN3')

modules <- list('Core' = core_genes, 'Maturity' = maturity_genes, 
                'ER Stress'= er_stress_genes, 'Exocrine'= exocrine_genes)

# Calculate module scores
lin2_obj <- AddModuleScore(
  lin2_obj,
  features = list(er_stress = er_stress_genes),
  name = "ER_Stress",
  assay = "RNA"
)
lin2_obj <- AddModuleScore(
  lin2_obj,
  features = list(exocrine = exocrine_genes),
  name = "Exocrine",
  assay = "RNA"
)
lin2_obj <- AddModuleScore(
  lin2_obj,
  features = list(Core = core_genes),
  name = "Core",
  assay = "RNA"
)
lin2_obj <- AddModuleScore(
  lin2_obj,
  features = list(Maturity = maturity_genes),
  name = "Maturity",
  assay = "RNA"
)

colnames(lin2_obj@meta.data)[64] <- "ER_Stress"
colnames(lin2_obj@meta.data)[65] <- "Exocrine"
colnames(lin2_obj@meta.data)[66] <- "Core"
colnames(lin2_obj@meta.data)[67] <- "Maturity"

p1 <- VlnPlot(
  lin2_obj, features = c("ER_Stress"), assay = "RNA",
  group.by = "monocle_clusters2", ncol = 1, pt.size = 0
) +
  scale_x_discrete(limits = c("1", "4", "3", "5", "2")) +
  scale_fill_manual(values = c(
    "1" = "#F8766D",
    "2" = "#C49A00",
    "3" = "#53B400",
    "4" = "#00C094",
    "5" = "#00B6EB"
  )) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, fill = 'white') +
  ylab("Module Score") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank()
  ) +
  NoLegend()

p2 <- VlnPlot(
  lin2_obj, features = c("Exocrine"), assay = "RNA",
  group.by = "monocle_clusters2", ncol = 1, pt.size = 0
) +
  scale_x_discrete(limits = c("1", "4", "3", "5", "2")) +
  scale_fill_manual(values = c(
    "1" = "#F8766D",
    "2" = "#C49A00",
    "3" = "#53B400",
    "4" = "#00C094",
    "5" = "#00B6EB"
  )) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, fill = 'white') +
  xlab("Cluster") +
  ylab("Module Score") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 14)
  ) +
  NoLegend()

p3 <- VlnPlot(
  lin2_obj, features = c("Core"), assay = "RNA",
  group.by = "monocle_clusters2", ncol = 1, pt.size = 0
) +
  scale_x_discrete(limits = c("1", "4", "3", "5", "2")) +
  scale_fill_manual(values = c(
    "1" = "#F8766D",
    "2" = "#C49A00",
    "3" = "#53B400",
    "4" = "#00C094",
    "5" = "#00B6EB"
  )) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, fill = 'white') +
  ylab("Module Score") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  NoLegend()

p4 <- VlnPlot(
  lin2_obj, features = c("Maturity"), assay = "RNA",
  group.by = "monocle_clusters2", ncol = 1, pt.size = 0
) +
  scale_x_discrete(limits = c("1", "4", "3", "5", "2")) +
  scale_fill_manual(values = c(
    "1" = "#F8766D",
    "2" = "#C49A00",
    "3" = "#53B400",
    "4" = "#00C094",
    "5" = "#00B6EB"
  )) +
  xlab("Cluster") +
  ylab("Module Score") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, fill = 'white') +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 14)
  ) +
  NoLegend()

# (p3/p4) | (p1/p2)

################################################################################
# Lineage-long pseudotime data for Figures 4M–4N
################################################################################

score_df <- lin2_obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(
    cell,
    monocle_clusters2,
    Core,
    Maturity,
    ER_Stress,
    Exocrine
  )

sce <- readRDS('sce_lin2.rds')

crv_df <- as.data.frame(
  colData(sce)$crv
) %>%
  tibble::rownames_to_column("cell")

lin2_1_df <- crv_df %>%
  dplyr::transmute(
    cell = cell,
    Lineage = "Lineage2_1",
    Pseudotime = pseudotime.Lineage1,
    weight = cellWeights.Lineage1
  )

lin2_2_df <- crv_df %>%
  dplyr::transmute(
    cell = cell,
    Lineage = "Lineage2_2",
    Pseudotime = pseudotime.Lineage2,
    weight = cellWeights.Lineage2
  )

lineage_long <- dplyr::bind_rows(
  lin2_1_df,
  lin2_2_df
) %>%
  dplyr::filter(
    weight > 0,
    !is.na(Pseudotime)
  ) %>%
  dplyr::left_join(
    score_df,
    by = "cell"
  ) %>%
  dplyr::mutate(
    Cluster = as.character(monocle_clusters2),
    Lineage = factor(
      Lineage,
      levels = c("Lineage2_1", "Lineage2_2")
    )
  )

################################################################################
# Colors
################################################################################

cluster_colors <- c(
  "1" = "#F8766D",
  "2" = "#C49A00",
  "3" = "#53B400",
  "4" = "#00C094",
  "5" = "#00B6EB"
)

smoother_colors <- c(
  "Lineage2_1" = "#F62424",
  "Lineage2_2" = "#ffb300"
)

################################################################################
# Module-score trends along pseudotime
################################################################################

genes <- c(
  "Core",
  "Maturity",
  "ER_Stress",
  "Exocrine"
)

l <- c(3, 2, 2, 2.5)

plt_list <- list()

for (i in seq_along(genes)) {
  
  gene <- genes[i]
  
  plt_data <- lineage_long %>%
    dplyr::transmute(
      cell = cell,
      Pseudotime = Pseudotime,
      Lineage = Lineage,
      weight = weight,
      Cluster = Cluster,
      score = .data[[gene]]
    )
  
  p <- ggplot(
    plt_data,
    aes(
      x = Pseudotime,
      y = score
    )
  ) +
    geom_point(
      aes(
        color = Cluster,
        alpha = weight
      ),
      size = 0.6
    ) +
    geom_smooth(
      data = subset(
        plt_data,
        Lineage == "Lineage2_1"
      ),
      aes(
        group = Lineage,
        weight = weight
      ),
      linewidth = 2,
      color = smoother_colors["Lineage2_1"],
      se = FALSE,
      method = "lm",
      formula = y ~ splines::ns(x, df = 3)
    ) +
    
    geom_smooth(
      data = subset(
        plt_data,
        Lineage == "Lineage2_2"
      ),
      aes(
        group = Lineage,
        weight = weight
      ),
      linewidth = 2,
      color = smoother_colors["Lineage2_2"],
      se = FALSE,
      method = "lm",
      formula = y ~ splines::ns(x, df = 3)
    ) +
    
    scale_color_manual(
      name = "Cluster",
      values = cluster_colors,
      breaks = c("1", "4", "3", "5", "2")
    ) +
    
    scale_alpha_continuous(
      range = c(0.05, 0.25),
      guide = "none"
    ) +
    
    coord_cartesian(
      ylim = c(-0.5, l[i])
    ) +
    
    labs(
      x = "Pseudotime",
      y = gene
    ) +
    ggpubr::theme_pubr(
      legend = "right"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(
        face = "bold",
        size = 14
      ),
      axis.text = element_text(
        face = "bold",
        size = 14
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 14
      ),
      legend.title = element_text(
        face = "bold",
        size = 14
      ),
      legend.text = element_text(
        face = "bold",
        size = 14
      )
    ) +
    
    guides(
      color = guide_legend(
        override.aes = list(
          size = 3,
          alpha = 1
        )
      )
    )
  
  plt_list[[gene]] <- p
}

P <- wrap_plots(
  plt_list,
  guides = "collect",
  ncol = 2
)

# Figure 4M-4N composite
(p3 | P[[1]] + NoLegend()) / (p4 | P[[2]] + NoLegend()) |
  (p1 | P[[3]] + NoLegend()) / (p2 | P[[4]] + NoLegend())

################################################################################
### 11. Dot plot C2 vs C5 (Figure S6H)
################################################################################

SCpubr::do_DotPlot(
  subset(lin2_obj, subset = monocle_clusters2 %in% c("2", "5")),
  features = modules,
  assay = "RNA",
  group.by = "monocle_clusters2",
  use_viridis = TRUE,
  legend.length = 3,
  legend.position = "right",
  viridis.direction = 1,
  font.size = 16,
  ylab = "Cluster",
  viridis.palette = "F",
  scale = TRUE
) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 14)
  )

################################################################################
### 12. Dot plot C5 male vs female (Figure S6I)
################################################################################

C5 <- subset(lin2_obj, subset = monocle_clusters2 == "5")
dim(C5)
table(C5$Gender, useNA = "ifany")

C5 <- C5[, C5$Gender %in% c("male", "female") & !is.na(C5$Gender)]
table(C5$Gender, useNA = "ifany")

markers <- list(
  "Exocrine markers" = c("PRSS1", "PRSS2", "CPA1", "CPB1", "CTRB1",
                         "CTRB2", "REG1A", "REG1B", "SPINK1"),
  "EIF2α markers" = c("ESR1", "PPP1R15A", "PPP1R15B", "HSPA5", "DDIT3", "ATF4",
                      "EIF2B2", "EIF2S1", "EIF2AK3")
)

# Figure S6I
SCpubr::do_DotPlot(
  C5,
  features = markers,
  assay = "RNA",
  group.by = "Gender",
  use_viridis = TRUE,
  legend.length = 3,
  legend.position = "right",
  viridis.direction = 1,
  font.size = 16,
  xlab = "Cluster",
  viridis.palette = "F",
  scale = TRUE
) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 14)
  )

################################################################################
### Figure 4O — donor-level C2/C5 trait associations
################################################################################

# Revised donor-level analysis addressing Reviewer 1 Comment 9 and
# Reviewer 2 Comment 15 is implemented in:
# Figure4O_C2_C5_donor_traits_clean.R

################################################################################
### Figure 4P — IF donor sample meta data
################################################################################

library(flextable)
library(dplyr)
library(officer)

# Figure 4P donor table
donor_table <- tibble::tribble(
  ~`Donor ID`, ~`Age (years)`, ~Gender, ~`BMI (kg/m²)`, ~`Disease status`, ~`Years of diabetes`, ~CPA1, ~PRSS1,
  "ND00202",   55,             "F",     25.00,          "ND",              "–",                  "✓",   "✓",
  "ND09258",   96,             "F",     21.46,          "ND",              "–",                  "✓",   "✓",
  "ND11183",   46,             "F",     22.90,          "ND",              "–",                  "✓",   "✓",
  "ND13069",   46,             "M",     27.10,          "ND",              "–",                  "✓",   "✓",
  "ND13819",   NA,             NA,      NA,             "ND",              "–",                  "✓",   "–",
  "ND139817",  NA,             NA,      NA,             "T2D",             "4",                  "✓",   "–",
  "ND15112",   70,             "F",     42.20,          "T2D",             "6",                  "✓",   "✓",
  "ND14070",   66,             "M",     30.40,          "T2D",             "18",                 "✓",   "✓",
  "ND09756",   62,             "M",     32.69,          "T2D",             "20",                 "✓",   "✓",
  "ND18553",   67,             "M",     26.40,          "T2D",             "30",                 "✓",   "✓",
  "ND08838",   82,             "F",     31.27,          "T2D",             "32",                 "✓",   "✓"
)

ft <- flextable::flextable(donor_table) %>%
  flextable::theme_box() %>%
  flextable::bg(bg = "#D9D9D9", part = "header") %>%
  flextable::bg(bg = "white", part = "body") %>%
  flextable::color(color = "black", part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::padding(padding = 3, part = "all") %>%
  flextable::border_outer(
    border = officer::fp_border(color = "#7F7F7F", width = 1)
  ) %>%
  flextable::border_inner_h(
    border = officer::fp_border(color = "#A6A6A6", width = 0.75)
  ) %>%
  flextable::border_inner_v(
    border = officer::fp_border(color = "#A6A6A6", width = 0.75)
  ) %>%
  flextable::autofit()

ft

###############################################################################
# SOURCE DATA EXPORT — Figure 4 and Figure S6
###############################################################################

source_dir <- "../SourceData/Figure4_S6"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# Figure 4B, 4E, 4F, S6A and S6C
# Shared UMAP coordinates + plotted annotations
###############################################################################

umap_source <- as.data.frame(
  SingleCellExperiment::reducedDims(cds)$UMAP
)

colnames(umap_source)[1:2] <- c("UMAP_1", "UMAP_2")

umap_source <- umap_source %>%
  tibble::rownames_to_column("cell_id")

cds_meta_source <- as.data.frame(
  SummarizedExperiment::colData(cds)
) %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    dplyr::any_of(c(
      "monocle_clusters2",
      "Study",
      "Condition",
      "Lineage2_1",
      "Lineage2_2"
    ))
  )

fig4_umap_source <- umap_source %>%
  dplyr::left_join(
    cds_meta_source,
    by = "cell_id"
  )

readr::write_tsv(
  fig4_umap_source,
  file.path(
    source_dir,
    "Fig4B_E_F_S6A_C_UMAP.tsv.gz"
  )
)


###############################################################################
# Slingshot curves shown in Figure 4B
###############################################################################

fig4_slingshot_source <- combined_lineages %>%
  dplyr::select(
    lineage,
    path_id,
    Dim_1,
    Dim_2
  )

readr::write_tsv(
  fig4_slingshot_source,
  file.path(
    source_dir,
    "Fig4B_slingshot_curves.tsv"
  )
)


###############################################################################
# Figures 4C and 4D
# Cell-weighted age and BMI distributions
###############################################################################

fig4cd_source <- df %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    monocle_clusters2,
    Age,
    BMI
  )

readr::write_tsv(
  fig4cd_source,
  file.path(
    source_dir,
    "Fig4C_D_age_BMI.tsv.gz"
  )
)


###############################################################################
# Figures S6B and S6D
# Study and condition composition
###############################################################################

s6b_source <- df %>%
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

s6d_source <- df %>%
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
  dplyr::bind_rows(
    s6b_source,
    s6d_source
  ),
  file.path(
    source_dir,
    "FigS6B_D_composition.tsv"
  )
)


###############################################################################
# Figures 4G and 4H
# INS expression versus pseudotime
###############################################################################

ins_l21 <- Seurat::GetAssayData(
  lin2_1,
  assay = "RNA",
  layer = "data"
)["INS", ]

ins_l22 <- Seurat::GetAssayData(
  lin2_2,
  assay = "RNA",
  layer = "data"
)["INS", ]

fig4g_source <- data.frame(
  cell_id = colnames(lin2_1),
  lineage = "Lineage2_1",
  pseudotime = lin2_1@meta.data[
    colnames(lin2_1),
    "Lineage2_1"
  ],
  cluster = lin2_1@meta.data[
    colnames(lin2_1),
    "monocle_clusters2"
  ],
  INS = as.numeric(
    ins_l21[colnames(lin2_1)]
  )
)

fig4h_source <- data.frame(
  cell_id = colnames(lin2_2),
  lineage = "Lineage2_2",
  pseudotime = lin2_2@meta.data[
    colnames(lin2_2),
    "Lineage2_2"
  ],
  cluster = lin2_2@meta.data[
    colnames(lin2_2),
    "monocle_clusters2"
  ],
  INS = as.numeric(
    ins_l22[colnames(lin2_2)]
  )
)

readr::write_tsv(
  dplyr::bind_rows(
    fig4g_source,
    fig4h_source
  ),
  file.path(
    source_dir,
    "Fig4G_H_INS_pseudotime.tsv.gz"
  )
)


###############################################################################
# Figures S6F and S6G
# Exact GSEA terms plotted
###############################################################################

s6f_ids <- c(
  "GO:1901700",
  "GO:0045595",
  "GO:0002250",
  "GO:0042127",
  "GO:0010038",
  "GO:0006873"
)

s6g_ids <- c(
  "GO:0032940",
  "GO:0009725",
  "GO:1901700",
  "GO:0010628",
  "GO:0042127",
  "GO:0016042"
)

s6f_source <- as.data.frame(gse@result) %>%
  dplyr::filter(ID %in% s6f_ids) %>%
  dplyr::mutate(panel = "S6F_Lineage2_1")

s6g_source <- as.data.frame(gse2_2@result) %>%
  dplyr::filter(ID %in% s6g_ids) %>%
  dplyr::mutate(panel = "S6G_Lineage2_2")

readr::write_tsv(
  dplyr::bind_rows(
    s6f_source,
    s6g_source
  ),
  file.path(
    source_dir,
    "FigS6F_G_GSEA.tsv"
  )
)


###############################################################################
# Figures 4I–4K
# Already generated by Figure4_IJK_tradeSeq_modules_clean.R
###############################################################################

ijk_dir <- "Figure4_IJK_tradeSeq_outputs"

ijk_files <- c(
  "Figure4I_module_trajectory_source_data.tsv",
  "Figure4J_module_trajectory_source_data.tsv",
  "Figure4K_module_trajectory_source_data.tsv",
  "Figure4I_tradeSeq_prebranch_modules.tsv",
  "Figure4J_tradeSeq_Lineage2_1_modules.tsv",
  "Figure4K_tradeSeq_Lineage2_2_modules.tsv",
  "Figure4I-K_module_annotations.tsv"
)

for (f in ijk_files) {
  
  src <- file.path(ijk_dir, f)
  
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


###############################################################################
# Figure 4L
# Pseudobulk DESeq2 volcano plot
###############################################################################

readr::write_tsv(
  res_filt_2,
  file.path(
    source_dir,
    "Fig4L_pseudobulk_DE.tsv"
  )
)


###############################################################################
# Figures 4M–4N
# Cell-level module scores used for violin plots
###############################################################################

fig4mn_violin_source <- lin2_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    monocle_clusters2,
    Core,
    Maturity,
    ER_Stress,
    Exocrine
  )

readr::write_tsv(
  fig4mn_violin_source,
  file.path(
    source_dir,
    "Fig4M_N_module_scores.tsv.gz"
  )
)


###############################################################################
# Figures 4M–4N
# Module-score trajectories versus pseudotime
###############################################################################

fig4mn_trajectory_source <- lineage_long %>%
  dplyr::select(
    cell,
    Lineage,
    Pseudotime,
    weight,
    Cluster,
    Core,
    Maturity,
    ER_Stress,
    Exocrine
  )

readr::write_tsv(
  fig4mn_trajectory_source,
  file.path(
    source_dir,
    "Fig4M_N_module_pseudotime.tsv.gz"
  )
)


###############################################################################
# Figure S6H
# Dot plot: C2 versus C5
###############################################################################

s6h_genes <- unique(unlist(modules))

s6h_obj <- subset(
  lin2_obj,
  subset = monocle_clusters2 %in% c("2", "5")
)

s6h_source <- Seurat::DotPlot(
  s6h_obj,
  features = s6h_genes,
  group.by = "monocle_clusters2"
)$data %>%
  as.data.frame() %>%
  dplyr::rename(
    cluster = id,
    gene = features.plot,
    pct_expressing = pct.exp,
    average_expression = avg.exp,
    scaled_average_expression = avg.exp.scaled
  )

readr::write_tsv(
  s6h_source,
  file.path(
    source_dir,
    "FigS6H_dotplot.tsv"
  )
)


###############################################################################
# Figure S6I
# Dot plot: male versus female C5 cells
###############################################################################

s6i_genes <- unique(unlist(markers))

s6i_source <- Seurat::DotPlot(
  C5,
  features = s6i_genes,
  group.by = "Gender"
)$data %>%
  as.data.frame() %>%
  dplyr::rename(
    gender = id,
    gene = features.plot,
    pct_expressing = pct.exp,
    average_expression = avg.exp,
    scaled_average_expression = avg.exp.scaled
  )

readr::write_tsv(
  s6i_source,
  file.path(
    source_dir,
    "FigS6I_dotplot.tsv"
  )
)


###############################################################################
# Figure 4P
# Donor metadata table
###############################################################################

readr::write_tsv(
  donor_table,
  file.path(
    source_dir,
    "Fig4P_IF_donor_metadata.tsv"
  )
)

###############################################################################
# sessionInfo()
###############################################################################



