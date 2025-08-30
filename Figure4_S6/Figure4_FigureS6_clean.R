###############################################################################
# Project: scRNA beta cell atlas — Lineage 2 reintegration and Figure 4 suite
# File: Figure4_FigureS6.R
# Author: Elhadi Iich
# Date: 2025-08-27
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

sce <- readRDS('sce_lin2_1.rds')
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
C4 <- df %>% drop_na(Age) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(
    median_value = median(Age, na.rm = TRUE),
    mad_value = mad(Age, na.rm = TRUE)
  )  %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = monocle_clusters2)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(),
    color = 'black',
    width = 0.8
  ) +
  geom_errorbar(
    aes(ymin = median_value, #- mad_value, 
        ymax = median_value + mad_value),
    width = 0.2,
    position = position_dodge(0.9)
  ) + 
  ggpubr::theme_pubr(border = FALSE, legend = 'none',base_size = 10) + 
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median Age")+
  theme(axis.title.y = element_text(face = 'bold', size = 14),
        axis.text =  element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(face = 'bold', size = 14)) + 
  scale_fill_manual(values = c('1'='#F8766D',
                               '2'='#C49A00',
                               '3'='#53B400',
                               '4'='#00C094',
                               '5'='#00B6EB'))

### Figure 4D
D4 <- df %>% drop_na(BMI) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(
    median_value = median(BMI, na.rm = TRUE),
    mad_value = mad(BMI, na.rm = TRUE)
  )  %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = monocle_clusters2)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(),
    color = 'black',
    width = 0.8
  ) +
  geom_errorbar(
    aes(ymin = median_value, #- mad_value, 
        ymax = median_value + mad_value),
    width = 0.2,
    position = position_dodge(0.9)
  ) + 
  ggpubr::theme_pubr(border = FALSE, legend = 'none',base_size = 10) + 
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median BMI")+
  theme(axis.title.y = element_text(face = 'bold', size = 14),
        axis.text =  element_text(size = 14, face = 'bold'), 
        axis.title.x = element_text(face = 'bold', size = 14))+ 
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
S6E <- df %>% drop_na(Hba1c) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(
    median_value = median(Hba1c, na.rm = TRUE),
    mad_value = mad(Hba1c, na.rm = TRUE)
  ) %>%
  mutate(
    Hba1c_class = case_when(
      median_value < 5.7 ~ "ND",
      median_value >= 5.7 & median_value < 6.5 ~ "pre-T2D",
      median_value >= 6.5 ~ "T2D"
    )
  ) %>%
  ggplot(aes(x = monocle_clusters2, 
             y = median_value, 
             fill = Hba1c_class)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(),
    color = 'black',
    width = 0.8
  ) +
  geom_errorbar(
    aes(ymin = median_value, #- mad_value, 
        ymax = median_value + mad_value),
    width = 0.2,
    position = position_dodge(0.9)
  ) + ylim(0, 9) +
  scale_fill_manual(values = c('#F8766D','#FEAF16FF','#619CFF'))+
  geom_text(
    aes(label = paste0(round(median_value, 1), '%')),
    size = 8, vjust = 2,
    position = position_dodge(0.8)) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median HbA1c (%)")+ 
  guides(fill = guide_legend(title = "HbA1c class"))+
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text =  element_text(face = 'bold', size = 18), 
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

S6F <- plot_cells(cds, show_trajectory_graph = F,
                  color_cells_by = 'Hba1c_class_median',
                  label_cell_groups=FALSE, cell_size = 1,
                  # group_cells_by = 'cluster',
                  group_label_size = 4)+ ggpubr::clean_theme() + NoLegend()+
  scale_color_manual(values = c('#F8766D','#FEAF16FF','#619CFF'))+
  guides(color = guide_legend(title = "HbA1c class"))

(S6E/S6F)

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
### 6. Heatmap Data Prep, Drawing & Legends (Figure S6G)
################################################################################

### Figure S6G

coldata <- data.frame(
  study        = sim$study,
  Cluster      = sim$monocle_clusters2,
  Hba1c_class2 = sim$Hba1c_class_median,
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
# ck1 <- clusterData_custom(exp = mat1_z, scaleData = TRUE, min.std = 0)
# ck2 <- clusterData_custom(exp = mat2_z, scaleData = TRUE, min.std = 0)

# saveRDS(ck1, '20250624_ck1_branched_heatmap.rds')
# saveRDS(ck2, '20250624_ck2_branched_heatmap.rds')

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
  HbA1c    = anno_simple(as.character(coldata[cells1_rev,"Hba1c_class2"]),
                         col = setNames(c("#F8766D","#FEAF16FF","#619CFF"),
                                        sort(unique(coldata[cells1_rev,"Hba1c_class2"])))),
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
  HbA1c    = anno_simple(as.character(coldata[cells2_o,"Hba1c_class2"]),
                         col = setNames(c("#F8766D","#FEAF16FF","#619CFF"),
                                        sort(unique(coldata[cells2_o,"Hba1c_class2"])))),
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

# Legend for HbA1c (shared color scheme)
HbA1c_classes <- sort(unique(c(coldata[cells1_rev, "Hba1c_class2"], coldata[cells2_o, "Hba1c_class2"])))
HbA1c_colors  <- setNames(c("#F8766D", "#FEAF16FF", "#619CFF")[seq_along(HbA1c_classes)], HbA1c_classes)

HbA1c_legend <- Legend(
  labels     = names(HbA1c_colors),
  title      = "HbA1c Class",
  legend_gp  = gpar(fill = HbA1c_colors)
)

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
### Figure S6H
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


### Figure S6I
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
### 7. Module Detection & GO Functions (Figures 4I–4K)
################################################################################

mat1_z <- readRDS('mat1_z.rds')
mat2_z <- readRDS('mat2_z.rds')


perc_expr <- read_tsv('Lineage2_sling_perc_expr.tsv')
colnames(perc_expr)[1:2] <- c("gene", "pct_expr")
dim(perc_expr)


sce <- readRDS('sce_lin2.rds')


coldata <- data.frame(
  study        = sim$study,
  Cluster      = sim$monocle_clusters2,
  Hba1c_class2 = sim$Hba1c_class_median,
  Lineage2_1   = sim$slingPseudotime_1,
  Lineage2_2   = sim$slingPseudotime_2
)
rownames(coldata) <- colnames(sim)

# Get pseudotime
pt1 <- slingPseudotime(sim)[, "Lineage1"]
pt2 <- slingPseudotime(sim)[, "Lineage2"]

# Get Slingshot curve weights
wts <- slingCurveWeights(sim)  # sds is your SlingshotDataSet object

# Cells that are shared between both lineages (i.e., near bifurcation)
shared_cells <- which(wts[, "Lineage1"] > 0.1 & wts[, "Lineage2"] > 0.1)

# Estimate bifurcation pseudotime as the max pseudotime in shared cells
pt_shared <- slingPseudotime(sim)[shared_cells, ]
bifurcation_pt1 <- max(pt_shared[, "Lineage1"], na.rm = TRUE)
bifurcation_pt2 <- max(pt_shared[, "Lineage2"], na.rm = TRUE)

shared_coords <- reducedDim(sim)[shared_cells, ]

pre_branch_cells <- which(pt1 <= bifurcation_pt1 & pt2 <= bifurcation_pt2)

post_branch_l1 <- which(pt1 > bifurcation_pt1)
post_branch_l2 <- which(pt2 > bifurcation_pt2)

# Get weights
pseudotime <- as.data.frame(colData(sce)$crv)

w1 <- pseudotime[, "cellWeights.Lineage1"]
w2 <- pseudotime[, "cellWeights.Lineage2"]

# Estimate bifurcation point: cells with similar weights
shared_cells <- which(abs(w1 - w2) < 0.2 & !is.na(pt1) & !is.na(pt2))
bifurcation_pt <- max(pmin(pt1[shared_cells], pt2[shared_cells]))

cat("Estimated bifurcation pseudotime:", round(bifurcation_pt, 2), "\n")
# Estimated bifurcation pseudotime: 13.59 

### Combine scaled expression matrices used in the module detection
shared_cells <- intersect(colnames(mat1_z), colnames(mat2_z))
mat1_unique <- mat1_z[, setdiff(colnames(mat1_z), shared_cells)]
mat2_unique <- mat2_z[, setdiff(colnames(mat2_z), shared_cells)]
shared_expr <- mat1_z[, shared_cells]
combined_expr_matrix <- cbind(mat1_unique, mat2_unique, shared_expr)

# Shared aging trajectory
cells_pre <- which(pt1 <= bifurcation_pt & pt2 <= bifurcation_pt)

# Sublineage-specific post-bifurcation
cells_post_l1 <- which(pt1 > bifurcation_pt)
cells_post_l2 <- which(pt2 > bifurcation_pt)

# Filter expression matrices
mat_pre <- combined_expr_matrix[, cells_pre]
mat_l1 <- combined_expr_matrix[, cells_post_l1]
mat_l2 <- combined_expr_matrix[, cells_post_l2]

pt_pre <- pt1[cells_pre]
pt_l1 <- pt1[cells_post_l1]
pt_l2 <- pt2[cells_post_l2]

################################################################################
make_elbow_plot <- function(expression_matrix, pseudotime, 
                            sd_threshold = 0.90, nbins = 100, max_pcs = 15) {
  
  # Prepare data (same as your main function)
  row_sds <- matrixStats::rowSds(expression_matrix, na.rm = TRUE)
  names(row_sds) <- rownames(expression_matrix)
  sd_thresh <- quantile(row_sds, sd_threshold, na.rm = TRUE)
  variable_genes <- names(row_sds)[row_sds > sd_thresh]
  
  # FIX: Only keep genes that actually exist in the matrix
  variable_genes <- intersect(variable_genes, rownames(expression_matrix))
  cat("Using", length(variable_genes), "variable genes (top", (1-sd_threshold)*100, "%)\n")
  
  mat_sub <- expression_matrix[variable_genes, ]
  
  # Bin cells by pseudotime
  bin_id <- cut(pseudotime, breaks = nbins, labels = FALSE)
  bin_list <- split(seq_along(bin_id), bin_id)
  bin_list <- bin_list[lengths(bin_list) > 0]
  
  binned_mat <- sapply(bin_list, function(ix) {
    rowMeans(mat_sub[, ix, drop = FALSE], na.rm = TRUE)
  })
  
  # Run PCA
  pca_result <- prcomp(binned_mat, center = TRUE, scale. = TRUE)
  
  # Get variance explained
  variance_explained <- summary(pca_result)$importance[2, ]
  cumulative_var <- summary(pca_result)$importance[3, ]
  
  # Create data frame
  pc_df <- data.frame(
    PC = 1:min(max_pcs, length(variance_explained)),
    VarianceExplained = variance_explained[1:min(max_pcs, length(variance_explained))],
    CumulativeVariance = cumulative_var[1:min(max_pcs, length(variance_explained))]
  )
  
  # Create elbow plot
  p1 <- ggplot(pc_df, aes(x = PC, y = VarianceExplained)) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_point(color = "steelblue", size = 3) +
    labs(title = "Scree Plot - Variance Explained by Each PC",
         x = "Principal Component", 
         y = "Proportion of Variance Explained") +
    theme_minimal() +
    scale_x_continuous(breaks = pc_df$PC) +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 <- ggplot(pc_df, aes(x = PC, y = CumulativeVariance)) +
    geom_line(color = "darkgreen", size = 1.2) +
    geom_point(color = "darkgreen", size = 3) +
    geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed", alpha = 0.6, color = "red") +
    labs(title = "Cumulative Variance Explained",
         x = "Principal Component", 
         y = "Cumulative Proportion") +
    theme_minimal() +
    scale_x_continuous(breaks = pc_df$PC) +
    annotate("text", x = max(pc_df$PC)*0.7, y = c(0.8, 0.9), 
             label = c("80%", "90%"), vjust = -0.5, color = "red") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Print recommendations
  pcs_80 <- which(cumulative_var >= 0.80)[1]
  pcs_90 <- which(cumulative_var >= 0.90)[1]
  
  cat("PCA Variance Summary:\n")
  cat("- For 80% variance:", pcs_80, "PCs\n")
  cat("- For 90% variance:", pcs_90, "PCs\n")
  
  # Fix the elbow calculation to avoid potential errors
  if(length(variance_explained) >= 3) {
    elbow_pc <- which.max(diff(diff(variance_explained[1:min(10, length(variance_explained))]))) + 1
    cat("- Elbow appears around PC", elbow_pc, "\n")
  }
  
  # Combine plots
  combined_plot <- p1 / p2
  
  return(list(
    variance_df = pc_df,
    elbow_plot = combined_plot,
    recommended_pcs_80 = pcs_80,
    recommended_pcs_90 = pcs_90
  ))
}

###

elbow_results_pt_pre <- make_elbow_plot(mat_pre, pt_pre)
elbow_results_pt_pre$elbow_plot

###
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
  
  bin_id <- cut(pseudotime, breaks = nbins, labels = FALSE)
  bin_list <- split(seq_along(bin_id), bin_id)
  bin_list <- bin_list[lengths(bin_list) > 0]
  cat("Created", length(bin_list), "non-empty pseudotime bins\n")
  
  binned_mat <- sapply(bin_list, function(ix) {
    rowMeans(mat_sub[, ix, drop = FALSE], na.rm = TRUE)
  })
  rownames(binned_mat) <- rownames(mat_sub)
  
  pca_result <- prcomp(binned_mat, center = TRUE, scale. = TRUE)
  coords <- pca_result$x[, 1:min(n_pcs, ncol(pca_result$x))]
  
  cat("PCA variance explained by first", ncol(coords), "PCs:", 
      round(sum(summary(pca_result)$importance[2, 1:ncol(coords)]) * 100, 1), "%\n")
  
  dist_mat <- dist(coords)
  hc <- hclust(dist_mat, method = "ward.D2")
  
  modules_raw <- dynamicTreeCut::cutreeDynamic(
    dendro = hc,
    distM = as.matrix(dist_mat),
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
    gene = names(modules_raw),
    module = modules_raw,
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
  
  df_long <- tidyr::pivot_longer(df_plot, cols = -pseudotime, 
                                 names_to = "Module", values_to = "Expression")
  df_long$Module <- paste0(lineage_name, "_M", df_long$Module)
  
  return(list(
    modules = module_list,
    module_df = module_df,
    binned_matrix = binned_mat,
    pca_result = pca_result,
    plot_data = df_long,
    mean_expr_modules = mean_expr_modules
  ))
}
###

results_pre <- identify_pseudotime_modules_pca(
  expression_matrix = mat_pre,
  pseudotime = pt_pre,
  lineage_name = "PreBifurcation",
  nbins = 50,
  min_genes = 200,
  perc_expr = 10,
  sd_threshold = 0.90,
  n_pcs = 3
)

###
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
  return(module_names)
}
###
# Function for GO enrichment analysis
perform_go_enrichment <- function(module_list, lineage_name = "Lineage", 
                                  org_db = org.Hs.eg.db, ont = "BP") {
  
  cat("Performing GO enrichment for", length(module_list), "modules in", lineage_name, "\n")
  
  ego_list <- list()
  
  for (i in seq_along(module_list)) {
    module_name <- names(module_list)[i]
    if(is.null(module_name)) module_name <- i
    
    cat("Processing module", module_name, "with", length(module_list[[i]]), "genes\n")
    
    ego <- enrichGO(
      gene = module_list[[i]],
      OrgDb = org_db,
      keyType = 'SYMBOL',
      ont = ont,
      pvalueCutoff = 0.05,
      pAdjustMethod = 'BH',
      minGSSize = 10,
      maxGSSize = 500
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
  
  return(ego_list)
}
###
ego_pre <- perform_go_enrichment(results_pre$modules, "pre")

###
module_labels <- get_module_names_from_go(ego_pre)
print(module_labels)

###
plot_module_trajectories <- function(results, title = NULL, smooth_method = "loess") {
  
  if(is.null(title)) {
    title <- paste("Module Expression Trajectories")
  }
  
  p <- ggplot(results$plot_data, aes(x = pseudotime, y = Expression, color = Module)) +
    geom_smooth(method = smooth_method, se = FALSE, size = 1.2) +
    labs(
      title = title,
      x = "Pseudotime", 
      y = "Average Expression (log-normalized)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  return(p)
}
###
library(geomtextpath)
plot_pre <- plot_module_trajectories(results_pre, "Pre-branching")
p_pre <- plot_pre +
  ggpubr::theme_pubr(legend = 'right') +
  scale_color_manual(
    labels = c("PreBifurcation_MX1" = "Endocrine Secretion",
               "PreBifurcation_MX2" = "Mitochondria + Cilium",
               "PreBifurcation_MX3" = "Chemotaxis + PLA2",
               "PreBifurcation_MX4" = "MAPK + Fluid Balance"),
    values = c("PreBifurcation_MX1" = "#E31A1C",
               "PreBifurcation_MX2" = "#FF7F00",
               "PreBifurcation_MX3" = "#999999",
               "PreBifurcation_MX4" = "#6A3D9A")) +
  ylab("Average Expression") + 
  ylim(-1.5,1.5)+
  geom_vline(xintercept = round(bifurcation_pt, 2), linetype = "dotted", color = "blue", size = 1.5) +
  geom_textvline(label = "Pre-bifurcation", xintercept = round(bifurcation_pt, 2), vjust = -0.2,
                 size = 5, fontface = 'bold') +
  theme(
    plot.title = element_text(face = 'bold', size = 16, hjust = 0.5),
    axis.title.y = element_text(face = 'bold', size = 14),
    axis.text.x = element_text(face = 'bold', size = 14),
    axis.text.y = element_text(face = 'bold', size = 14),
    axis.title.x = element_text(face = 'bold', size = 15),
    legend.title = element_text(face = 'bold', size = 14),
    legend.text = element_text(face = 'bold', size = 14)
  )

###

results_l1 <- identify_pseudotime_modules_pca(
  expression_matrix = mat_l1,
  pseudotime = pt_l1,
  lineage_name = "Lineage1_Post",
  nbins = 10,
  min_genes = 200,
  perc_expr = 10,
  sd_threshold = 0.90,
  n_pcs = 3
)

ego_lin2_1 <- perform_go_enrichment(results_l1$modules, "post_l2_1")

###
module_labels <- get_module_names_from_go(ego_lin2_1)
print(module_labels)

###

plot_lin2_1 <- plot_module_trajectories(results_l1, "lineage2_1")
p2_1 <- plot_lin2_1+
  ggpubr::theme_pubr(legend = 'right') +
  scale_color_manual(
    labels = c("Lineage1_Post_MX1" = "Secretory morphogenesis +\nepithelial remodeling*",
               "Lineage1_Post_MX2" = "Epithelial Development",
               "Lineage1_Post_MX3" = "MAPK Signaling",
               "Lineage1_Post_MX4" = "Ion Homeostasis + Immune Migration"),
    values = c("Lineage1_Post_MX1" = "#E31A1C",
               "Lineage1_Post_MX2" = "#FF7F00",
               "Lineage1_Post_MX3" = "#6A3D9A",
               "Lineage1_Post_MX4" = "#1F78B4")
  ) + 
  xlim(13.2,18) +
  ylim(-1.5,1.5)+
  geom_vline(xintercept = round(bifurcation_pt, 2), linetype="dotted", 
             color = "blue", size=1.5) + 
  geom_textvline(label = "Post-branching", 
                 xintercept = round(bifurcation_pt, 2), vjust = -0.3,
                 size = 5, fontface = 'bold') +
  ylab(NULL) + 
  theme(plot.title = element_text(face = 'bold', size = 16,hjust = 0.5),
        axis.title.y = element_text(face = 'bold', size = 14),
        axis.text.x =  element_text(face = 'bold', size = 14),
        axis.text.y = element_text(face = 'bold', size = 14),
        axis.title.x = element_text(face = 'bold', size = 15),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14))

###

results_l2 <- identify_pseudotime_modules_pca(
  expression_matrix = mat_l2,
  pseudotime = pt_l2,
  lineage_name = "Lineage2_Post",
  nbins = 10,
  min_genes = 100,
  perc_expr = 10,
  sd_threshold = 0.90,
  n_pcs = 3
)

ego_lin2_2 <- perform_go_enrichment(results_l2$modules, "post_l2_2")

###
module_labels <- get_module_names_from_go(ego_lin2_2)
print(module_labels)

###

plot_lin2_2 <- plot_module_trajectories(results_l2, "lineage2_2")
p2_2 <- plot_lin2_2+
  ylim(-1.5,1.5) +
  ggpubr::theme_pubr(legend = 'right') +
  scale_color_manual(
    labels = c("Lineage2_Post_MX1" = "Proteostasis +\ninnate immune*",
               "Lineage2_Post_MX2" = "Wound healing",
               "Lineage2_Post_MX3" = "ER Stress Response",
               "Lineage2_Post_MX4" = "Cytoplasmic Translation"),
    values = c("Lineage2_Post_MX1" = "#999999",#"#E31A1C",
               "Lineage2_Post_MX2" = "#FF7F00",
               "Lineage2_Post_MX3" = "#1F78B4",
               "Lineage2_Post_MX4" = "#B15928")
  ) + 
  xlim(13.2,18.3) +
  geom_vline(xintercept = round(bifurcation_pt, 2), linetype="dotted", 
             color = "blue", size=1.5) + 
  geom_textvline(label = "Post-branching", 
                 xintercept = round(bifurcation_pt, 2), vjust = -0.3,
                 size = 5, fontface = 'bold') +
  ylab(NULL) + 
  theme(plot.title = element_text(face = 'bold', size = 16, hjust = 0.5),
        axis.title.y = element_text(face = 'bold', size = 14),
        axis.text.x =  element_text(face = 'bold', size = 14),
        axis.text.y = element_text(face = 'bold', size = 14),
        axis.title.x = element_text(face = 'bold', size = 15),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14))

### Figure 4I/J/K
# Cairo::CairoPNG('module_GO_terms.png', width = 23*300,height = 3.5*300,res = 300)
p_pre | p2_1 | p2_2
dev.off()

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
                               round(res_filt_2$padj, 2) <= 0.05),
                            'UP',
                            ifelse((
                              res_filt_2$log2FoldChange < -1 &
                                round(res_filt_2$padj, 2) <= 0.05
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
### 10.                 module scores (Figures 4M-4N)
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
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
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
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
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
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
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
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  theme(
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 14)
  ) +
  NoLegend()

# (p3/p4) | (p1/p2)

# Extract expression data
exp <- lin2_obj@meta.data

# Extract pseudotime values and cell weights
pseudotime <- as.data.frame(colData(sce)$crv)

# Unified pseudotime assignment
pseudotime$UnifiedPseudotime <- rep(NA, nrow(pseudotime))
for (jj in 1:2) {
  lineage_column <- paste0("cellWeights.Lineage", jj)
  pseudotime_column <- paste0("pseudotime.Lineage", jj)
  pseudotime$UnifiedPseudotime[pseudotime[, lineage_column] == 1] <-
    pseudotime[pseudotime[, lineage_column] == 1, pseudotime_column]
}

# Normalize pseudotime
pseudotime$AdjustedPseudotime <- pseudotime$UnifiedPseudotime
for (lineage in 1:2) {
  lineage_mask <- pseudotime[, paste0("cellWeights.Lineage", lineage)] == 1
  min_pseudotime <- min(pseudotime$UnifiedPseudotime[lineage_mask], na.rm = TRUE)
  max_pseudotime <- max(pseudotime$UnifiedPseudotime[lineage_mask], na.rm = TRUE)
  pseudotime$AdjustedPseudotime[lineage_mask] <-
    (pseudotime$UnifiedPseudotime[lineage_mask] - min_pseudotime) /
    (max_pseudotime - min_pseudotime)
}

# Assign lineage membership
pseudotime$Membership <- apply(
  pseudotime[, c("cellWeights.Lineage1", "cellWeights.Lineage2")], 1,
  function(weights) {
    lineages <- c("Lineage2_1", "Lineage2_2")
    paste0(lineages[weights > 0], collapse = ", ")
  }
)

# Assign dominant lineage
pseudotime$DominantLineage <- apply(
  pseudotime[, c("cellWeights.Lineage1", "cellWeights.Lineage2")], 1,
  function(weights) {
    lineages <- c("Lineage2_1", "Lineage2_2")
    if (all(weights == 0)) return(NA)
    lineages[which.max(weights)]
  }
)

pseudotime$Cluster <- lin2_obj@meta.data$monocle_clusters2
pseudotime$Hba1c_class <- lin2_obj@meta.data$Hba1c_class_median

# Colors
point_colors <- c("ND" = "#F8766D", "pre-T2D" = "#FEAF16FF", "T2D" = "#619CFF")
smoother_colors <- c("Lineage2_1" = "#F62424", "Lineage2_2" = "#ffb300")

# Module trends along pseudotime
genes <- c("Core", "Maturity", "ER_Stress", "Exocrine")
l <- c(3, 2, 2, 2.5)

plt_list <- list()
i <- 1
for (gene in genes) {
  plt_data <- data.frame(
    AdjustedPseudotime = pseudotime$AdjustedPseudotime,
    Membership = pseudotime$Membership,
    DominantLineage = pseudotime$DominantLineage,
    Cluster = pseudotime$Cluster,
    HbA1c_class = pseudotime$Hba1c_class,
    gene = exp[, gene]
  )

  p <- ggplot(plt_data, aes(x = AdjustedPseudotime, y = gene)) +
    geom_point(aes(color = HbA1c_class), size = 0.6, alpha = 0.2) +
    geom_smooth(
      data = subset(plt_data, grepl("Lineage2_1", Membership)),
      aes(group = "Lineage2_1"),
      size = 2, color = smoother_colors["Lineage2_1"],
      se = FALSE, method = "lm",
      formula = y ~ splines::ns(x, df = 3)
    ) +
    geom_smooth(
      data = subset(plt_data, grepl("Lineage2_2", Membership)),
      aes(group = "Lineage2_2"),
      size = 2, color = smoother_colors["Lineage2_2"],
      se = FALSE, method = "lm",
      formula = y ~ splines::ns(x, df = 3)
    ) +
    coord_cartesian(ylim = c(-0.5, l[i])) +
    scale_color_manual(
      name = "HbA1c class",
      values = c("ND" = "#F8766D", "pre-T2D" = "#FEAF16FF", "T2D" = "#619CFF")
    ) +
    labs(x = "Pseudotime", y = gene) +
    ggpubr::theme_pubr(legend = "right") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 14),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 14)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  plt_list[[gene]] <- p
  i <- i + 1
}

P <- wrap_plots(plt_list, guides = "collect", ncol = 2)

# Figure 4M-4N composite
(p3 | P[[1]] + NoLegend()) / (p4 | P[[2]] + NoLegend()) |
  (p1 | P[[3]] + NoLegend()) / (p2 | P[[4]] + NoLegend())

################################################################################
### 11. Dot plot C2 vs C5 (Figure S6J)
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
### 12. Dot plot C5 male vs female (Figure S6K)
################################################################################

C5 <- subset(lin2_obj, subset = monocle_clusters2 == "5")
dim(C5)
table(C5$Gender, useNA = "ifany")

C5 <- C5[, C5$Gender %in% c("male", "female") & !is.na(C5$Gender)]
table(C5$Gender, useNA = "ifany")

markers <- list(
  "Exocrine markers" = c("PRSS1", "PRSS2", "CPA1", "CPB1", "CTRB1",
                         "CTRB2", "REG1A", "REG1B", "SPINK1"),
  "EIF2A markers" = c("ESR1", "PPP1R15A", "PPP1R15B", "HSPA5", "DDIT3", "ATF4",
                      "EIF2B2", "EIF2S1", "EIF2AK3")
)

# Figure S6K
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
### 13. Meta-data association with C2 and C5 (Figure 4O)
################################################################################

metadata <- lin2_obj@meta.data[, c(
  "orig.ident", "monocle_clusters2", "Age", "BMI", "Hba1c", "Gender", "Ethnicity"
)]

# Aggregate proportions per donor
cluster_props <- table(metadata$orig.ident, metadata$monocle_clusters2) %>%
  prop.table(margin = 1) %>%
  as.data.frame()
names(cluster_props) <- c("orig.ident", "cluster", "proportion")

# Merge with metadata
donor_data <- merge(
  cluster_props,
  unique(metadata[, c("orig.ident", "Age", "BMI", "Hba1c", "Gender", "Ethnicity")]),
  by = "orig.ident"
)

# Correlations
cor.test(
  donor_data$BMI[donor_data$cluster == "2"],
  donor_data$proportion[donor_data$cluster == "2"]
)
cor.test(
  donor_data$BMI[donor_data$cluster == "5"],
  donor_data$proportion[donor_data$cluster == "5"]
)
cor.test(
  donor_data$Age[donor_data$cluster == "5"],
  donor_data$proportion[donor_data$cluster == "5"]
)
cor.test(
  donor_data$Age[donor_data$cluster == "2"],
  donor_data$proportion[donor_data$cluster == "2"]
)
cor.test(
  donor_data$Hba1c[donor_data$cluster == "5"],
  donor_data$proportion[donor_data$cluster == "5"]
)
cor.test(
  donor_data$Hba1c[donor_data$cluster == "2"],
  donor_data$proportion[donor_data$cluster == "2"]
)

# Gender t tests
t.test(proportion ~ Gender, data = donor_data[donor_data$cluster == "2", ])
t.test(proportion ~ Gender, data = donor_data[donor_data$cluster == "5", ])

# Ethnicity Kruskal Wallis
kruskal.test(proportion ~ Ethnicity, data = donor_data[donor_data$cluster == "5", ])
kruskal.test(proportion ~ Ethnicity, data = donor_data[donor_data$cluster == "2", ])

# Figure 4O table
stats_table <- data.frame(
  Metadata = c(
    "BMI", "BMI", "Age", "Age", "HbA1c", "HbA1c",
    "Gender (male)", "Gender (male)", "Ethnicity", "Ethnicity"
  ),
  Cluster = c("C2", "C5", "C2", "C5", "C2", "C5", "C2", "C5", "C2", "C5"),
  Test = c(
    "Pearson", "Pearson", "Pearson", "Pearson", "Pearson", "Pearson",
    "t-test", "t-test", "Kruskal-Wallis", "Kruskal-Wallis"
  ),
  Statistic = c("0.236", "-0.186", "0.086", "0.074", "0.178", "0.216",
                "t=0.741", "t=-2.346", "chi^2=0.892", "chi^2=1.398"),
  P_value = c("0.0285", "0.0868", "0.4146", "0.4858", "0.1977", "0.1169",
              "0.4609", "0.0219", "0.8273", "0.7060"),
  CI_95 = c(
    "0.026 to 0.427", "-0.383 to 0.027", "-0.121 to 0.286", "-0.133 to 0.274",
    "-0.094 to 0.426", "-0.055 to 0.457", "-0.056 to 0.122", "-0.119 to -0.010",
    "-", "-"
  ),
  stringsAsFactors = FALSE
)

stats_table %>%
  flextable() %>%
  autofit() %>%
  merge_v(j = "Metadata") %>%
  colformat_double(big.mark = ",", digits = 0, j = "Metadata") %>%
  colformat_num(j = c("Statistic", "P_value", "CI_95"), digits = 3) %>%
  bg(i = ~ Cluster == "C2", bg = "#f0f0f0", part = "body") %>%
  bg(i = ~ Cluster == "C5", bg = "#FFFFFF", part = "body") %>%
  bg(bg = "#f0f0f0", part = "header") %>%
  bold(bold = TRUE, part = "header") %>%
  hline(part = "all", border = fp_border(color = "black")) %>%
  align(j = "Metadata", align = "center", part = "body")

###############################################################################
# sessionInfo()
###############################################################################


