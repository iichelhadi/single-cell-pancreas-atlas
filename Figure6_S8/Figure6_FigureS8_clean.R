###############################################################################
# Project: scRNA beta cell atlas — Lineage 3 (Figures 6B–6J, S8A–S8H)
# File: Figure6_FigureS8.R
# Author: Elhadi Iich
# Date: 2025-08-27
#
# Description:
#   Generates Lineage 3 figures and analyses (Monocle3, Slingshot, DESeq2, GSEA,
#   ComplexHeatmap). This version preserves the original code verbatim.
#   Only section banners and comments were added for clarity on GitHub.
#
# Inputs (as used in the script):
#   - lin3.beta.reintegrated.rds
#   - cds_lin3.rds
#   - slingshot_sim_lin3.rds
#   - ck1_branched_heatmap.rds
#   - lin3_pseudotime_start_end_association_filt.tsv
#   - beta_lin3_T2D_pseudobulk_cts.tsv beta_lin3_T2D_pseudobulk_coldata.tsv
#
###############################################################################

################################################################################
### 1) Environment Setup & Library Loading
################################################################################

suppressPackageStartupMessages({
  # Data manipulation and visualization
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)

  # String handling
  library(stringi)
  library(stringr)

  # Single-cell analysis
  library(Seurat)
  library(patchwork)
  library(monocle3)
  library(slingshot)

  # Functional enrichment
  library(clusterProfiler)
  library(org.Hs.eg.db)

  # Table/report generation
  library(flextable)
  library(officer)

  # Additional tools
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(data.table)
  library(viridis)
  library(parallel)
  library(reticulate)
})

set.seed(123)
reticulate::py_module_available(module = 'leidenalg')
reticulate::import('leidenalg')

palette <- c(
  "#5A5156FF", "#E4E1E3FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", "#3283FEFF",
  "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", "#F56642", "#DEA0DFFF",
  "#AA0EDFFF", "#F8A19FFF", "#F62424", "#C4451CFF", "#1C8356FF", "#85660DFF",
  "#B10DA1FF", "#FBE426FF", "#1CBE4FFF"
)

################################################################################
### 2) Data Loading
################################################################################
data.integrated <- readRDS('20250806_lin3.beta.reintegrated.rds')
df <- as.data.frame(data.integrated@meta.data)
cds <- readRDS('20250630_cds_lin3.rds')

desired_order <- c('1', '2', '3')

################################################################################
### 3) Re-clustering (Figures 6B–6D)
################################################################################
# Figure 6B
# B6 <- plot_cells(
#   cds, show_trajectory_graph = F,
#   color_cells_by = 'monocle_clusters2',
#   label_cell_groups = T, cell_size = 1,
#   group_label_size = 10
# ) +
#   ggpubr::clean_theme() +
#   scale_color_manual(values = c('1'='#F8766D','2'='#C49A00','3'='#53B400','4'='#00C094','5'='#00B6EB'))
# B6

# Figure 6C
C6 <- df %>%
  drop_na(Age) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(median_value = median(Age, na.rm = TRUE),
            mad_value = mad(Age, na.rm = TRUE)) %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = monocle_clusters2)) +
  geom_bar(stat = "identity", position = position_dodge(), color = 'black', width = 0.8) +
  geom_errorbar(aes(ymin = median_value, ymax = median_value + mad_value),
                width = 0.2, position = position_dodge(0.9)) +
  ggpubr::theme_pubr(border = FALSE, legend = 'none', base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median Age") +
  theme(axis.title.y = element_text(face = 'bold', size = 14),
        axis.text = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(face = 'bold', size = 14)) +
  scale_fill_manual(values = c('1'='#F8766D','2'='#C49A00','3'='#53B400','4'='#00C094','5'='#00B6EB'))

# Figure 6D
D6 <- df %>%
  drop_na(BMI) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(median_value = median(BMI, na.rm = TRUE),
            mad_value = mad(BMI, na.rm = TRUE)) %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = monocle_clusters2)) +
  geom_bar(stat = "identity", position = position_dodge(), color = 'black', width = 0.8) +
  geom_errorbar(aes(ymin = median_value, ymax = median_value + mad_value),
                width = 0.2, position = position_dodge(0.9)) +
  ggpubr::theme_pubr(border = FALSE, legend = 'none', base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median BMI") +
  theme(axis.title.y = element_text(face = 'bold', size = 14),
        axis.text = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(face = 'bold', size = 14)) +
  scale_fill_manual(values = c('1'='#F8766D','2'='#C49A00','3'='#53B400','4'='#00C094','5'='#00B6EB'))

(C6 / D6)

################################################################################
### 4) Supplementary Figure S8A–S8F
################################################################################
# S8A
S8A <- plot_cells(
  cds, show_trajectory_graph = F,
  color_cells_by = 'Study',
  label_cell_groups = FALSE, cell_size = 1, group_label_size = 4
) + ggpubr::clean_theme() + NoLegend() +
  scale_color_manual(values = palette[-c(12, 13, 14)])

# S8B
S8B <- df %>%
  group_by(monocle_clusters2, Study) %>% dplyr::count() %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  mutate(per = paste0(round(100 * n / sum(n), 1), '%')) %>%
  ggplot(aes(x = monocle_clusters2, y = n, fill = Study)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top', base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  ylab('Proportion') + xlab("Cluster") +
  scale_fill_manual(values = palette[-c(12, 13, 14)]) +
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 18),
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 14))

# S8C
S8C <- plot_cells(
  cds, show_trajectory_graph = F,
  color_cells_by = 'Condition',
  label_cell_groups = FALSE, cell_size = 1, group_label_size = 4
) + ggpubr::clean_theme() + NoLegend()

# S8D
S8D <- df %>%
  group_by(monocle_clusters2, Condition) %>% dplyr::count() %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  mutate(per = paste0(round(100 * n / sum(n), 1), '%')) %>%
  ggplot(aes(x = monocle_clusters2, y = n, fill = Condition)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top', base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  ylab('Proportion') + xlab("Cluster") +
  geom_text(aes(label = per), size = 8, position = position_fill(vjust = 0.5)) +
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 18),
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

# S8A–D combo
(S8A / S8B) | (S8C / S8D)

# S8E
S8E <- df %>%
  drop_na(Hba1c) %>%
  ggplot(aes(monocle_clusters2, Hba1c, fill = monocle_clusters2)) +
  geom_boxplot() +
  ggpubr::theme_pubr(border = FALSE, legend = 'none') +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "HbA1c (%)") +
  scale_fill_manual(values = c('1'='#F8766D','2'='#C49A00','3'='#53B400','4'='#00C094','5'='#00B6EB')) +
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 18),
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

# S8F
S8F <- df %>%
  drop_na(Hba1c) %>%
  mutate(monocle_clusters2 = as.character(monocle_clusters2)) %>%
  group_by(monocle_clusters2) %>%
  summarise(median_value = median(Hba1c, na.rm = TRUE),
            mad_value = mad(Hba1c, na.rm = TRUE)) %>%
  mutate(Hba1c_class = case_when(
    median_value < 5.7 ~ "ND",
    median_value >= 5.7 & median_value < 6.5 ~ "pre-T2D",
    median_value >= 6.5 ~ "T2D"
  )) %>%
  ggplot(aes(x = monocle_clusters2, y = median_value, fill = Hba1c_class)) +
  geom_bar(stat = "identity", position = position_dodge(), color = 'black', width = 0.8) +
  geom_errorbar(aes(ymin = median_value, ymax = median_value + mad_value),
                width = 0.2, position = position_dodge(0.9)) +
  ylim(0, 9) +
  scale_fill_manual(values = c('#F8766D', '#FEAF16FF', '#619CFF')) +
  geom_text(aes(label = paste0(round(median_value, 1), '%')),
            size = 8, vjust = 2, position = position_dodge(0.8)) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top', base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median HbA1c (%)") +
  guides(fill = guide_legend(title = "HbA1c class")) +
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 18),
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

(S8E / S8F)

################################################################################
### 5) Slingshot Curves & UMAP Labels
################################################################################
sim <- readRDS("slingshot_sim_lin3.rds")
lineages <- slingshot::slingCurves(SlingshotDataSet(sim))

lineage_data <- lapply(names(lineages), function(lineage_name) {
  curve_df <- as.data.frame(lineages[[lineage_name]]$s)
  colnames(curve_df) <- c("Dim-1", "Dim-2")
  curve_df$lineage <- lineage_name
  curve_df$path_id <- paste(lineage_name, seq_len(nrow(curve_df)), sep = "_")
  curve_df
})
combined_lineages <- do.call(rbind, lineage_data)

E6 <- plot_cells(
  cds, show_trajectory_graph = FALSE,
  color_cells_by = 'Lineage3',
  cell_size = 1, group_label_size = 10
) + ggpubr::clean_theme() +
  theme(axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 18),
        axis.title.x = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))
E6

# UMAP + labels
umap_coords <- reducedDims(cds)$UMAP
umap_df <- data.frame(
  Dim_1 = umap_coords[, 1],
  Dim_2 = umap_coords[, 2],
  cluster = clusters(cds, reduction_method = "UMAP")
)
label_coords <- aggregate(cbind(Dim_1, Dim_2) ~ cluster, data = umap_df, FUN = mean)
colnames(combined_lineages)[colnames(combined_lineages) == "Dim-1"] <- "Dim_1"
colnames(combined_lineages)[colnames(combined_lineages) == "Dim-2"] <- "Dim_2"

colData(cds)$Lineage3 <- sim@colData$slingPseudotime_1

# “Figure 6B” (trajectory overlay on cluster UMAP)
B6 <- plot_cells(
  cds, show_trajectory_graph = FALSE,
  color_cells_by = 'cluster',
  label_groups_by_cluster = FALSE,
  label_cell_groups = FALSE,
  cell_size = 1, group_label_size = 10
) +
  geom_path(data = combined_lineages,
            aes(x = Dim_1, y = Dim_2, color = lineage, group = lineage),
            size = 2) +
  scale_color_manual(values = c('1'='#F8766D','2'='#C49A00','3'='#53B400','Lineage1'='#E8DD31')) +
  geom_text(data = label_coords, aes(x = Dim_1, y = Dim_2, label = cluster), size = 10)

B6 + ggpubr::clean_theme() +
  theme(legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

# Figure 6F (replotted)
E6 <- plot_cells(
  cds, show_trajectory_graph = FALSE,
  color_cells_by = 'Lineage3',
  cell_size = 1, group_label_size = 10
) + ggpubr::clean_theme() +
  theme(legend.text = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14))

# “Figure 6F” — gene vs pseudotime scatter/smooth
exp3 <- data.integrated@assays$RNA@data
exp3 <- exp3 %>% t()
dim(exp3)  # [1]  8882 18829

meta <- data.integrated@meta.data[, c('Lineage3', 'monocle_clusters2')]
meta <- as.data.frame(meta)
colnames(meta)[1:2] <- c('Lineage3', 'Clusters')
all(rownames(meta) %in% rownames(exp3))

plt_data <- cbind(meta, exp3)
mycolors <- c('#F8766D', '#C49A00', '#53B400')
genes <- 'INS'

plt_list <- list(); i <- 1
for (gene in genes) {
  print(gene); print(i)
  p <- ggscatter(
    data = plt_data, x = 'Lineage3', y = gene,
    color = 'Clusters', size = 0.6
  ) +
    geom_smooth(se = F, color = 'black', method = 'lm',
                formula = 'y~ splines::ns(x, df=3)') +
    ggpubr::theme_pubr() +
    scale_color_manual(values = mycolors) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    xlab('Pseudotime') +
    theme(legend.position = 'right',
          axis.title.y = element_text(face = 'bold', size = 16),
          axis.title.x = element_text(face = 'bold', size = 16),
          axis.text = element_text(face = 'bold', size = 16),
          legend.title = element_text(face = 'bold', size = 14),
          legend.text = element_text(face = 'bold', size = 14, hjust = 1, vjust = 1))
  plt_list[[gene]] <- p; i <- i + 1
}
wrap_plots(plt_list, guides = 'collect')

################################################################################
### 6) Heatmap (6G) & Legends + GSEA (6H)
################################################################################
# (association test code commented out in original)

# 6G heatmap
coldata <- data.frame(
  study = sim$study,
  Cluster = sim$monocle_clusters2,
  Hba1c_class2 = sim$Hba1c_class_median,
  Lineage3 = sim$slingPseudotime_1
)
rownames(coldata) <- colnames(sim)

# External helpers
source("custom_ClusterData.R")
source("custom_visCluster.R")

genes_of_interest <- c('INS','MAFA','ABCC8','IAPP','G6PC2','PAX6',
                       'HSPA1A','HSPA1B','HSPA5','MT1E','JUND',
                       'GCG','CLU','HADH','ASCL2','FFAR4','S100A10')

mat1_raw <- data.integrated@assays$RNA@data
cells1 <- intersect(colnames(mat1_raw), rownames(coldata))
ord1 <- order(coldata[cells1, "Lineage3"])
cells1_o <- cells1[ord1]
pt1 <- coldata[cells1_o, "Lineage3"]
mat1_ord <- mat1_raw[, cells1_o]

ck1 <- readRDS('ck1_branched_heatmap.rds')

mat1_common <- ck1$wide.res %>%
  dplyr::select(-cluster) %>%
  column_to_rownames("gene") %>%
  as.matrix()

pt1_fun <- colorRamp2(seq(min(pt1), max(pt1), length.out = 100), plasma(100))

col_annot1 <- HeatmapAnnotation(
  HbA1c = anno_simple(as.character(coldata[cells1_o, "Hba1c_class2"]),
                      col = setNames(c("#F8766D", "#FEAF16FF", "#619CFF"),
                                     sort(unique(coldata[cells1_o, "Hba1c_class2"])))),
  Cluster = anno_simple(as.character(coldata[cells1_o, "Cluster"]),
                        col = setNames(c("#F8766D", "#C49A00", "#53B400", "#00C094"),
                                       sort(unique(coldata[cells1_o, "Cluster"])))),
  Pseudotime = anno_simple(pt1, col = pt1_fun),
  annotation_height = unit(c(4, 4, 4), "mm"),
  annotation_label = "",
  annotation_name_gp = gpar(fontsize = 12, fontface = 'bold')
)

pos_interest <- which(rownames(mat1_common) %in% genes_of_interest)
row_mark <- rowAnnotation(mark = anno_mark(at = pos_interest,
                                           labels = rownames(mat1_common)[pos_interest],
                                           labels_gp = gpar(fontsize = 14)))

col_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "#00FFFF", "#FF0000"))
ht1 <- Heatmap(
  mat1_common, name = "Expression (z)", col = col_fun,
  top_annotation = col_annot1, right_annotation = row_mark,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  show_row_dend = FALSE, show_heatmap_legend = TRUE, border = TRUE
)

draw(ht1, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()

# Legends
HbA1c_classes <- sort(unique(coldata[cells1_o, "Hba1c_class2"]))
HbA1c_colors <- setNames(c("#F8766D", "#FEAF16FF", "#619CFF")[seq_along(HbA1c_classes)], HbA1c_classes)
HbA1c_legend <- Legend(labels = names(HbA1c_colors), title = "HbA1c Class", legend_gp = gpar(fill = HbA1c_colors))

clusters <- sort(unique(coldata[cells1_o, "Cluster"]))
cluster_colors <- setNames(c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB")[seq_along(clusters)], clusters)
Cluster_legend <- Legend(labels = names(cluster_colors), title = "Cluster", legend_gp = gpar(fill = cluster_colors))

pt_min <- min(pt1, na.rm = TRUE); pt_max <- max(pt1, na.rm = TRUE)
pseudotime_fun <- colorRamp2(seq(pt_min, pt_max, length.out = 100), plasma(100))
Pseudotime_legend <- Legend(col_fun = pseudotime_fun, title = "Pseudotime", legend_height = unit(3, "cm"),
                            at = round(c(0, 5, 10, 15, 20), 1))

expr_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "#00FFFF", "#FF0000"))
Expression_legend <- Legend(col_fun = expr_fun, title = "Expression (z)", legend_height = unit(2, "cm"))

all_legends <- packLegend(HbA1c_legend, Cluster_legend, Pseudotime_legend, Expression_legend, direction = "vertical")
draw(all_legends, just = "left")
dev.off()

# 6H GSEA
pseudotime_start_end_association_filt <- read_tsv(
  'lin3_pseudotime_start_end_association_filt.tsv'
)

genelist <- pseudotime_start_end_association_filt$logFClineage1
names(genelist) <- pseudotime_start_end_association_filt$feature_id
genelist <- genelist[abs(genelist) > 1]
genelist <- genelist[order(genelist, decreasing = T)]
length(genelist)  # 756

library(org.Hs.eg.db)
gse <- clusterProfiler::gseGO(
  geneList = genelist, ont = 'BP',
  OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
  eps = 0, minGSSize = 25, maxGSSize = 500,
  pvalueCutoff = 0.05, pAdjustMethod = 'BH'
)
gse <- clusterProfiler::simplify(gse)

gse@result %>%
  as.data.frame() %>%
  filter(ID %in% c('GO:0033993','GO:0042127','GO:0007186','GO:0006950','GO:0071396','GO:1901700')) %>%
  ggplot(aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
  ylab('Term') + geom_bar(stat = 'identity') + ggpubr::theme_pubr() +
  scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
  theme(axis.title.y = element_text(face = 'bold', size = 14),
        axis.text = element_text(size = 16),
        axis.text.y = element_text(face = 'bold', size = 14),
        axis.title.x = element_text(face = 'bold', size = 15),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14, angle = 45, hjust = 1, vjust = 1))

################################################################################
### 7) Pseudobulk DESeq2 (Figure 6I)
################################################################################
library(DESeq2)

cts <- read.table('beta_lin3_T2D_pseudobulk_cts.tsv',
                  header = T, sep = '\t', row.names = 1)
colnames(cts) <- gsub(colnames(cts), pattern = '[.]', replacement = '-')

coldata <- read.table('beta_lin3_T2D_pseudobulk_coldata.tsv',
                      header = T, sep = '\t', row.names = 1)
coldata$monocle_clusters <- as.character(coldata$monocle_clusters)

perc_expr <- read_tsv('Lineage3_sling_perc_expr.tsv')
dim(perc_expr)
colnames(perc_expr) <- c('feature_id', 'perc_expr')
summary(perc_expr$perc_expr)

cts <- cts[rownames(cts) %in% perc_expr$feature_id[perc_expr$perc_expr >= 5], ]

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ study + monocle_clusters)

size_factors <- coldata$counts / median(coldata$counts)
names(size_factors) <- colnames(dds)
sizeFactors(dds) <- size_factors

dds <- DESeq(dds, parallel = F)
which(mcols(dds)$betaConv) %>% length()
dds <- dds[which(mcols(dds)$betaConv), ]

resultsNames(dds)
res <- results(dds, contrast = c('monocle_clusters', 'C3', 'C2'), alpha = 0.05)
summary(res)

res_filt <- res %>% as.data.frame() %>% filter(complete.cases(.))
res_filt$genes <- rownames(res_filt)
res_filt <- left_join(res_filt, perc_expr, by = c('genes' = 'feature_id'))

res_filt_2 <- res_filt
res_filt_2$change <- ifelse(
  (res_filt_2$log2FoldChange > 1 & round(res_filt_2$padj, 2) <= 0.05), 'UP',
  ifelse((res_filt_2$log2FoldChange < -1 & round(res_filt_2$padj, 2) <= 0.05), 'DOWN', 'NS')
)

gene_sel <- c('TTR','TMED6','FXYD2','ASCL2','TM4SF4','S100A10','SST','PPY','GCG','GC',
              'NKX6-1','INS','IAPP','G6PC2','MAFA','SLC30A8','TMED6','LDHA','LDHB',
              'DLK1','HADH','SAMD11','ABCC8','PDX1','SERPINA1','CRYBA2','CLU',
              'TMEM176A','TMEM176B','TUBB2A','PCSK1','VIM')

ggplot(res_filt_2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(res_filt_2, !genes %in% gene_sel), aes(size = perc_expr), color = "grey") +
  geom_point(data = subset(res_filt_2, genes %in% gene_sel), aes(size = perc_expr, color = change)) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) +
  ggrepel::geom_text_repel(data = subset(res_filt_2, genes %in% gene_sel),
                           aes(label = genes), max.overlaps = 60, box.padding = 0.5, fontface = "bold") +
  guides(color = 'none') + xlim(-6, 6) +
  ggpubr::theme_pubr(legend = 'right') +
  scale_size_continuous(name = "Percent Expression") +
  xlab('Log2 Fold Change') + ylab('Log10 padj') +
  theme(axis.title.y = element_text(face = 'bold', size = 14),
        axis.text = element_text(face = 'bold', size = 14),
        axis.title.x = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14))

################################################################################
### 8) Dotplot (Figure S8G)
################################################################################
beta_genes <- c("INS","PDX1","NKX6-1","PCSK1","MAFA","G6PC2","IAPP","SLC30A8","GAD2","ABCC8")
alpha_genes <- c("GCG","IRX2","ARX","TM4SF4","TTR")
delta_genes <- c("SST","HHEX","GHSR")
progenitor_genes <- c("SOX9","HES1","ALDH1A1","ID2","SOX4")
lipotoxicity_genes <- c("CD36","SREBF1","FASN","MAPK14","CLU","UCP2","ACACA","PDK4","PPARA","EHHADH","FFAR1","FFAR4")

genes_l <- list('Beta'=beta_genes, 'Alpha'=alpha_genes, 'Delta'=delta_genes, 'Progenitor'=progenitor_genes)

SCpubr::do_DotPlot(
  data.integrated, assay = 'RNA', features = genes_l,
  viridis.palette = 'F', viridis.direction = 1, use_viridis = T, flip = F,
  legend.length = 5, group.by = 'monocle_clusters2', scale = T
) +
  scale_y_discrete(limits = c('1','2','3')) +
  ylab('Cluster') +
  theme(axis.title.y = element_text(face = 'bold', size = 14),
        axis.text.x = element_text(face = 'bold', size = 14),
        axis.text.y = element_text(face = 'bold', size = 14),
        axis.title.x = element_text(face = 'bold', size = 14),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14, angle = 0, hjust = 1, vjust = 1))

################################################################################
### 9) Module Analysis (Figure 6J)
################################################################################
core_genes <- c("INS","PDX1","NKX6-1","PCSK1")
maturity_genes <- c("MAFA","G6PC2","IAPP","SLC30A8","ABCC8",'UCN3')
poly_endocrine <- c("GCG","MAFB","TM4SF4","IRX2","SST","HHEX","RBP4","PPY","TTR","VGF","PAX6")
dediff <- c("ALDH1A1","ID2","ID3","SOX4","SOX9","HES1","ONECUT2","KRT19","PROM1","DKK3")
lipotoxic_stress <- c("HADH","CLU","HSP90AA1","HSP90B1","HSPA5","TXNIP","DDIT3","ATF4",
                      "SREBF1","SREBF2","SCD","FADS1","FADS2","PLIN2","APOE","ANGPTL4","GDF15","FKBP11")
exo <- c("PRSS1","PRSS2","REG1A","REG1B","REG3A","SPINK1","CTRB2","CTRB1","CTRC","CPA3")

data.integrated <- AddModuleScore(data.integrated, features = list(Core = core_genes), name = "Core", assay = "RNA")
data.integrated <- AddModuleScore(data.integrated, features = list(maturity = maturity_genes), name = "Maturity", assay = "RNA")
data.integrated <- AddModuleScore(data.integrated, features = list(poly_endocrine = poly_endocrine), name = "Poly_Endocrine", assay = "RNA")
data.integrated <- AddModuleScore(data.integrated, features = list(dediff = dediff), name = "Dedifferentiation", assay = "RNA")
data.integrated <- AddModuleScore(data.integrated, features = list(lipotoxic_stress = lipotoxic_stress), name = "Lipotoxic_stress", assay = "RNA")
data.integrated <- AddModuleScore(data.integrated, features = list(exocrine = exo), name = "Exocrine", assay = "RNA")

colnames(data.integrated@meta.data)[63] <- 'Core'
colnames(data.integrated@meta.data)[64] <- 'Maturity'
colnames(data.integrated@meta.data)[65] <- 'Poly_Endocrine'
colnames(data.integrated@meta.data)[66] <- 'Dedifferentiation'
colnames(data.integrated@meta.data)[67] <- 'Lipotoxic Stress'
colnames(data.integrated@meta.data)[68] <- 'Exocrine'

monocle_theme_opts <- function() {
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25, color = "black")) +
    theme(axis.line.y = element_line(size = 0.25, color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = 'white')) +
    theme(legend.key = element_blank())
}

point_colors <- c("ND"="#F8766D","pre-T2D"="#FEAF16FF","T2D"="#619CFF")
plt_data <- data.integrated@meta.data %>% as.data.frame()
mycolors <- c('#F8766D','#C49A00','#53B400')
genes <- c('Core','Maturity','Poly_Endocrine','Lipotoxic Stress','Dedifferentiation','Exocrine')

plt_list <- list(); i <- 1
for (gene in genes) {
  print(gene); print(i)
  p <- ggscatter(
    data = plt_data, x = 'Lineage3', y = gene,
    color = 'monocle_clusters2', size = 0.6
  ) +
    geom_smooth(se = F, color = 'black', method = 'lm', formula = 'y~ splines::ns(x, df=2)') +
    ggpubr::theme_pubr() +
    scale_color_manual(values = mycolors) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    xlab('Pseudotime') + labs(color = 'Clusters') +
    theme(legend.position = 'right',
          axis.title.y = element_text(face = 'bold', size = 16),
          axis.title.x = element_text(face = 'bold', size = 16),
          axis.text = element_text(face = 'bold', size = 16),
          legend.title = element_text(face = 'bold', size = 14),
          legend.text = element_text(face = 'bold', size = 14, hjust = 1, vjust = 1))
  plt_list[[gene]] <- p; i <- i + 1
}
wrap_plots(plt_list, guides = 'collect', nrow = 2)

################################################################################
### 9) Meta-data correlation (Figure S8H)
################################################################################
metadata <- data.integrated@meta.data[, c("orig.ident","monocle_clusters2","Age","BMI","Hba1c","Gender","Ethnicity")]

cluster_props <- table(metadata$orig.ident, metadata$monocle_clusters2) %>% prop.table(margin = 1) %>% as.data.frame()
names(cluster_props) <- c("orig.ident", "cluster", "proportion")

donor_data <- merge(cluster_props,
                    unique(metadata[, c("orig.ident","Age","BMI","Hba1c","Gender","Ethnicity")]),
                    by = "orig.ident")

cor.test(donor_data$BMI[donor_data$cluster == "3"], donor_data$proportion[donor_data$cluster == "3"])
cor.test(donor_data$Age[donor_data$cluster == "3"], donor_data$proportion[donor_data$cluster == "3"])
cor.test(donor_data$Hba1c[donor_data$cluster == "3"], donor_data$proportion[donor_data$cluster == "3"])
t.test(proportion ~ Gender, data = donor_data[donor_data$cluster == "3", ])
kruskal.test(proportion ~ Ethnicity, data = donor_data[donor_data$cluster == "3", ])

stats_table <- data.frame(
  Metadata = c("BMI","Age","HbA1c","Gender","Ethnicity"),
  Cluster = c("C3","C3","C3","C3","C3"),
  Test = c("Pearson","Pearson","Pearson","t-test","Kruskal-Wallis"),
  Statistic = c("r=0.375","r=0.055","r=0.191","t=-0.500","χ²=7.910"),
  P_value = c("0.002","0.644","0.170","0.619","0.019"),
  CI_95 = c("0.151 to 0.562","-0.176 to 0.280","-0.083 to 0.439","-0.166 to 0.099","-")
)

ft <- stats_table %>%
  flextable() %>%
  merge_v(j = "Metadata") %>%
  colformat_double(x = ., big.mark = ",", digits = 0, j = "Metadata") %>%
  colformat_num(j = c("Statistic","P_value","CI_95"), digits = 3) %>%
  width(j = c("Metadata","Cluster","Test","CI_95"), width = c(1.5,1,1.6,2), unit = "in") %>%
  bg(i = ~ Cluster == "C2", bg = "#f0f0f0", part = "body") %>%
  bg(i = ~ Cluster == "C5", bg = "#FFFFFF", part = "body") %>%
  bg(bg = "#f0f0f0", part = "header") %>%
  bold(bold = TRUE, part = "header") %>%
  hline(part = "all", border = fp_border(color = "black")) %>%
  align(j = "Metadata", align = "center", part = "body")

ft

###############################################################################
# sessionInfo()
###############################################################################
