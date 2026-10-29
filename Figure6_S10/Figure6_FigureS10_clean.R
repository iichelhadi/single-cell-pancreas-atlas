###############################################################################
# Project: scRNA beta cell atlas — Lineage 3 (Figures 6B–6J, S10A–S10G)
# File: Figure6_FigureS10.R
# Author: Elhadi Iich
# Date: 2026-09-10
#
# Description:
#   Generates Lineage 3 figures and analyses using Monocle3, Slingshot,
#   DESeq2, GSEA and ComplexHeatmap, with revised donor-level analyses
#   and source-data export for the revised manuscript.
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
data.integrated <- readRDS('lin3.beta.reintegrated.rds')
df <- as.data.frame(data.integrated@meta.data)
cds <- readRDS('cds_lin3.rds')

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
  ) +
  scale_fill_manual(values = c('1'='#F8766D','2'='#C49A00','3'='#53B400'))

# Figure 6D
D6 <- df %>%
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
  scale_fill_manual(values = c('1'='#F8766D','2'='#C49A00','3'='#53B400'))

(C6 / D6)

################################################################################
### 4) Supplementary Figure S10A–S10F
################################################################################
# S10A
S10A <- plot_cells(
  cds, show_trajectory_graph = F,
  color_cells_by = 'Study',
  label_cell_groups = FALSE, cell_size = 1, group_label_size = 4
) + ggpubr::clean_theme() + NoLegend() +
  scale_color_manual(values = palette[-c(12, 13, 14)])

# S10B
S10B <- df %>%
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

# S10C
S10C <- plot_cells(
  cds, show_trajectory_graph = F,
  color_cells_by = 'Condition',
  label_cell_groups = FALSE, cell_size = 1, group_label_size = 4
) + ggpubr::clean_theme() + NoLegend()

# S10D
S10D <- df %>%
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

# S10A–D combo
(S10A / S10B) | (S10C / S10D)

################################################################################
# S10E Reconstruct harmonized donor IDs for donor-level HbA1c analysis
################################################################################

annotation_dir <- "revision_out"

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

################################################################################
# Create one metadata row per harmonized donor/sample for S10G
################################################################################

donor_meta <- df %>%
  dplyr::select(
    donor_id,
    Study,
    orig.ident,
    Age,
    BMI,
    Hba1c,
    Gender,
    Ethnicity,
    Condition
  ) %>%
  dplyr::distinct()

donor_meta_check <- donor_meta %>%
  dplyr::count(donor_id) %>%
  dplyr::filter(n > 1)

print(donor_meta_check)

################################################################################
# Figure S10E — donor-level HbA1c distribution
################################################################################

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

S10E <- ggplot(
  hba1c_donor_cluster,
  aes(
    x = factor(monocle_clusters2, levels = desired_order),
    y = Hba1c,
    fill = factor(monocle_clusters2, levels = desired_order)
  )
) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(
    width = 0.15,
    size = 1.2
  ) +
  geom_hline(
    yintercept = 5.7,
    linetype = "dashed"
  ) +
  ggpubr::theme_pubr(legend = "none") +
  labs(
    x = "Cluster",
    y = "Donor HbA1c (%)"
  ) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  )

print(S10E)

hba1c_cluster_summary <- hba1c_donor_cluster %>%
  dplyr::group_by(monocle_clusters2) %>%
  dplyr::summarise(
    n_donor_samples = dplyr::n(),
    median_HbA1c = median(Hba1c, na.rm = TRUE),
    MAD_HbA1c = mad(Hba1c, na.rm = TRUE),
    min_HbA1c = min(Hba1c, na.rm = TRUE),
    max_HbA1c = max(Hba1c, na.rm = TRUE),
    .groups = "drop"
  )

print(hba1c_cluster_summary, n = Inf)

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

clusters <- sort(unique(coldata[cells1_o, "Cluster"]))
cluster_colors <- setNames(c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB")[seq_along(clusters)], clusters)
Cluster_legend <- Legend(labels = names(cluster_colors), title = "Cluster", legend_gp = gpar(fill = cluster_colors))

pt_min <- min(pt1, na.rm = TRUE); pt_max <- max(pt1, na.rm = TRUE)
pseudotime_fun <- colorRamp2(seq(pt_min, pt_max, length.out = 100), plasma(100))
Pseudotime_legend <- Legend(col_fun = pseudotime_fun, title = "Pseudotime", legend_height = unit(3, "cm"),
                            at = round(c(0, 5, 10, 15, 20), 1))

expr_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "#00FFFF", "#FF0000"))
Expression_legend <- Legend(col_fun = expr_fun, title = "Expression (z)", legend_height = unit(2, "cm"))

all_legends <-
  packLegend(Cluster_legend,
             Pseudotime_legend,
             Expression_legend,
             direction = "vertical")

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
  (res_filt_2$log2FoldChange > 1 & res_filt_2$padj <= 0.05), 'UP',
  ifelse((res_filt_2$log2FoldChange < -1 & res_filt_2$padj <= 0.05), 'DOWN', 'NS')
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
### 8) Dotplot (Figure S10F)
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
# 6. Figure S10G
#    Donor-level C3 occupancy
#
# Occupancy =
#   donor's number of Lineage3 C3 cells
#   ----------------------------------
#   donor's total number of Lineage3 cells
#
# Donors with zero C3 cells MUST be retained with occupancy = 0.
################################################################################

donor_totals <- df %>%
  dplyr::filter(
    !is.na(donor_id),
    !is.na(monocle_clusters2)
  ) %>%
  dplyr::count(
    donor_id,
    name = "total_lineage3_cells"
  )

donor_c3_counts <- df %>%
  dplyr::filter(
    !is.na(donor_id),
    monocle_clusters2 == "3"
  ) %>%
  dplyr::count(
    donor_id,
    name = "C3_cells"
  )

donor_props <- donor_totals %>%
  dplyr::left_join(
    donor_c3_counts,
    by = "donor_id"
  ) %>%
  dplyr::mutate(
    C3_cells = tidyr::replace_na(
      C3_cells,
      0L
    ),
    proportion = C3_cells /
      total_lineage3_cells
  ) %>%
  dplyr::left_join(
    donor_meta,
    by = "donor_id"
  ) %>%
  dplyr::arrange(
    donor_id
  )


################################################################################
# 7. Statistical test functions
################################################################################

continuous_test <- function(
    dat,
    variable
) {
  
  d <- dat %>%
    dplyr::filter(
      !is.na(
        .data[[variable]]
      ),
      !is.na(proportion)
    )
  
  if (nrow(d) < 3) {
    
    return(
      tibble::tibble(
        trait = variable,
        test_method = "Spearman correlation",
        n_donors = nrow(d),
        statistic = NA_real_,
        statistic_name = "rho",
        p_value = NA_real_,
        group_summary = "insufficient donors"
      )
    )
  }
  
  fit <- stats::cor.test(
    d[[variable]],
    d$proportion,
    method = "spearman",
    exact = FALSE
  )
  
  tibble::tibble(
    trait = variable,
    test_method = "Spearman correlation",
    n_donors = nrow(d),
    statistic = unname(
      fit$estimate
    ),
    statistic_name = "rho",
    p_value = fit$p.value,
    group_summary = NA_character_
  )
}

sex_test <- function(dat) {
  
  d <- dat %>%
    dplyr::mutate(
      Gender = tolower(
        trimws(Gender)
      )
    ) %>%
    dplyr::filter(
      Gender %in% c(
        "male",
        "female"
      ),
      !is.na(proportion)
    )
  
  if (
    dplyr::n_distinct(
      d$Gender
    ) < 2 ||
    nrow(d) < 3
  ) {
    
    return(
      tibble::tibble(
        trait = "Gender",
        test_method = "Wilcoxon rank-sum",
        n_donors = nrow(d),
        statistic = NA_real_,
        statistic_name = "W",
        p_value = NA_real_,
        group_summary =
          "insufficient donors or single sex"
      )
    )
  }
  
  fit <- stats::wilcox.test(
    proportion ~ Gender,
    data = d,
    exact = FALSE
  )
  
  summary_by_group <- d %>%
    dplyr::group_by(
      Gender
    ) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_proportion = median(
        proportion,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  tibble::tibble(
    trait = "Gender",
    test_method = "Wilcoxon rank-sum",
    n_donors = dplyr::n_distinct(
      d$donor_id
    ),
    statistic = unname(
      fit$statistic
    ),
    statistic_name = "W",
    p_value = fit$p.value,
    group_summary = paste0(
      summary_by_group$Gender,
      ": n=",
      summary_by_group$n,
      ", median=",
      round(
        summary_by_group$median_proportion,
        4
      ),
      collapse = "; "
    )
  )
}

ethnicity_test <- function(dat) {
  
  d <- dat %>%
    dplyr::mutate(
      Ethnicity = trimws(
        Ethnicity
      )
    ) %>%
    dplyr::filter(
      !is.na(Ethnicity),
      Ethnicity != "",
      !is.na(proportion)
    )
  
  if (
    dplyr::n_distinct(
      d$Ethnicity
    ) < 2 ||
    nrow(d) < 3
  ) {
    
    return(
      tibble::tibble(
        trait = "Ethnicity",
        test_method = "Kruskal-Wallis",
        n_donors = nrow(d),
        statistic = NA_real_,
        statistic_name = "chi-square",
        p_value = NA_real_,
        group_summary =
          "insufficient donors or single group"
      )
    )
  }
  
  fit <- stats::kruskal.test(
    proportion ~ Ethnicity,
    data = d
  )
  
  tibble::tibble(
    trait = "Ethnicity",
    test_method = "Kruskal-Wallis",
    n_donors = dplyr::n_distinct(
      d$donor_id
    ),
    statistic = unname(
      fit$statistic
    ),
    statistic_name = "chi-square",
    p_value = fit$p.value,
    group_summary = NA_character_
  )
}

################################################################################
# 8. Run S10G tests
################################################################################

results_S10G <- dplyr::bind_rows(
  continuous_test(
    donor_props,
    "BMI"
  ),
  continuous_test(
    donor_props,
    "Age"
  ),
  continuous_test(
    donor_props,
    "Hba1c"
  ),
  sex_test(
    donor_props
  ),
  ethnicity_test(
    donor_props
  )
) %>%
  dplyr::mutate(
    cluster = "C3",
    .before = 1
  ) %>%
  dplyr::mutate(
    p_adjusted_BH = stats::p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(
    match(
      trait,
      c(
        "BMI",
        "Age",
        "Hba1c",
        "Gender",
        "Ethnicity"
      )
    )
  )

cat("\nFigure S10G results:\n")

print(
  results_S10G,
  n = Inf
)

################################################################################
# 9. Create S10G publication table
################################################################################

stats_table_S10G <- results_S10G %>%
  dplyr::mutate(
    Metadata = dplyr::recode(
      trait,
      BMI = "BMI",
      Age = "Age",
      Hba1c = "HbA1c",
      Gender = "Sex",
      Ethnicity = "Ethnicity"
    ),
    Cluster = cluster,
    Test = test_method,
    Statistic = dplyr::case_when(
      statistic_name == "rho" ~ paste0(
        "rho=",
        sprintf(
          "%.3f",
          statistic
        )
      ),
      
      statistic_name == "W" ~ paste0(
        "W=",
        sprintf(
          "%.0f",
          statistic
        )
      ),
      
      statistic_name ==
        "chi-square" ~ paste0(
          "chi²=",
          sprintf(
            "%.3f",
            statistic
          )
        )
    ),
    P_value = sprintf(
      "%.4f",
      p_value
    ),
    P_adjusted_BH = sprintf(
      "%.4f",
      p_adjusted_BH
    )
  ) %>%
  dplyr::select(
    Metadata,
    Cluster,
    Test,
    n_donors,
    Statistic,
    P_value,
    P_adjusted_BH
  )

print(
  stats_table_S10G,
  n = Inf
)

figS10G_table <- stats_table_S10G %>%
  flextable::flextable() %>%
  flextable::set_header_labels(
    Metadata = "Trait",
    Cluster = "Cluster",
    Test = "Test",
    n_donors = "n donors",
    Statistic = "Statistic",
    P_value = "P",
    P_adjusted_BH = "BH-adjusted P"
  ) %>%
  flextable::autofit() %>%
  flextable::bg(
    i = seq(
      2,
      nrow(
        stats_table_S10G
      ),
      by = 2
    ),
    bg = "#f0f0f0",
    part = "body"
  ) %>%
  flextable::bg(
    bg = "#f0f0f0",
    part = "header"
  ) %>%
  flextable::bold(
    bold = TRUE,
    part = "header"
  ) %>%
  flextable::align(
    align = "center",
    part = "all"
  ) %>%
  flextable::hline(
    part = "all",
    border = officer::fp_border(
      color = "black",
      width = 0.5
    )
  )

figS10G_table

################################################################################
# 10. Save source-data outputs
################################################################################

dir.create(
  "revision_out",
  showWarnings = FALSE
)

readr::write_tsv(
  hba1c_donor_cluster,
  "revision_out/S10E_donor_HbA1c_by_cluster.tsv"
)

readr::write_tsv(
  hba1c_cluster_summary,
  "revision_out/S10E_donor_HbA1c_summary.tsv"
)

readr::write_tsv(
  donor_props,
  "revision_out/S10G_C3_donor_level_occupancy.tsv"
)

readr::write_tsv(
  donor_meta,
  "revision_out/S10G_C3_donor_metadata.tsv"
)

readr::write_tsv(
  results_S10G,
  "revision_out/S10G_C3_donor_trait_associations.tsv"
)

###############################################################################
# SOURCE DATA EXPORT — Figure 6 and Figure S10
###############################################################################

source_dir <- "../SourceData/Figure6_S10"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# Figure 6B, Figure 6E, Figure S10A and Figure S10C
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
      "Lineage3"
    ))
  )

fig6_umap_source <- umap_source %>%
  dplyr::left_join(
    cds_meta_source,
    by = "cell_id"
  )

readr::write_tsv(
  fig6_umap_source,
  file.path(
    source_dir,
    "Fig6B_E_S10A_C_UMAP.tsv.gz"
  )
)


###############################################################################
# Figure 6B — Slingshot trajectory curve
# Reconstruct directly from saved Slingshot object
###############################################################################

sim_source <- readRDS("slingshot_sim_lin3.rds")

curves_source <- slingshot::slingCurves(
  slingshot::SlingshotDataSet(sim_source)
)

fig6_slingshot_source <- lapply(
  names(curves_source),
  function(lineage_name) {
    
    curve_df <- as.data.frame(
      curves_source[[lineage_name]]$s
    )
    
    colnames(curve_df)[1:2] <- c(
      "Dim_1",
      "Dim_2"
    )
    
    curve_df$lineage <- lineage_name
    curve_df$path_id <- paste0(
      lineage_name,
      "_",
      seq_len(nrow(curve_df))
    )
    
    curve_df
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::select(
    lineage,
    path_id,
    Dim_1,
    Dim_2
  )

readr::write_tsv(
  fig6_slingshot_source,
  file.path(
    source_dir,
    "Fig6B_slingshot_curve.tsv"
  )
)

###############################################################################
# Figures 6C and 6D
# Cell-weighted age and BMI distributions
###############################################################################

fig6cd_source <- data.integrated@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    monocle_clusters2,
    Age,
    BMI
  )

readr::write_tsv(
  fig6cd_source,
  file.path(
    source_dir,
    "Fig6C_D_age_BMI.tsv.gz"
  )
)


###############################################################################
# Figures S10B and S10D
# Study and condition composition
###############################################################################

s10_meta <- data.integrated@meta.data %>%
  as.data.frame()

s10b_source <- s10_meta %>%
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

s10d_source <- s10_meta %>%
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
    s10b_source,
    s10d_source
  ),
  file.path(
    source_dir,
    "FigS10B_D_composition.tsv"
  )
)


###############################################################################
# Figure 6F
# INS expression versus Lineage3 pseudotime
###############################################################################

ins_expr <- Seurat::GetAssayData(
  data.integrated,
  assay = "RNA",
  layer = "data"
)["INS", ]

fig6f_source <- data.integrated@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::transmute(
    cell_id,
    pseudotime = Lineage3,
    cluster = monocle_clusters2,
    INS = as.numeric(
      ins_expr[cell_id]
    )
  ) %>%
  dplyr::filter(
    !is.na(pseudotime)
  )

readr::write_tsv(
  fig6f_source,
  file.path(
    source_dir,
    "Fig6F_INS_pseudotime.tsv.gz"
  )
)


###############################################################################
# Figure 6H — GSEA
# Export only the terms actually displayed
###############################################################################

fig6h_ids <- c(
  "GO:0033993",
  "GO:0042127",
  "GO:0007186",
  "GO:0006950",
  "GO:0071396",
  "GO:1901700"
)

fig6h_source <- as.data.frame(
  gse@result
) %>%
  dplyr::filter(
    ID %in% fig6h_ids
  )

readr::write_tsv(
  fig6h_source,
  file.path(
    source_dir,
    "Fig6H_GSEA.tsv"
  )
)


###############################################################################
# Figure 6I — pseudobulk DESeq2 volcano plot
###############################################################################

readr::write_tsv(
  res_filt_2,
  file.path(
    source_dir,
    "Fig6I_pseudobulk_DE.tsv"
  )
)

###############################################################################
# Figure S10F — dot plot
###############################################################################

s10f_genes <- unique(c(
  beta_genes,
  alpha_genes,
  delta_genes,
  progenitor_genes
))

s10f_source <- Seurat::DotPlot(
  data.integrated,
  features = s10f_genes,
  group.by = "monocle_clusters2"
)$data %>%
  as.data.frame() %>%
  dplyr::rename(
    cluster = id,
    gene = features.plot,
    pct_expressing = pct.exp,
    average_expression = avg.exp,
    scaled_average_expression = avg.exp.scaled
  ) %>%
  dplyr::mutate(
    gene_group = dplyr::case_when(
      gene %in% beta_genes ~ "Beta",
      gene %in% alpha_genes ~ "Alpha",
      gene %in% delta_genes ~ "Delta",
      gene %in% progenitor_genes ~ "Progenitor"
    )
  )

readr::write_tsv(
  s10f_source,
  file.path(
    source_dir,
    "FigS10F_dotplot.tsv"
  )
)


###############################################################################
# Figure 6J
# Cell-level module scores plotted against pseudotime
###############################################################################

fig6j_source <- data.integrated@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    Lineage3,
    monocle_clusters2,
    Core,
    Maturity,
    Poly_Endocrine,
    `Lipotoxic Stress`,
    Dedifferentiation,
    Exocrine
  ) %>%
  dplyr::rename(
    pseudotime = Lineage3,
    cluster = monocle_clusters2
  )

readr::write_tsv(
  fig6j_source,
  file.path(
    source_dir,
    "Fig6J_module_scores_pseudotime.tsv.gz"
  )
)


###############################################################################
# Figures S10E and S10G
#
# These are already generated in revision_out.
###############################################################################

existing_source_files <- c(
  "S10E_donor_HbA1c_by_cluster.tsv",
  "S10E_donor_HbA1c_summary.tsv",
  "S10G_C3_donor_level_occupancy.tsv",
  "S10G_C3_donor_metadata.tsv",
  "S10G_C3_donor_trait_associations.tsv"
)

for (f in existing_source_files) {
  
  src <- file.path(
    "revision_out",
    f
  )
  
  if (file.exists(src)) {
    
    file.copy(
      src,
      file.path(
        source_dir,
        f
      ),
      overwrite = TRUE
    )
    
  } else {
    
    warning(
      "Missing existing source-data file: ",
      src
    )
  }
}


###############################################################################
# sessionInfo()
###############################################################################
