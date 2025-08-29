################################################################################
# Project: scRNA beta cell atlas — Lineage 2 reintegration and Figure 5 suite
# File: Figure5_FigureS7.R
# Author: Elhadi Iich
# Date: 2025-08-27
#
# Description:
#   Pseudobulk DE (C3 vs C4, C2, C5), volcano plots with highlights, DEG Venn,
#   unified/lineage-specific pseudotime gene trends (S7A–C), HbA1c & lineage
#   stratification UMAPs, SCENIC regulon activity heatmaps and circos, and
#   stress/lipotoxicity module dot plots.
#
# Inputs:
#   - *_cts.tsv, *_coldata.tsv
#   - Lineage2_sling_perc_expr.tsv
#   - DEGsVsC3.RData
#   - lin2.beta.reintegrated.rds
#   - sce_lin2.rds
#   - /beta/lineage2_reintegration_test/SCENIC/Cluster3/int/*Rds
################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(Seurat)
  library(DESeq2)
  library(patchwork)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggpubr)
  library(tibble)
  library(tidyverse)         
  library(ggVennDiagram)
  library(SCpubr)
  library(SCENIC)
  library(AUCell)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(Cairo)
  library(splines)           
})

set.seed(123)

# ---------------------------------------------------------------------------
# Helper: clean volcano plotting (keeps your styling)
plot_volcano_with_highlights <- function(df, gene_sel, xlim=NULL, max_overlaps=40) {
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = subset(df, !genes %in% gene_sel),
               aes(size = perc_expr), color = "grey") +
    geom_point(data = subset(df, genes %in% gene_sel),
               aes(size = perc_expr, color = change)) +
    scale_color_manual(values = c("UP" = "red", "DOWN" = "blue")) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 4) +
    ggrepel::geom_text_repel(
      data = subset(df, genes %in% gene_sel),
      aes(label = genes),
      max.overlaps = max_overlaps,
      box.padding = 0.5,
      fontface = "bold"
    ) +
    guides(color = 'none') +
    ggpubr::theme_pubr(legend = 'right') +
    guides(size = guide_legend(override.aes = list(color = "black"))) +
    scale_size_continuous(name = "Percent Expression") +
    xlab('Log2 Fold Change') + ylab('Log10 padj') +
    theme(axis.title.y = element_text(face = 'bold', size = 14),
          plot.title   = element_text(face = 'bold', size = 14),
          axis.text    = element_text(face = 'bold', size = 14),
          axis.title.x = element_text(face = 'bold', size = 14),
          legend.title = element_text(face = 'bold', size = 14),
          legend.text  = element_text(face = 'bold', size = 14))
  if (!is.null(xlim)) p <- p + xlim(xlim)
  p
}

# =============================================================================
# Section 1 — Pseudobulk DESeq2: C3 vs C4
# =============================================================================
cts <- read.table(
  'beta_C3_v_C4_pseudobulk_cts.tsv',
  header = TRUE, sep = '\t', row.names = 1
)
colnames(cts) <- gsub('[.]','-', colnames(cts))

coldata <- read.table(
  'beta_C3_v_C4_pseudobulk_coldata.tsv',
  header = TRUE, sep = '\t', row.names = 1
)
coldata$monocle_clusters2 <- paste0('C', coldata$monocle_clusters2)

perc_expr <- read_tsv('Lineage2_sling_perc_expr.tsv')
colnames(perc_expr) <- c('feature_id','perc_expr')

cts <- cts[rownames(cts) %in% perc_expr$feature_id[perc_expr$perc_expr >= 5],]

coldata$study <- as.factor(coldata$study)
coldata$monocle_clusters2 <- as.factor(coldata$monocle_clusters2)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ study + monocle_clusters2)

# Size factors proportional to cell counts
sf <- coldata$counts / median(coldata$counts); names(sf) <- colnames(dds); sizeFactors(dds) <- sf

dds <- DESeq(dds, parallel = FALSE)
dds <- dds[which(mcols(dds)$betaConv),]

res <- results(dds, contrast = c('monocle_clusters2', 'C3', 'C4'), alpha = 0.05)

res_filt <- as.data.frame(res) %>%
  filter(complete.cases(.)) %>%
  rownames_to_column('genes') %>%
  left_join(perc_expr, by = c('genes'='feature_id'))


res_filt_2 <- res_filt
res_filt_2$change <- ifelse(res_filt_2$log2FoldChange > 1 & round(res_filt_2$padj, 2) <= 0.05, 'UP',
                            ifelse(res_filt_2$log2FoldChange < -1 & round(res_filt_2$padj, 2) <= 0.05, 'DOWN', ''))

gene_sel <- c('TTR','TMED6','PPP1R1A','FXYD2','ASCL2','TM4SF4',
              'CD82','ELF3','PDLIM2','GRAMD2B','SGIP1','S100A6','NPY','INHBA',
              'IGFBP5','CDKN1A','CDKN2A','SSTR1','ETV1','SIX3','NKX6-1',
              'INS','IAPP','G6PC2','MAFA','SLC30A8','FSTL4')

p_5A <- plot_volcano_with_highlights(res_filt_2, gene_sel, max_overlaps = 20)
print(p_5A)

# =============================================================================
# Section 2 — Pseudobulk DESeq2: C2 vs C3
# =============================================================================
cts <- read.table(
  'beta_C2_v_C3_pseudobulk_cts.tsv',
  header = TRUE, sep = '\t', row.names = 1
)
colnames(cts) <- gsub('[.]','-', colnames(cts))

coldata <- read.table(
  'beta_C2_v_C3_pseudobulk_coldata.tsv',
  header = TRUE, sep = '\t', row.names = 1
)
coldata$monocle_clusters2 <- paste0('C', coldata$monocle_clusters2)

perc_expr <- read_tsv('Lineage2_sling_perc_expr.tsv')
colnames(perc_expr) <- c('feature_id','perc_expr')

cts <- cts[rownames(cts) %in% perc_expr$feature_id[perc_expr$perc_expr >= 5],]

coldata$study <- as.factor(coldata$study)
coldata$monocle_clusters2 <- as.factor(coldata$monocle_clusters2)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ study + monocle_clusters2)

sf <- coldata$counts / median(coldata$counts); names(sf) <- colnames(dds); sizeFactors(dds) <- sf

dds <- DESeq(dds, parallel = FALSE)
dds <- dds[which(mcols(dds)$betaConv),]

res <- results(dds, contrast = c('monocle_clusters2', 'C2', 'C3'), alpha = 0.05)

res_filt <- as.data.frame(res) %>%
  filter(complete.cases(.)) %>%
  rownames_to_column('genes') %>%
  left_join(perc_expr, by = c('genes'='feature_id'))


res_filt_2 <- res_filt
res_filt_2$change <- ifelse(res_filt_2$log2FoldChange > 1 & round(res_filt_2$padj, 2) <= 0.05, 'UP',
                            ifelse(res_filt_2$log2FoldChange < -1 & round(res_filt_2$padj, 2) <= 0.05, 'DOWN', ''))

gene_sel <- c('NKX6-1','SLC30A8','NEUROD1','DNAJC15',
              'INS','IAPP','G6PC2','MAFA','PAX6',
              'ABCC8','PCSK1','GAD2','FXYD2','LMO1',
              'ATF3','DDIT3','PPP1R15A','SLC2A1','YOD1')

p_5B <- plot_volcano_with_highlights(res_filt_2, gene_sel, max_overlaps = 60)
print(p_5B)

# =============================================================================
# Section 3 — Pseudobulk DESeq2: C3 vs C5
# =============================================================================
cts <- read.table(
  'beta_C3_v_C5_pseudobulk_cts.tsv',
  header = TRUE, sep = '\t', row.names = 1
)
colnames(cts) <- gsub('[.]','-', colnames(cts))

coldata <- read.table(
  'beta_C3_v_C5_pseudobulk_coldata.tsv',
  header = TRUE, sep = '\t', row.names = 1
)
coldata$monocle_clusters2 <- paste0('C', coldata$monocle_clusters2)

perc_expr <- read_tsv('Lineage2_sling_perc_expr.tsv')
colnames(perc_expr) <- c('feature_id','perc_expr')

cts <- cts[rownames(cts) %in% perc_expr$feature_id[perc_expr$perc_expr >= 5],]

coldata$study <- as.factor(coldata$study)
coldata$monocle_clusters2 <- as.factor(coldata$monocle_clusters2)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ study + monocle_clusters2)

sf <- coldata$counts / median(coldata$counts); names(sf) <- colnames(dds); sizeFactors(dds) <- sf

dds <- DESeq(dds, parallel = FALSE)
dds <- dds[which(mcols(dds)$betaConv),]

res <- results(dds, contrast = c('monocle_clusters2', 'C5', 'C3'), alpha = 0.05)

res_filt <- as.data.frame(res) %>%
  filter(complete.cases(.)) %>%
  rownames_to_column('genes') %>%
  left_join(perc_expr, by = c('genes'='feature_id'))

res_filt_2 <- res_filt
res_filt_2$change <- ifelse(res_filt_2$log2FoldChange > 1 & round(res_filt_2$padj, 2) <= 0.05, 'UP',
                            ifelse(res_filt_2$log2FoldChange < -1 & round(res_filt_2$padj, 2) <= 0.05, 'DOWN', ''))

gene_sel <- c('NKX6-1','SLC30A8','NEUROD1','PDK4',
              'INS','IAPP','G6PC2','MAFA','PAX6',
              'ABCC8','PCSK1','GAD2','FXYD2','GLP1R',
              'REG1A','SPINK1','PRSS1','CTRB2','REG1B',
              'CPA1','S100A10','TM4SF4')

p_5C <- plot_volcano_with_highlights(res_filt_2, gene_sel, max_overlaps = 60)
print(p_5C)

# =============================================================================
# Section 4 — Figure 5D: DEG Venn (C4vC3, C3vC2, C3vC5)
# =============================================================================
load('DEGsVsC3.RData')
# Expecting objects: C4vC3, C3vC2, C3vC5 with $genes
x <- list(
  C4vC3 = C4vC3$genes,
  C3vC2 = C3vC2$genes,
  C3vC5 = C3vC5$genes
)

p_5D <- ggVennDiagram(
  x,
  set_color = "black",
  lwd = 0.8,
  label_size = 6,
  set_size = 8,
  lty = 1,
  label = "both"
) + scale_fill_gradient(low = "#FFFFFF", high = "#FEAF16FF") +
  scale_color_manual(values = c("black", "black","black")) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))+  # Adds spacing on x-axis
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +  # Adds spacing on y-axis
  theme(legend.title = element_text(size = 16, face = "bold"),
        legend.text  = element_text(size = 14))
print(p_5D)

# =============================================================================
# Section 5 — Figure 5E & S7A–C: Pseudotime gene trends & lineage overlays
# =============================================================================
lin2_obj_full <- readRDS("lin2.beta.reintegrated.rds")
exp <- t(lin2_obj_full@assays$RNA@data)

# Load SCE for full pseudotime (for Fig 5E)
sce <- readRDS('sce_lin2.rds')

# Unified pseudotime (Lineage2_1 vs Lineage2_2)
pseudotime <- as.data.frame(colData(sce)$crv)
pseudotime$UnifiedPseudotime <- NA_real_
for (jj in 1:2) {
  lc <- paste0("cellWeights.Lineage", jj)
  pc <- paste0("pseudotime.Lineage", jj)
  mask <- pseudotime[[lc]] == 1
  pseudotime$UnifiedPseudotime[mask] <- pseudotime[[pc]][mask]
}

pseudotime$AdjustedPseudotime <- pseudotime$UnifiedPseudotime
for (lineage in 1:2) {
  lc <- paste0("cellWeights.Lineage", lineage)
  mask <- pseudotime[[lc]] == 1
  min_pt <- min(pseudotime$UnifiedPseudotime[mask], na.rm = TRUE)
  max_pt <- max(pseudotime$UnifiedPseudotime[mask], na.rm = TRUE)
  pseudotime$AdjustedPseudotime[mask] <- (pseudotime$UnifiedPseudotime[mask] - min_pt) / (max_pt - min_pt)
}

pseudotime$Membership <- apply(pseudotime[, c("cellWeights.Lineage1","cellWeights.Lineage2")], 1, function(w) {
  lineages <- c("Lineage2_1","Lineage2_2")
  paste0(lineages[w > 0], collapse = ", ")
})

pseudotime$DominantLineage <- apply(pseudotime[, c("cellWeights.Lineage1","cellWeights.Lineage2")], 1, function(w) {
  lineages <- c("Lineage2_1","Lineage2_2"); if (all(w == 0)) return(NA_character_); lineages[which.max(w)]
})

pseudotime$Cluster     <- lin2_obj_full@meta.data$monocle_clusters2
pseudotime$HbA1c_class <- lin2_obj_full@meta.data$Hba1c_class_median

# ------------ S7 gene panels (A/B/C), you can switch which block is active
# # S7A
# genes <- c('INS','TMED6','TTR','PPP1R1A','CD82','FXYD2','ASCL2'); l <- c(9,2,6,4,2,5,2)
# # S7B
genes <- c('S100A6','IGFBP5','FSTL4','SGIP1','INHBA','SIX3','NPY','CDKN1A','CDKN2A')
l     <- c(5,4.5,2.5,2,4,3,6,5,3)
# # S7C
# genes <- c('SLC30A8','GAD2','NEUROD1','MAFA','ABCC8','DNAJC15','PCSK1','LMO1')
# l     <- c(4,5,4,4,5,3,4,2.5)

smoother_colors <- c("Lineage2_1" = "#F62424", "Lineage2_2" = "#f56642")
hba1c_colors    <- c("ND" = "#F8766D", "pre-T2D" = "#FEAF16FF", "T2D" = "#619CFF")

plt_list <- list(); i <- 1
for (gene in genes) {
  plt_data <- data.frame(
    AdjustedPseudotime = pseudotime$AdjustedPseudotime,
    Membership         = pseudotime$Membership,
    DominantLineage    = pseudotime$DominantLineage,
    Cluster            = pseudotime$Cluster,
    HbA1c_class        = pseudotime$HbA1c_class,
    gene               = exp[, gene]
  )
  p <- ggplot(plt_data, aes(x = AdjustedPseudotime, y = gene)) +
    geom_point(aes(color = HbA1c_class), size = 0.6, alpha = 0.2) +
    geom_smooth(data = subset(plt_data, grepl("Lineage2_1", Membership)),
                aes(group = "Lineage2_1"), size = 2, color = smoother_colors["Lineage2_1"],
                se = FALSE, method = 'lm', formula = y ~ splines::ns(x, df = 3)) +
    geom_smooth(data = subset(plt_data, grepl("Lineage2_2", Membership)),
                aes(group = "Lineage2_2"), size = 2, color = smoother_colors["Lineage2_2"],
                se = FALSE, method = 'lm', formula = y ~ splines::ns(x, df = 3)) +
    coord_cartesian(ylim = c(0, l[i])) +
    scale_color_manual(name = "HbA1c class", values = hba1c_colors) +
    labs(x = "Pseudotime", y = gene) +
    ggpubr::theme_pubr(legend = 'right') +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(face = 'bold', size = 14),
          plot.title   = element_text(face = 'bold', size = 14),
          axis.text    = element_text(face = 'bold', size = 14),
          axis.title.x = element_text(face = 'bold', size = 14),
          legend.title = element_text(face = 'bold', size = 14),
          legend.text  = element_text(face = 'bold', size = 14)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  plt_list[[gene]] <- p
  i <- i + 1
}
print(wrap_plots(plt_list, guides = 'collect', ncol = 4))

# =============================================================================
# Section 6 — Figures 5F & 5G: DotPlots by HbA1c, plus lineage overlays in C3
# =============================================================================
# HbA1c-class dot plot
# p_dot <- SCpubr::do_DotPlot(
#   lin2_obj_full,
#   features = c('PRG4','FSTL4','IGFBP5','CDKN2A'),
#   assay = 'RNA',
#   group.by = 'Hba1c_class_median',
#   use_viridis = TRUE,
#   legend.length = 3,
#   legend.position = 'right',
#   viridis.direction = 1,
#   font.size = 16,
#   ylab = 'HbA1c\n class',
#   viridis.palette = 'F',
#   scale = TRUE
# ) +
#   theme(axis.title.y = element_text(face = 'bold', size = 14),
#         axis.text.x  = element_text(face = 'bold', size = 14),
#         axis.text.y  = element_text(face = 'bold', size = 14),
#         axis.title.x = element_text(face = 'bold', size = 14),
#         legend.title = element_text(face = 'bold', size = 14),
#         legend.text  = element_text(face = 'bold', size = 14))
# print(p_dot)

# C3-only pseudotime + lineage grouping for UMAP overlays
colData(sce)$Clusters <- lin2_obj_full@meta.data$monocle_clusters2
lin2_obj_c3 <- subset(lin2_obj_full, subset = monocle_clusters2 == '3')

sce_c3 <- sce[, colData(sce)$Clusters == "3"]
pseudotime_c3 <- as.data.frame(colData(sce_c3)$crv)
pseudotime_c3$Clusters <- colData(sce_c3)$Clusters

# Dominant-lineage unified pseudotime in C3
pseudotime_c3$UnifiedPseudotime <- NA_real_
for (jj in 1:2) {
  lc <- paste0("cellWeights.Lineage", jj)
  pc <- paste0("pseudotime.Lineage", jj)
  mask <- pseudotime_c3[[lc]] > 0.5
  pseudotime_c3$UnifiedPseudotime[mask] <- pseudotime_c3[[pc]][mask]
}

pseudotime_c3$AdjustedPseudotime <- NA_real_
for (lineage in 1:2) {
  lc <- paste0("cellWeights.Lineage", lineage)
  mask <- pseudotime_c3[[lc]] > 0.5
  min_pt <- min(pseudotime_c3$UnifiedPseudotime[mask], na.rm = TRUE)
  max_pt <- max(pseudotime_c3$UnifiedPseudotime[mask], na.rm = TRUE)
  pseudotime_c3$AdjustedPseudotime[mask] <- (pseudotime_c3$UnifiedPseudotime[mask] - min_pt) / (max_pt - min_pt)
}
# exact-1 fallback (kept from your code)
for (lineage in 1:2) {
  lc <- paste0("cellWeights.Lineage", lineage)
  mask <- pseudotime_c3[[lc]] == 1
  min_pt <- min(pseudotime_c3$UnifiedPseudotime[mask], na.rm = TRUE)
  max_pt <- max(pseudotime_c3$UnifiedPseudotime[mask], na.rm = TRUE)
  pseudotime_c3$AdjustedPseudotime[mask] <- (pseudotime_c3$UnifiedPseudotime[mask] - min_pt) / (max_pt - min_pt)
}

pseudotime_c3$Membership <- apply(pseudotime_c3[, c("cellWeights.Lineage1","cellWeights.Lineage2")], 1, function(w) {
  lineages <- c("Lineage2_1","Lineage2_2")
  paste0(lineages[w > 0], collapse = ", ")
})
pseudotime_c3$DominantLineage <- apply(pseudotime_c3[, c("cellWeights.Lineage1","cellWeights.Lineage2")], 1, function(w) {
  lineages <- c("Lineage2_1","Lineage2_2"); if (all(w == 0)) return(NA_character_); lineages[which.max(w)]
})

# LineageGroup thresholding
pseudotime_c3$LineageGroup <- NA_character_
pseudotime_c3$LineageGroup[abs(pseudotime_c3$cellWeights.Lineage1 - pseudotime_c3$cellWeights.Lineage2) < 0.35] <- "Transitional"
pseudotime_c3$LineageGroup[pseudotime_c3$cellWeights.Lineage1 > pseudotime_c3$cellWeights.Lineage2 + 0.35] <- "Lineage2_1"
pseudotime_c3$LineageGroup[pseudotime_c3$cellWeights.Lineage2 > pseudotime_c3$cellWeights.Lineage1 + 0.35] <- "Lineage2_2"

# Add C3 metadata to the C3 object only (mapping via rownames)
# Ensure rownames present
if (is.null(rownames(pseudotime_c3))) rownames(pseudotime_c3) <- colnames(sce_c3)
lin2_obj_c3 <- AddMetaData(lin2_obj_c3, metadata = pseudotime_c3)

p1 <- DimPlot(lin2_obj_c3, reduction = 'umap_MNN', group.by = 'LineageGroup') +
  scale_color_manual(values = c("Lineage2_1" = "#00B0F6",
                                "Lineage2_2" = "#E76BF3",
                                "Transitional" = "#9590FF")) +
  ggpubr::theme_pubr(legend = 'bottom') +
  ggpubr::font("title", size = 16, color = "black", face = "bold") +
  NoAxes() + NoGrid() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.text = element_text(face = 'bold', size = 14))

lin2_obj_c3$Membership2 <- ifelse(lin2_obj_c3$Membership == "Lineage2_1, Lineage2_2", "Both", lin2_obj_c3$Membership)

p2 <- DimPlot(lin2_obj_c3, reduction = 'umap_MNN', group.by = 'Membership2') +
  scale_color_manual(values = c("Lineage2_1" = "#00B0F6",
                                "Lineage2_2" = "#E76BF3",
                                "Both"       = "#9590FF")) +
  ggpubr::theme_pubr(legend = 'bottom') +
  ggpubr::font("title", size = 16, color = "black", face = "bold") +
  NoAxes() + NoGrid() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.text = element_text(face = 'bold', size = 14))

print(p2 / p1)

# =============================================================================
# Section 7 — SCENIC regulon activity heatmaps + circos (Cluster 3)
# =============================================================================
regulonAUC <- readRDS('3.4_regulonAUC.Rds')
cellInfo   <- readRDS('cellInfo.Rds')

regulonAUC@colData$Hba1c_class   <- cellInfo$Hba1c_class_median
regulonAUC@colData$Condition     <- cellInfo$Condition
regulonAUC@colData$Membership    <- as.character(cellInfo$Membership)
regulonAUC@colData$DominantLineage <- as.character(cellInfo$DominantLineage)
regulonAUC@colData$LineageGroup  <- as.character(cellInfo$LineageGroup)

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byLineageGroup <- sapply(
  split(rownames(cellInfo), regulonAUC@colData$LineageGroup),
  function(cells) rowMeans(getAUC(regulonAUC)[, cells])
)
regulonActivity_byLineageGroup_Scaled <- t(scale(t(regulonActivity_byLineageGroup), center = TRUE, scale = TRUE))
df <- regulonActivity_byLineageGroup_Scaled[grep('extended', rownames(regulonActivity_byLineageGroup_Scaled), invert = TRUE), ]

ht <- ComplexHeatmap::Heatmap(
  df, name = "Regulon activity",
  cluster_rows = TRUE, cluster_columns = FALSE, column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_order = c('Lineage2_1','Transitional','Lineage2_2')
)
draw(ht)

# Circular heatmap (circos)
mat <- as.matrix(df[, c("Lineage2_1","Transitional","Lineage2_2")])
row_dend <- hclust(dist(mat)); mat <- mat[row_dend$order, ]
col_fun <- colorRamp2(c(min(mat), 0, max(mat)), c("blue","white","red"))

# Cairo::CairoPNG('20250627_circos.png', width = 8, height = 8, units = "in", dpi = 300)
circos.clear()
circos.par(start.degree = 90, gap.degree = 40, cell.padding = c(0,0,0,0))

l1    <- c("ATF4 (28g)", "ATF3 (36g)", "FOS (22g)", "JUN (21g)","HNF4A (12g)")
l2    <- c("CREB5 (27g)", "ETV4 (16g)", "NR1D1 (10g)", "PPARD (11g)", "NFKB2 (18g)")
trans <- c("HNF1A (19g)", "KLF6 (27g)", "TP53 (18g)", "EGR3 (15g)", "EGR4 (20g)")

all_genes   <- rownames(mat)
color_vector <- ifelse(all_genes %in% l1, "#00B0F6",
                       ifelse(all_genes %in% trans, "#9590FF",
                              ifelse(all_genes %in% l2, "#E76BF3", 'black')))
font_vector  <- ifelse(all_genes %in% c(l1,l2,trans), 2, 1)

circos.heatmap(
  mat,
  rownames.font = font_vector,
  rownames.col  = color_vector,
  col           = col_fun,
  track.height  = 0.3, rownames.cex = 0.8,
  rownames.side = "inside",
  dend.side     = "outside",
  cell.border   = "white",
  cell.lwd      = 1.5
)

ht_track <- get.current.track.index()
heatmap_sector <- get.all.sector.index()[1]
n1 <- 5; n2 <- 3; n3 <- 5
circos.trackPlotRegion(
  heatmap_sector, track.index = ht_track + 1, ylim = c(0,1), track.height = 0.1,
  panel.fun = function(x, y) {
    circos.text(x = n1/4, y = 1.5, labels = "Lineage2_1", facing = "downward", niceFacing = TRUE,
                adj = c(1.4, -20), cex = 0.8, font = 2)
    circos.text(x = n2/4, y = 1.5, labels = "Transitional", facing = "downward", niceFacing = TRUE,
                adj = c(1.3, -16), cex = 0.8, font = 2)
    circos.text(x = n3/4, y = 1.5, labels = "Lineage2_2", facing = "downward", niceFacing = TRUE,
                adj = c(1.4, -12.3), cex = 0.8, font = 2)
  },
  bg.border = NA
)
dev.off()

lgd <- Legend(title = "Regulon activity", col_fun = col_fun)
draw(lgd, x = unit(0.9,"npc"), y = unit(0.1,"npc"), just = c("right","bottom"))


###############################################################################
# sessionInfo()
###############################################################################
