################################################################################
# Project: scRNA beta cell atlas — Lineage 2 reintegration and Figure 5 suite
# File: Figure5_FigureS8.R
# Author: Elhadi Iich
# Date: 2026-09-10
#
# Description:
#   Pseudobulk DE (C3 vs C4, C2, C5), volcano plots with highlights, DEG Venn,
#   unified/lineage-specific pseudotime gene trends (S8A–C), HbA1c & lineage
#   stratification UMAPs, SCENIC regulon activity heatmaps and circos, and
#   stress/lipotoxicity module dot plots.
#
# Inputs:
#   - *_cts.tsv, *_coldata.tsv
#   - Lineage2_sling_perc_expr.tsv
#   - DEGsVsC3.RData
#   - lin2.beta.reintegrated.rds
#   - sce_lin2.rds
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
  dplyr::filter(complete.cases(.)) %>%
  rownames_to_column('genes') %>%
  left_join(perc_expr, by = c('genes'='feature_id'))


res_filt_2 <- res_filt

res_filt_2$change <- ifelse(res_filt_2$log2FoldChange > 1 & res_filt_2$padj <= 0.05, 'UP',
                            ifelse(res_filt_2$log2FoldChange < -1 & res_filt_2$padj <= 0.05, 'DOWN', ''))

fig5A_source <- res_filt_2

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
  dplyr::filter(complete.cases(.)) %>%
  rownames_to_column('genes') %>%
  left_join(perc_expr, by = c('genes'='feature_id'))


res_filt_2 <- res_filt
res_filt_2$change <- ifelse(res_filt_2$log2FoldChange > 1 & res_filt_2$padj <= 0.05, 'UP',
                            ifelse(res_filt_2$log2FoldChange < -1 & res_filt_2$padj <= 0.05, 'DOWN', ''))

fig5B_source <- res_filt_2

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
res_filt_2$change <- ifelse(res_filt_2$log2FoldChange > 1 & res_filt_2$padj <= 0.05, 'UP',
                            ifelse(res_filt_2$log2FoldChange < -1 & res_filt_2$padj <= 0.05, 'DOWN', ''))

fig5C_source <- res_filt_2

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

## Fig 5D: up in C3 vs C4, down in C2 vs C3, down in C5 vs C3 (padj <= 0.05, |log2FC| > 1)
C4vC3 <- fig5A_source %>% dplyr::filter(change == 'UP')
C3vC2 <- fig5B_source %>% dplyr::filter(change == 'DOWN')
C3vC5 <- fig5C_source %>% dplyr::filter(change == 'DOWN')
save(C4vC3, C3vC2, C3vC5, file = 'DEGsVsC3.RData')
Reduce(intersect, list(C4vC3$genes, C3vC2$genes, C3vC5$genes))   # expect CDKN2A, FSTL4, IGFBP5, PRG4

load('DEGsVsC3.RData')
# Expecting objects: C4vC3, C3vC2, C3vC5 with $genes
x <- list(
  C4vC3 = C4vC3$genes,
  C3vC2 = C3vC2$genes,
  C3vC5 = C3vC5$genes
)

p_5D <- ggVennDiagram::ggVennDiagram(
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

################################################################################
### 12. Unified Pseudotime Trajectory Plots (Figure 5E and Figure S8A–S8C)
################################################################################

# Extract expression data and pseudotime
lin2_obj_full <- readRDS("lin2.beta.reintegrated.rds")

exp <- lin2_obj_full@assays$RNA@data %>% t()
gc()

sce <- readRDS("sce_lin2.rds")

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

colnames(pt) <- c(
  "Lineage2_1",
  "Lineage2_2"
)

colnames(wt) <- c(
  "Lineage2_1",
  "Lineage2_2"
)

# Cluster identity for each cell
cluster <- setNames(
  paste0(
    "C",
    as.character(lin2_obj_full$monocle_clusters2)
  ),
  colnames(lin2_obj_full)
)


cluster_colors <- c(
  "C1" = "#F8766D",
  "C2" = "#C49A00",
  "C3" = "#53B400",
  "C4" = "#00C094",
  "C5" = "#00B6EB"
)

smoother_colors <- c(
  "Lineage2_1" = "#F62424",
  "Lineage2_2" = "#f56642"
)


################################################################################
### Figure 5E
################################################################################

# genes <- c(
#   "PRG4",
#   "FSTL4",
#   "IGFBP5",
#   "CDKN2A"
# )
# l <- c(2,2,3,2)

################################################################################
### Figure S8A
################################################################################

# genes <- c(
#   "INS",
#   "TMED6",
#   "TTR",
#   "PPP1R1A",
#   "CD82",
#   "FXYD2",
#   "ASCL2"
# )
# l <- c(8.5,2,6,4,2,5,2)

################################################################################
### Figure S8B
################################################################################

# genes <- c(
#   "S100A6",
#   "IGFBP5",
#   "FSTL4",
#   "SGIP1",
#   "INHBA",
#   "SIX3",
#   "NPY",
#   "CDKN1A",
#   "CDKN2A"
# )
# l <- c(5,4,2.5,2,3,3,5,4,3)

################################################################################
### Figure S8C
################################################################################

genes <- c(
  "SLC30A8",
  "GAD2",
  "NEUROD1",
  "MAFA",
  "ABCC8",
  "DNAJC15",
  "PCSK1",
  "LMO1"
)
l <- c(4,5,4,4,5,3,4,2.5)


################################################################################
### Plot
################################################################################

plt_list <- list()
i <- 1
for (i in seq_along(genes)) {
  
  gene <- genes[i]
  
  gene_expr <- exp[, gene]
  
  plot_data <- lapply(
    seq_len(ncol(pt)),
    function(j) {
      
      keep <- !is.na(pt[, j]) & wt[, j] > 0
      
      cells_keep <- rownames(pt)[keep]
      
      tibble(
        cell = cells_keep,
        lineage = colnames(pt)[j],
        pseudotime = pt[keep, j],
        weight = wt[keep, j],
        Cluster = factor(
          cluster[cells_keep],
          levels = c("C1", "C2", "C3", "C4", "C5")
        ),
        expr = as.numeric(
          gene_expr[cells_keep]
        )
      )
    }
  ) %>%
    bind_rows()
  
  p <- ggplot(
    plot_data,
    aes(
      x = pseudotime,
      y = expr
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
        plot_data,
        lineage == "Lineage2_1"
      ),
      aes(
        group = lineage,
        weight = weight
      ),
      size = 2,
      color = smoother_colors["Lineage2_1"],
      se = FALSE,
      method = "lm",
      formula = y ~ splines::ns(x, df = 2)
    ) +
    
    geom_smooth(
      data = subset(
        plot_data,
        lineage == "Lineage2_2"
      ),
      aes(
        group = lineage,
        weight = weight
      ),
      size = 2,
      color = smoother_colors["Lineage2_2"],
      se = FALSE,
      method = "lm",
      formula = y ~ splines::ns(x, df = 3)
    ) +
    
    scale_color_manual(
      name = "Cluster",
      values = cluster_colors
    ) +
    
    scale_alpha_continuous(
      range = c(0.02, 0.25),
      guide = "none"
    ) +
    
    labs(
      x = "Pseudotime",
      y = gene
    ) +
    
    theme_pubr(
      legend = "right"
    ) +
    
    theme(
      axis.title.y = element_text(
        face = "bold",
        size = 16
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 16
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(
        face = "bold",
        size = 16
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
    coord_cartesian(ylim = c(0, l[i])) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = 3,
          alpha = 0.5
        )
      )
    )
  
  plt_list[[gene]] <- p
}

print(
  wrap_plots(
    plt_list,
    guides = "collect",
    nrow = 2
  )
)


# =============================================================================
# Section 6 — Figures 5F–5H: revised C3 lineage grouping and SCENIC analysis
# =============================================================================

# C3 lineage-group assignment and the revised Figure 5H SCENIC analysis
# are generated in dedicated scripts:
#
#   Figure5H_define_C3_lineage_groups.R
#   Figure5H_SCENIC_consensus_analysis.R
#   Figure5H_SCENIC_figures.R
#
# The revised analysis uses Slingshot lineage weights >= 0.5 to define
# Lineage2_1, Shared and Lineage2_2, together with a 50-run pySCENIC
# consensus workflow.

###############################################################################
# SOURCE DATA EXPORT — Figure 5 and Figure S8
###############################################################################

source_dir <- "../SourceData/Figure5_S8"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# Figures 5A–5C
# Pseudobulk DESeq2 volcano plots
###############################################################################

fig5A_source$comparison <- "C3_vs_C4"
fig5B_source$comparison <- "C2_vs_C3"
fig5C_source$comparison <- "C5_vs_C3"

readr::write_tsv(
  dplyr::bind_rows(
    fig5A_source,
    fig5B_source,
    fig5C_source
  ),
  file.path(
    source_dir,
    "Fig5A_C_pseudobulk_DE.tsv"
  )
)


###############################################################################
# Figure 5D
# DEG Venn diagram
###############################################################################

fig5d_source <- dplyr::bind_rows(
  data.frame(
    comparison = "C4vC3",
    gene = as.character(C4vC3$genes)
  ),
  data.frame(
    comparison = "C3vC2",
    gene = as.character(C3vC2$genes)
  ),
  data.frame(
    comparison = "C3vC5",
    gene = as.character(C3vC5$genes)
  )
)

readr::write_tsv(
  fig5d_source,
  file.path(
    source_dir,
    "Fig5D_DEG_sets.tsv"
  )
)


###############################################################################
# Figure 5E and Figures S8A–S8C
#
# Cell-level expression values underlying the pseudotime scatter/line plots.
# Reconstruct independently so this does not depend on the current value of
# the 'genes', 'pt', 'wt', or 'exp' objects.
###############################################################################

trajectory_genes <- unique(c(
  
  # Figure 5E
  "PRG4",
  "FSTL4",
  "IGFBP5",
  "CDKN2A",
  
  # Figure S8A
  "INS",
  "TMED6",
  "TTR",
  "PPP1R1A",
  "CD82",
  "FXYD2",
  "ASCL2",
  
  # Figure S8B
  "S100A6",
  "IGFBP5",
  "FSTL4",
  "SGIP1",
  "INHBA",
  "SIX3",
  "NPY",
  "CDKN1A",
  "CDKN2A",
  
  # Figure S8C
  "SLC30A8",
  "GAD2",
  "NEUROD1",
  "MAFA",
  "ABCC8",
  "DNAJC15",
  "PCSK1",
  "LMO1"
))


# Reload trajectory object
sce_source <- readRDS("sce_lin2.rds")

crv_source <- as.data.frame(
  SummarizedExperiment::colData(sce_source)$crv
)

pt_source <- as.matrix(
  crv_source[, c(
    "pseudotime.Lineage1",
    "pseudotime.Lineage2"
  ), drop = FALSE]
)

wt_source <- as.matrix(
  crv_source[, c(
    "cellWeights.Lineage1",
    "cellWeights.Lineage2"
  ), drop = FALSE]
)

rownames(pt_source) <- colnames(sce_source)
rownames(wt_source) <- colnames(sce_source)

colnames(pt_source) <- c(
  "Lineage2_1",
  "Lineage2_2"
)

colnames(wt_source) <- c(
  "Lineage2_1",
  "Lineage2_2"
)


# Pull only genes actually plotted
trajectory_expr <- Seurat::GetAssayData(
  lin2_obj_full,
  assay = "RNA",
  layer = "data"
)[trajectory_genes, , drop = FALSE]


# Cluster identities
cluster_source <- setNames(
  paste0(
    "C",
    as.character(lin2_obj_full$monocle_clusters2)
  ),
  colnames(lin2_obj_full)
)


trajectory_source <- lapply(
  trajectory_genes,
  function(gene) {
    
    lapply(
      seq_len(ncol(pt_source)),
      function(j) {
        
        keep <- !is.na(pt_source[, j]) &
          wt_source[, j] > 0
        
        cells <- intersect(
          rownames(pt_source)[keep],
          colnames(trajectory_expr)
        )
        
        data.frame(
          cell_id = cells,
          gene = gene,
          lineage = colnames(pt_source)[j],
          pseudotime = pt_source[cells, j],
          weight = wt_source[cells, j],
          cluster = cluster_source[cells],
          expression = as.numeric(
            trajectory_expr[gene, cells]
          )
        )
      }
    ) %>%
      dplyr::bind_rows()
  }
) %>%
  dplyr::bind_rows()


# Annotate which panel(s) each gene belongs to
trajectory_source <- trajectory_source %>%
  dplyr::mutate(
    panel = dplyr::case_when(
      gene %in% c("IGFBP5", "FSTL4", "CDKN2A"
      ) ~ "Figure5E_and_S8B",
         
      gene %in% c(
        "PRG4", "FSTL4", "IGFBP5", "CDKN2A"
      ) ~ "Figure5E",
      
      gene %in% c(
        "INS", "TMED6", "TTR", "PPP1R1A",
        "CD82", "FXYD2", "ASCL2"
      ) ~ "FigureS8A",
      
      gene %in% c(
        "S100A6", "SGIP1", "INHBA", "SIX3",
        "NPY", "CDKN1A"
      ) ~ "FigureS8B",
      
      gene %in% c(
        "SLC30A8", "GAD2", "NEUROD1", "MAFA",
        "ABCC8", "DNAJC15", "PCSK1", "LMO1"
      ) ~ "FigureS8C"
     
  )


readr::write_tsv(
  trajectory_source,
  file.path(
    source_dir,
    "Fig5E_S8A_C_pseudotime_expression.tsv.gz"
  )
)

###############################################################################
# SOURCE DATA COLLECTION — revised Figure 5H SCENIC analysis
###############################################################################

source_dir <- "../SourceData/Figure5_S8"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

scenic_source_files <- c(
  "Figure5H_regulon_donor_data.tsv",
  "Figure5H_pairwise_mixed_model.tsv",
  "Figure5H_heatmap_source_data.tsv",
  "Figure5H_regulon_threshold_sensitivity.tsv",
  "Figure5H_regulon_robustness_summary.tsv",
  "Figure5H_sig70_regulons_across_thresholds.tsv",
  "Table_S16_SCENIC_regulon_robustness.tsv"
)

for (f in scenic_source_files) {
  
  if (file.exists(f)) {
    
    file.copy(
      f,
      file.path(source_dir, f),
      overwrite = TRUE
    )
    
  } else {
    
    warning("Missing SCENIC source-data file: ", f)
    
  }
}

###############################################################################
# sessionInfo()
###############################################################################

