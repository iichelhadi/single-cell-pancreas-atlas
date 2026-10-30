################################################################################
### R2-3: INS expression and sequencing-depth sensitivity
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(splines)
  library(SingleCellExperiment)
  library(slingshot)
  library(monocle3)
})
set.seed(123)

################################################################################
### 1. Recreate the Figure 2 Slingshot pseudotime
################################################################################

data.integrated <- readRDS('beta.integrated.rds')
cds <- readRDS('cds_all.rds')

## Figure 2 stores pseudotime in the Slingshot/Monocle objects, not
## necessarily in the Seurat metadata. Match cells explicitly.
stopifnot(setequal(colnames(data.integrated), colnames(cds)))
umap <- reducedDims(cds)$UMAP
stopifnot(setequal(colnames(data.integrated), rownames(umap)))
umap <- umap[colnames(data.integrated), , drop = FALSE]

sim <- SingleCellExperiment(assays = list(counts = GetAssayData(
  data.integrated, assay = 'SCT', slot = 'counts'
)))
reducedDims(sim)$UMAP <- umap
sim$monocle_clusters <- as.character(data.integrated$monocle_clusters)

## Same UMAP, clusters and root used in Figure2_clean.R.
sim <- slingshot(
  sim,
  clusterLabels = 'monocle_clusters',
  reducedDim = 'UMAP',
  start.clus = '3'
)
pt <- slingPseudotime(sim)
stopifnot(ncol(pt) == 3)

################################################################################
### 2. INS expression and sequencing depth
################################################################################

DefaultAssay(data.integrated) <- 'RNA'
ins_expr <- FetchData(data.integrated,
                      vars = c('INS', 'nCount_RNA', 'study'),
                      slot = 'data')
stopifnot(setequal(rownames(ins_expr), rownames(pt)))
pt <- pt[rownames(ins_expr), , drop = FALSE]
ins_expr$Lineage_1 <- pt[, 1]
ins_expr$Lineage_2 <- pt[, 2]
ins_expr$Lineage_3 <- pt[, 3]

raw_cor <- cor(ins_expr$INS,
               ins_expr$nCount_RNA,
               method = 'spearman',
               use = 'complete.obs')
write.table(
  data.frame(comparison = 'INS_vs_nCount_RNA', spearman_rho = raw_cor),
  'INS_raw_depth_correlation.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

################################################################################
### 3. Adjust INS for sequencing depth within each lineage
################################################################################

test_depth_adjusted <- function(df, lineage_col) {
  dat <- df |> filter(!is.na(.data[[lineage_col]]), !is.na(INS), !is.na(nCount_RNA))
  fit <- lm(INS ~ log1p(nCount_RNA), data = dat)
  dat$INS_depth_adjusted <- residuals(fit)
  data.frame(
    lineage = lineage_col,
    n_cells = nrow(dat),
    raw_rho = cor(dat$INS, dat[[lineage_col]], method = 'spearman'),
    depth_adjusted_rho = cor(dat$INS_depth_adjusted, dat[[lineage_col]], method = 'spearman')
  )
}

depth_results <- bind_rows(
  test_depth_adjusted(ins_expr, 'Lineage_1'),
  test_depth_adjusted(ins_expr, 'Lineage_2'),
  test_depth_adjusted(ins_expr, 'Lineage_3')
)
write.table(
  depth_results,
  'INS_depth_adjusted_correlations.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

################################################################################
### 4. Study- and depth-adjusted pseudotime models
################################################################################

model_results <- lapply(c('Lineage_2', 'Lineage_3'), function(lineage) {
  dat <- ins_expr[!is.na(ins_expr[[lineage]]), ]
  dat$pt <- dat[[lineage]]
  reduced <- lm(INS ~ log1p(nCount_RNA) + study, data = dat)
  full <- lm(INS ~ ns(pt, df = 3) + log1p(nCount_RNA) + study, data = dat)
  comparison <- anova(reduced, full)
  data.frame(
    lineage = lineage,
    F_statistic = comparison$F[2],
    p_value = comparison$`Pr(>F)`[2]
  )
})
model_results <- bind_rows(model_results)
write.table(
  model_results,
  'INS_depth_study_adjusted_models.tsv',
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
