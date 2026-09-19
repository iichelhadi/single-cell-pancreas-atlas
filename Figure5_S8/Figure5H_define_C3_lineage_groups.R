################################################################################
# Project: scRNA beta cell atlas — Figure 5H C3 lineage-group definition
# File: Figure5H_define_C3_lineage_groups.R
# Author: Elhadi Iich
#
# Description:
#   Defines C3 lineage groups from Slingshot lineage weights for revised
#   Figure 5H. Cells are assigned to Lineage2_1, Shared, or Lineage2_2
#   using a lineage-weight threshold of 0.5.
################################################################################



suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
})

set.seed(123)

################################################################################
# 1. Paths
################################################################################

LIN2_OBJECT <-
  "./revision/Figure4_S6/lin2.beta.reintegrated.rds"

SCE_OBJECT <-
  "./revision/Figure4_S6/sce_lin2.rds"

OUT_FILE <-
  "cellInfo_slingshotBranchID.Rds"

################################################################################
# 2. Load objects
################################################################################

lin2_obj <- readRDS(LIN2_OBJECT)
sce <- readRDS(SCE_OBJECT)

cat("Full Lineage 2 Seurat object:", ncol(lin2_obj), "cells\n")
cat("Full Slingshot SCE:", ncol(sce), "cells\n")

################################################################################
# 3. Identify C3 cells
################################################################################

colData(sce)$Clusters <- lin2_obj$monocle_clusters2

c3_cells <- colnames(sce)[
  as.character(colData(sce)$Clusters) == "3"
]

cat("C3 cells:", length(c3_cells), "\n")

lin2_c3 <- subset(
  lin2_obj,
  cells = c3_cells
)

sce_c3 <- sce[, c3_cells]

################################################################################
# 4. Extract Slingshot pseudotime and lineage weights
################################################################################

cellInfo <- as.data.frame(
  colData(sce_c3)$crv
)

rownames(cellInfo) <- colnames(sce_c3)

cat("\nColumns extracted from Slingshot:\n")
print(colnames(cellInfo))

required <- c(
  "pseudotime.Lineage1",
  "pseudotime.Lineage2",
  "cellWeights.Lineage1",
  "cellWeights.Lineage2"
)

w1 <- cellInfo$cellWeights.Lineage1
w2 <- cellInfo$cellWeights.Lineage2

################################################################################
# 5. Define branch membership
################################################################################

BRANCH_THRESHOLD <- 0.5

cellInfo$LineageGroup <- case_when(
  
  w1 >= BRANCH_THRESHOLD &
    w2 < BRANCH_THRESHOLD ~ "Lineage2_1",
  
  w1 >= BRANCH_THRESHOLD &
    w2 >= BRANCH_THRESHOLD ~ "Shared",
  
  w1 < BRANCH_THRESHOLD &
    w2 >= BRANCH_THRESHOLD ~ "Lineage2_2",
  
  TRUE ~ NA_character_
)

cellInfo$LineageGroup <- factor(
  cellInfo$LineageGroup,
  levels = c(
    "Lineage2_1",
    "Shared",
    "Lineage2_2"
  )
)

################################################################################
# 6. Sanity checks
################################################################################

cat("\nC3 lineage groups:\n")
print(
  table(
    cellInfo$LineageGroup,
    useNA = "ifany"
  )
)
# Lineage2_1     Shared Lineage2_2 
# 765        637       1360 


################################################################################
# 7. Add useful continuous lineage information
################################################################################

cellInfo$weight_difference <-
  w1 - w2

cellInfo$absolute_weight_difference <-
  abs(w1 - w2)

cellInfo$DominantLineage <- case_when(
  
  w1 > w2 ~ "Lineage2_1",
  
  w2 > w1 ~ "Lineage2_2",
  
  TRUE ~ "Equal"
)

################################################################################
# 8. Add relevant biological metadata
################################################################################

meta_c3 <- lin2_c3@meta.data[
  rownames(cellInfo),
  ,
  drop = FALSE
]


for (x in c(
  "Condition",
  "Hba1c_class_median",
  "Study",
  "orig.ident",
  "sample_name"
)) {
  
  if (x %in% colnames(meta_c3)) {
    cellInfo[[x]] <- meta_c3[[x]]
  }
}

################################################################################
# 9. Diagnostic weight plot
################################################################################

p_weights <- ggplot(
  cellInfo,
  aes(
    x = cellWeights.Lineage1,
    y = cellWeights.Lineage2,
    color = LineageGroup
  )
) +
  geom_point(
    size = 0.8,
    alpha = 0.4
  ) +
  geom_vline(
    xintercept = BRANCH_THRESHOLD,
    linetype = 2
  ) +
  geom_hline(
    yintercept = BRANCH_THRESHOLD,
    linetype = 2
  ) +
  scale_color_manual(
    values = c(
      "Lineage2_1" = "#00B0F6",
      "Shared" = "#9590FF",
      "Lineage2_2" = "#E76BF3"
    )
  ) +
  labs(
    x = "Slingshot weight: Lineage2_1",
    y = "Slingshot weight: Lineage2_2",
    color = NULL
  ) +
  theme_classic()

print(p_weights)

################################################################################
# 10. C3 UMAP using new branch definition
################################################################################

lin2_c3 <- AddMetaData(
  lin2_c3,
  metadata = cellInfo["LineageGroup"]
)

p_umap <- DimPlot(
  lin2_c3,
  reduction = "umap_MNN",
  group.by = "LineageGroup"
) +
  scale_color_manual(
    values = c(
      "Lineage2_1" = "#00B0F6",
      "Shared" = "#9590FF",
      "Lineage2_2" = "#E76BF3"
    )
  ) +
  ggpubr::theme_pubr(legend = 'bottom') +
  ggpubr::font("title", size = 16, color = "black", face = "bold") +
  NoAxes() + NoGrid() + ggtitle('Lineage group') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.text = element_text(face = 'bold', size = 14))

print(p_umap)

################################################################################
# 11. Source data for Figure 5F-G
################################################################################

SOURCE_DIR <-
  "./revision/SourceData/Figure5_S8"

dir.create(
  SOURCE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Extract UMAP coordinates used in the plot
umap_c3 <- Embeddings(
  lin2_c3,
  reduction = "umap_MNN"
)

umap_c3 <- umap_c3[
  rownames(cellInfo),
  ,
  drop = FALSE
]


# Build source-data table
fig5fg_source <- data.frame(
  cell_id = rownames(cellInfo),
  
  Lineage2_1_weight =
    cellInfo$cellWeights.Lineage1,
  
  Lineage2_2_weight =
    cellInfo$cellWeights.Lineage2,
  
  Lineage2_1_pseudotime =
    cellInfo$pseudotime.Lineage1,
  
  Lineage2_2_pseudotime =
    cellInfo$pseudotime.Lineage2,
  
  LineageGroup =
    as.character(cellInfo$LineageGroup),
  
  UMAP_1 =
    umap_c3[, 1],
  
  UMAP_2 =
    umap_c3[, 2],
  
  stringsAsFactors = FALSE
)


write.table(
  fig5fg_source,
  file = gzfile(
    file.path(
      SOURCE_DIR,
      "Fig5F_G_C3_lineage_groups.tsv.gz"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nSaved Figure 5F-G source data:",
  nrow(fig5fg_source),
  "cells\n"
)

print(
  table(fig5fg_source$LineageGroup)
)

################################################################################
# 12. Save clean metadata for Figure 5H
################################################################################

saveRDS(
  cellInfo,
  OUT_FILE
)

################################################################################
# sessionInfo()
################################################################################

