###############################################################################
# Script: lin2C5_SCENIC_regulon_RIS_plot.R
# Purpose:
#   - Load SCENIC regulon AUCs and regulonTargetsInfo
#   - Build TF→target graph weighted by NES
#   - Combine regulon AUC (AUC) and target NES into a composite
#     Regulon Importance Score (RIS)
#   - Visualize the GRN highlighting high-confidence TFs and high-RIS hubs
#
# Notes:
#   * Original logic preserved verbatim. Only comments/headers added.
#   * Paths are unchanged and point to your HPC directories.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(Seurat)
  library(Cairo)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(123)

###############################################################################
# 1) Load regulon AUC matrix and tidy into (regulons x cells)
###############################################################################

AUCmatrix <- readRDS("L2_C5_regulonAUC.Rds")
# AUCmatrix <- readRDS("L3_C3_regulonAUC.Rds")

# Extract numeric AUC matrix and transpose to cells x regulons, then to regulons x cells
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names = FALSE)
RegulonName_AUC <- colnames(AUCmatrix)

# Keep only non-duplicated "extended" regulons (SCENIC helper)
RegulonName_AUC <- unname(onlyNonDuplicatedExtended(RegulonName_AUC))
AUCmatrix <- AUCmatrix[, RegulonName_AUC]
dim(AUCmatrix)

# Clean up regulon names to plain "TF" or "TF_extended" forms
RegulonName_AUC <- gsub(' \\(', '_', RegulonName_AUC)
RegulonName_AUC <- gsub('\\)', '', RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
dim(AUCmatrix)

# Flip back to regulons x cells
AUCmatrix <- data.frame(t(AUCmatrix))
dim(AUCmatrix)

# Harmonize cell names (replace '.' with '-')
colnames(AUCmatrix) <- gsub(pattern = '[.]', replacement = '-', x = colnames(AUCmatrix))

# Mean AUC per regulon (aggregates across cells)
meanAUC <- rowMeans(AUCmatrix)

# Collapse “TF_extended” → “TF” (drop the suffix after first underscore)
names(meanAUC) <- str_split_fixed(names(meanAUC), pattern = '_', n = 2)[, 1]

###############################################################################
# 2) Load regulon targets with NES and build TF→Target graph
###############################################################################

regulonTargetsInfo <- readRDS('L2_C5_regulonTargetsInfo.Rds')
# regulonTargetsInfo <- readRDS('L3_C3_regulonTargetsInfo.Rds')

# Edge list: TF, Target, high/low confidence, and NES (edge weight)
edges <- data.frame(
  TF = regulonTargetsInfo$TF,
  TargetGene = regulonTargetsInfo$gene,
  InteractionType = ifelse(regulonTargetsInfo$highConfAnnot == "TRUE", "High Confidence", "Low Confidence"),
  NES = as.numeric(regulonTargetsInfo$NES),
  stringsAsFactors = FALSE
)

# Keep only TFs that have meanAUC (ensures the two metrics intersect later)
edges <- edges[edges$TF %in% names(meanAUC), ]

library(igraph)

# Directed TF→Target graph; NES used as edge weight
graph <- graph_from_data_frame(edges, directed = TRUE)
E(graph)$weight <- edges$NES

###############################################################################
# 3) Combine mean AUC (per TF) with mean NES (per TF) → RIS
###############################################################################

# Average NES per TF across its targets
meanNES <- tapply(edges$NES, edges$TF, mean, na.rm = TRUE)

# Intersect TF sets present in both meanAUC and meanNES
common_regulons <- intersect(names(meanAUC), names(meanNES))
print(common_regulons)

# Align and z-score each metric before combining
meanAUC_common <- meanAUC[common_regulons]
meanNES_common <- meanNES[common_regulons]
normAUC <- as.vector(scale(meanAUC_common))
normNES <- as.vector(scale(meanNES_common))

# Regulon Importance Score (RIS) = z(AUC) + z(NES)
regulonMetrics <- data.frame(
  Regulon = common_regulons,
  normAUC = normAUC,
  normNES = normNES,
  stringsAsFactors = FALSE
)
regulonMetrics$RIS <- regulonMetrics$normAUC + regulonMetrics$normNES
print(regulonMetrics)

###############################################################################
# 4) Annotate graph vertices: confidence, RIS, labels, sizes
###############################################################################

# Mark TFs with high-confidence interactions present at least once
high_confidence_TFs <- edges %>%
  dplyr::filter(InteractionType == "High Confidence") %>%
  dplyr::pull(TF) %>%
  unique()

V(graph)$HighConfidence <- ifelse(V(graph)$name %in% high_confidence_TFs, "High", "Low")
V(graph)$alpha <- ifelse(V(graph)$name %in% high_confidence_TFs, 0.9, 0.2)

# Attach RIS to TF vertices (non-TFs get a tiny value to keep them small)
V(graph)$RIS <- ifelse(
  V(graph)$name %in% regulonMetrics$Regulon,
  regulonMetrics$RIS[match(V(graph)$name, regulonMetrics$Regulon)],
  1e-7
)

# Label only high-confidence TFs in the top 20% by RIS
threshold <- quantile(V(graph)$RIS, 0.8)
(V(graph)$name %in% high_confidence_TFs & V(graph)$RIS > threshold) %>% sum()

V(graph)$Label <- ifelse(
  V(graph)$name %in% high_confidence_TFs & V(graph)$RIS > threshold,
  V(graph)$name, NA
)

unique(V(graph)$Label)

# Mark TFs for point sizing
V(graph)$isTF <- ifelse(V(graph)$Label %in% high_confidence_TFs, TRUE, FALSE)
V(graph)$adjustedSize <- ifelse(V(graph)$isTF, V(graph)$RIS, 1e-7)

###############################################################################
# 5) Plot GRN with ggraph (KK layout), scaled by RIS, colored by confidence
###############################################################################

library(ggraph)

ggraph(graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), color = "grey60") +
  geom_node_point(
    aes(size = adjustedSize, fill = HighConfidence, alpha = alpha),
    shape = 21, color = "black"
  ) +
  geom_node_text(
    aes(label = Label, color = 'black'),
    size = 4, fontface = "bold", repel = TRUE,
    box.padding = 0.5, segment.size = 0.2
  ) +
  scale_fill_manual(values = c("High" = "#0baedd",
                               # "High" = "#53B400",
                               "Low" = "black"), guide = "none") +
  scale_color_identity() +
  scale_size_continuous(name = "Regulon Importance Score (RIS)", range = c(2, 15)) +
  scale_edge_alpha_continuous(name = "NES", range = c(0.2, 1)) +
  guides(alpha = "none") +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.box.margin = margin(10, 10, 10, 10)
  )

# Optional export:
# write_graph(graph, file = "graph_lin2C5.graphml", format = "graphml")
# write_graph(graph, file = "graph_lin3C3.graphml", format = "graphml")