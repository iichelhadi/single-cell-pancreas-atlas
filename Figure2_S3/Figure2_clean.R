###############################################################################
# Project: scRNA beta cell atlas Fig 2 and lineage analysis
# File: Fig2.R
# Author: Elhadi Iich
# Date: 2026-09-07
#
# Description:
#   Generates Fig 2 panels and lineage analyses using Monocle3 and Slingshot.
#
# Inputs:
#   - beta.integrated.rds
#   - cds_all.rds
#
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(reticulate)
  library(monocle3)
  library(slingshot)
})

set.seed(123)

cores = 1
resolution = 0.00004

palette <- c("#5A5156FF", "#E4E1E3FF", "#F6222EFF" ,
             "#FE00FAFF", "#16FF32FF" ,"#3283FEFF",
             "#FEAF16FF", "#B00068FF", "#1CFFCEFF",
             "#90AD1CFF", "#2ED9FFFF" ,"#DEA0DFFF",
             "#AA0EDFFF", "#F8A19FFF", "#325A9BFF" ,
             "#C4451CFF", "#1C8356FF" ,"#85660DFF",
             "#B10DA1FF", "#FBE426FF", "#1CBE4FFF")



data.integrated <-
  readRDS('beta.integrated.rds')

cds <- readRDS('cds_all.rds')

# cds@clusters@listData$UMAP$clusters <- data.integrated$monocle_clusters

# colData(cds)$monocle_clusters <- data.integrated$monocle_clusters

df <- data.integrated@meta.data %>% as.data.frame()
df$monocle_clusters <- as.character(df$monocle_clusters) 

rownames(df) <- rownames(data.integrated@meta.data)

desired_order <- df %>%
  group_by( monocle_clusters,Condition) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  mutate(per = round(100 *n/sum(n),1)) %>% 
  filter(Condition == 'T2D') %>%
  arrange(per) %>% pull(monocle_clusters)

# Initialize an empty list
numbers <- list()


df %>%
  group_by( monocle_clusters,Gender) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  mutate(per = paste0(round(100 *n/sum(n),1),'%')) %>%
  ggplot(aes(x = monocle_clusters, y = n, fill = Gender)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) +
  scale_x_discrete(limits = desired_order) +
  geom_text(aes( label = per),position = position_fill(vjust = 0.5)) +
  ylab('Proportion') +xlab(NULL)

cds_subset_INS <- cds[rownames(cds)=='INS',]
cds_subset_GCG <- cds[rownames(cds)=='GCG',]
cds_subset_PPY <- cds[rownames(cds)=='PPY',]
cds_subset_SST <- cds[rownames(cds)=='SST',]

########################################################################################

##Fig 2A ####
A2 <- plot_cells(cds, show_trajectory_graph = F, 
                 color_cells_by = 'cluster',
                 label_cell_groups=T, cell_size = 1,
                 group_cells_by = 'cluster', 
                 group_label_size = 12)+
  ggpubr::clean_theme()
A2

## Fig 2B-1 - Age distribution ####

B2_1 <- df %>%
  drop_na(Age) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  ggplot(
    aes(
      x = monocle_clusters,
      y = Age,
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
    x = NULL,
    y = "Age (years)"
  ) +
  theme(
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16)
  )

## Fig 2B-2 - BMI distribution ####

B2_2 <- df %>%
  drop_na(BMI) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
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

##Fig 2A-2B ####
A2/(B2_1/B2_2)

##Fig 2C ####

C2 <- plot_cells(cds, show_trajectory_graph = F, 
                 color_cells_by = 'study',
                 label_cell_groups=FALSE, cell_size = 1,
                 group_cells_by = 'cluster', 
                 group_label_size = 4)+
  ggpubr::clean_theme()+NoLegend()+
  scale_color_manual(values = palette[-c(12,13)])

##Fig 2D ####

D2 <- df %>%
  group_by( monocle_clusters,study) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>% 
  group_by(monocle_clusters) %>% 
  mutate(per = paste0(round(100 *n/sum(n),1),'%')) %>%
  ggplot(aes(x = monocle_clusters, y = n, fill = study)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
  scale_x_discrete(limits = desired_order) + 
  # geom_text(aes( label = per),position = position_fill(vjust = 0.5)) +
  ylab('Proportion') +xlab("Cluster")+
  scale_fill_manual(values = palette[-c(12,13)])+
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 16),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 7),
    legend.title  = element_text(face = 'bold', size = 14)
  )

##Fig 2E ####

E2 <- plot_cells(cds, show_trajectory_graph = F, 
                 color_cells_by = 'Condition',
                 label_cell_groups=FALSE, cell_size = 1,
                 group_cells_by = 'cluster', 
                 group_label_size = 4)+
  ggpubr::clean_theme()+NoLegend()

##Fig 2F ####
F2 <- df %>%
  group_by( monocle_clusters,Condition) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>% 
  group_by(monocle_clusters) %>% 
  mutate(per = paste0(round(100 *n/sum(n),1),'%')) %>%
  ggplot(aes(x = monocle_clusters, y = n, fill = Condition)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
  scale_x_discrete(limits = desired_order)+
  ylab('Proportion') +xlab('Cluster') + 
  geom_text(aes( label = per),size = 4, fontface = 'bold',
            position = position_fill(vjust = 0.5))+
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 16),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 10),
    legend.title  = element_text(face = 'bold', size = 14)
  )

## Figs 2C-F ####
(C2/D2)|(E2/F2)

#### Fig 2G ####

G2_1 <- plot_genes_violin(
  cds_subset_INS,
  group_cells_by = "monocle_clusters"
) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_discrete(limits = desired_order) +
  ggpubr::theme_pubr() +
  theme(
    strip.background = element_blank(), 
    axis.line.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14)
  ) +
  xlab(NULL) +
  ylab("Log10(expression)")


G2_2 <- plot_genes_violin(cds_subset_GCG, group_cells_by = "monocle_clusters") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_discrete(limits = desired_order) +
  ggpubr::theme_pubr() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14)
  ) +
  xlab("Cluster") +
  ylab("Log10(expression)")


G2_3 <- plot_genes_violin(cds_subset_PPY, group_cells_by = "monocle_clusters") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_discrete(limits = desired_order) +
  ggpubr::theme_pubr() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14)
  ) +
  xlab("Cluster") +
  ylab("Log10(expression)")


G2_4 <- plot_genes_violin(cds_subset_SST, group_cells_by = "monocle_clusters") +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  scale_x_discrete(limits = desired_order) +
  ggpubr::theme_pubr() +
  theme(
    strip.background = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 14)
  ) +
  xlab(NULL) +
  ylab("Log10(expression)")

##Fig 2G
(G2_1/G2_2|G2_4/G2_3)


##Fig 2H

umap_cds <- reducedDims(cds)$UMAP

data.integrated[['umap_MNN']] <-
  CreateDimReducObject(embeddings = umap_cds,
                       key = "umap_MNN_")

DefaultAssay(data.integrated) <- 'SCT'

H2_1 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('INS', 'GCG'),
  pt.size = 0.5,
  order = T,
  blend = T, 
  cols = c('black','blue','yellow')
)

H2_2 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('INS', 'SST'),
  pt.size = 0.5,
  order = T,
  blend = T, 
  cols = c('black','blue','yellow')
)

H2_3 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('INS', 'PPY'),
  pt.size = 0.5,
  order = T,
  blend = T, 
  cols = c('black','blue','yellow')
)

(H2_1/H2_2/H2_3)

################################################################################
################################################################################

################################################################################
################################# Slingshot ####################################
################################################################################

sim <-
  SingleCellExperiment(assays =
                         List(counts =
                                as.matrix(data.integrated@assays$SCT@counts)))

umap_cds <- reducedDims(cds)$UMAP
reducedDims(sim) = SimpleList(UMAP = umap_cds)

# metadata
colData(sim)$study = df$study
colData(sim)$monocle_clusters =  df$monocle_clusters
colData(sim)$Condition = df$Condition

sim <- slingshot(sim, 
                 clusterLabels = 'monocle_clusters',  
                 reducedDim = 'UMAP',  
                 start.clus= "3",  
                 end.clus = NULL   
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

##Fig 2J
umap_plot <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = 'cluster',
  label_groups_by_cluster = F,
  cell_size = 1,label_cell_groups = T,
  group_cells_by = 'cluster',
  group_label_size = 14
)

# Create a trajectory plot
trajectory_plot <- umap_plot + geom_path(
  data = combined_lineages,
  aes(
    x = `Dim-1`,
    y = `Dim-2`,
    color = lineage,
    group = lineage
    
  ),
  size = 2
)+ scale_color_manual(values = c('1'="#F8766D",'2'= "#C49A00",'3'= "#53B400",
                                 '4'= "#00C094",'5'= "#00B6EB" ,'6'="#A58AFF" ,
                                 '7'="#FB61D7",'Lineage2'='red',
                                 'Lineage1'='blue','Lineage3'="#D89000", 
                                 'Lineage4'="#39B600", 'Lineage5'="#9590FF"))

J2 <- trajectory_plot+
  ggpubr::clean_theme()



colData(cds)$Lineage1 <- sim@colData$slingPseudotime_1
colData(cds)$Lineage2 <- sim@colData$slingPseudotime_2
colData(cds)$Lineage3 <- sim@colData$slingPseudotime_3

# N2_1 <- plot_cells(
#   cds,
#   show_trajectory_graph = FALSE,
#   color_cells_by = 'Lineage1',
#   cell_size = 1,
#   group_label_size = 10
# )+ggpubr::clean_theme()+
#   theme(
#     legend.text   = element_text(face = 'bold', size = 14),
#     legend.title  = element_text(face = 'bold', size = 14)
#   )

# N2_2 <- plot_cells(
#   cds,
#   show_trajectory_graph = FALSE,
#   color_cells_by = 'Lineage2',
#   cell_size = 1,
#   group_label_size = 10
# )+ggpubr::clean_theme()+
#   theme(
#     legend.text   = element_text(face = 'bold', size = 14),
#     legend.title  = element_text(face = 'bold', size = 14)
#   )


# N2_3 <- plot_cells(
#   cds,
#   show_trajectory_graph = FALSE,
#   color_cells_by = 'Lineage3',
#   cell_size = 1,
#   group_label_size = 10
# )+ggpubr::clean_theme()+
#   theme(
#     legend.text   = element_text(face = 'bold', size = 14),
#     legend.title  = element_text(face = 'bold', size = 14)
#   )


# M2 + N2_1 + N2_2 + N2_3

#######################################################################################
####################### Lineage subsets generation #################################
#######################################################################################

lin1_obj <- subset(data.integrated, cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage1)])
dim(lin1_obj)
saveRDS(lin1_obj, 'Beta_T2D_lineage1.rds')
lin1_obj <- SetIdent(lin1_obj, value = 'assigned_label')
DefaultAssay(lin1_obj) <- 'RNA'
lin1_obj %>% 
  scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
  write.table(., 'Lineage1_sling_perc_expr.tsv',row.names = T,
              col.names = NA, sep = '\t')


lin2_obj <- subset(data.integrated, 
                   cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage2)])
dim(lin2_obj)
saveRDS(lin2_obj, 'Beta_T2D_lineage2.rds')
lin2_obj <- SetIdent(lin2_obj, value = 'assigned_label')
DefaultAssay(lin2_obj) <- 'RNA'
lin2_obj %>% 
  scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
  write.table(., 'Lineage2_sling_perc_expr.tsv',
              row.names = T,col.names = NA, sep = '\t')


lin3_obj <- subset(data.integrated, 
                   cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage3)])
dim(lin3_obj)
saveRDS(lin3_obj, 'Beta_T2D_lineage3.rds')
lin3_obj <- SetIdent(lin3_obj, value = 'assigned_label')
DefaultAssay(lin3_obj) <- 'RNA'
lin3_obj %>% 
  scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
  write.table(., 'Lineage3_sling_perc_expr.tsv',
              row.names = T,col.names = NA, sep = '\t')

################################################################################
################################################################################

################################################################################
## Fig 2K — INS expression along lineage-specific pseudotime
################################################################################

DefaultAssay(data.integrated) <- "RNA"

ins_expr <- GetAssayData(
  data.integrated,
  assay = "RNA",
  slot = "data"
)["INS", ]

pt <- slingPseudotime(sim)
wt <- slingCurveWeights(sim)

colnames(pt) <- paste0("Lineage", seq_len(ncol(pt)))
colnames(wt) <- paste0("Lineage", seq_len(ncol(wt)))

plot_data <- lapply(
  seq_len(ncol(pt)),
  function(i) {
    
    keep <- !is.na(pt[, i]) & wt[, i] > 0
    
    x <- pt[keep, i]
    
    tibble(
      cell = rownames(pt)[keep],
      lineage = paste0("Lineage", i),
      pseudotime = (x - min(x)) / (max(x) - min(x)),
      INS = as.numeric(ins_expr[rownames(pt)[keep]])
    )
  }
) %>%
  bind_rows()

point_colors <- c(
  "Lineage1" = "#89CFF0",
  "Lineage2" = "#FF9999",
  "Lineage3" = "#FFB347"
)

smoother_colors <- c(
  "Lineage1" = "blue",
  "Lineage2" = "red",
  "Lineage3" = "#d89000"
)

K2 <- ggplot(
  plot_data,
  aes(
    x = pseudotime,
    y = INS
  )
) +
  geom_point(
    aes(color = lineage),
    size = 0.6,
    alpha = 0.6
  ) +
  geom_smooth(
    data = subset(plot_data, lineage == "Lineage1"),
    aes(group = lineage),
    method = "lm",
    formula = y ~ splines::ns(x, df = 2),
    se = FALSE,
    linewidth = 2,
    color = smoother_colors["Lineage1"]
  ) +
  geom_smooth(
    data = subset(plot_data, lineage == "Lineage2"),
    aes(group = lineage),
    method = "lm",
    formula = y ~ splines::ns(x, df = 2),
    se = FALSE,
    linewidth = 2,
    color = smoother_colors["Lineage2"]
  ) +
  geom_smooth(
    data = subset(plot_data, lineage == "Lineage3"),
    aes(group = lineage),
    method = "lm",
    formula = y ~ splines::ns(x, df = 2),
    se = FALSE,
    linewidth = 2,
    color = smoother_colors["Lineage3"]
  ) +
  scale_color_manual(
    name = "Lineages",
    values = point_colors
  ) +
  labs(
    x = "Pseudotime",
    y = "INS"
  ) +
  ggpubr::theme_pubr(
    legend = "right"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3)
    )
  )

K2

###############################################################################
# SOURCE DATA EXPORT — Fig 2
###############################################################################

source_dir <- "../SourceData/Fig2_S3"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Figs 2A, 2B, 2C, 2E, 2G and 2H
#
# Shared cell-level data:
# UMAP coordinates, metadata, age/BMI and the four genes plotted.
###############################################################################

umap_source <- as.data.frame(
  SingleCellExperiment::reducedDims(cds)$UMAP
)

colnames(umap_source)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

umap_source <- umap_source %>%
  tibble::rownames_to_column("cell_id")


meta_source <- data.integrated@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    monocle_clusters,
    study,
    Condition,
    Age,
    BMI
  )


fig2_genes <- c(
  "INS",
  "GCG",
  "SST",
  "PPY"
)

fig2_expr <- Seurat::GetAssayData(
  data.integrated,
  assay = "RNA",
  layer = "data"
)[fig2_genes, , drop = FALSE]

fig2_expr <- as.data.frame(
  t(as.matrix(fig2_expr))
) %>%
  tibble::rownames_to_column("cell_id")


fig2_cell_source <- umap_source %>%
  dplyr::left_join(
    meta_source,
    by = "cell_id"
  ) %>%
  dplyr::left_join(
    fig2_expr,
    by = "cell_id"
  )

readr::write_tsv(
  fig2_cell_source,
  file.path(
    source_dir,
    "Fig2A_B_C_E_G_H_cell_data.tsv.gz"
  )
)


###############################################################################
# Figs 2D and 2F
# Study and condition composition
###############################################################################

fig2d_source <- meta_source %>%
  dplyr::count(
    monocle_clusters,
    study,
    name = "cell_count"
  ) %>%
  dplyr::group_by(
    monocle_clusters
  ) %>%
  dplyr::mutate(
    proportion = cell_count /
      sum(cell_count),
    variable = "Study",
    category = as.character(study)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = monocle_clusters,
    variable,
    category,
    cell_count,
    proportion
  )


fig2f_source <- meta_source %>%
  dplyr::count(
    monocle_clusters,
    Condition,
    name = "cell_count"
  ) %>%
  dplyr::group_by(
    monocle_clusters
  ) %>%
  dplyr::mutate(
    proportion = cell_count /
      sum(cell_count),
    variable = "Condition",
    category = as.character(Condition)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = monocle_clusters,
    variable,
    category,
    cell_count,
    proportion
  )

readr::write_tsv(
  dplyr::bind_rows(
    fig2d_source,
    fig2f_source
  ),
  file.path(
    source_dir,
    "Fig2D_F_composition.tsv"
  )
)


###############################################################################
# Fig 2J
# Slingshot trajectory curves
###############################################################################

fig2j_slingshot_source <- combined_lineages %>%
  dplyr::rename(
    UMAP_1 = `Dim-1`,
    UMAP_2 = `Dim-2`
  ) %>%
  dplyr::select(
    lineage,
    path_id,
    UMAP_1,
    UMAP_2
  )

readr::write_tsv(
  fig2j_slingshot_source,
  file.path(
    source_dir,
    "Fig2J_slingshot_curves.tsv"
  )
)


###############################################################################
# Fig 2K
# INS expression along normalized lineage-specific pseudotime
###############################################################################

fig2k_source <- plot_data %>%
  dplyr::rename(
    cell_id = cell
  ) %>%
  dplyr::select(
    cell_id,
    lineage,
    pseudotime,
    INS
  )

readr::write_tsv(
  fig2k_source,
  file.path(
    source_dir,
    "Fig2K_INS_pseudotime.tsv.gz"
  )
)


###############################################################################
# COLLECT SOURCE DATA — Fig 2I and Fig S3F
###############################################################################

outdir <- "Figure2_S3"
source_dir <- "../SourceData/Fig2_S3"

dir.create(
  source_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

source_files <- c(
  "Fig2I_HbA1c_donor_distribution_source_data.tsv",
  "Fig2I_HbA1c_cluster_summary.tsv",
  "FigS3F_donor_cluster_occupancy_source_data.tsv",
  "FigS3F_HbA1c_cluster_occupancy_associations.tsv",
  "FigS3F_HbA1c_threshold_occupancy_sensitivity.tsv"
)

for (f in source_files) {
  
  src <- file.path(outdir, f)
  
  if (file.exists(src)) {
    
    file.copy(
      src,
      file.path(source_dir, f),
      overwrite = TRUE
    )
    
  } else {
    
    warning("Missing source-data file: ", src)
  }
}


###############################################################################
# sessionInfo()
###############################################################################
