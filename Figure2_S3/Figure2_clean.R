###############################################################################
# Project: scRNA beta cell atlas Figure 2 and lineage analysis
# File: Figure2.R
# Author: Elhadi Iich
# Date: 2025-08-27
#
# Description:
#   Generates Figure 2 panels and lineage analyses using Monocle3 and Slingshot.
#   All original code is preserved verbatim. Only comments and section banners
#   were added for readability on GitHub.
#
# Inputs:
#   - beta.integrated.rds
#   - cds_all.rds
#   - sce_all.rds
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

cds@clusters@listData$UMAP$clusters <- data.integrated$monocle_clusters

colData(cds)$monocle_clusters <- data.integrated$monocle_clusters

df <- data.integrated@meta.data %>% as.data.frame()
df$monocle_clusters <- as.character(df$monocle_clusters) 
df2 <- df %>% as.data.frame()%>% group_by(monocle_clusters) %>%
  summarise(
    mean_value = round(mean(Hba1c, na.rm = TRUE),digits = 2),
    median_value = median(Hba1c, na.rm = T)) %>%
  mutate(
    Hba1c_class = case_when(
      mean_value < 5.7 ~ "ND",
      mean_value >= 5.7 & mean_value < 6.5 ~ "pre-T2D",
      mean_value >= 6.5 ~ "T2D"
    ),
    Hba1c_class2 = case_when(
      median_value < 5.7 ~ "ND",
      median_value >= 5.7 & median_value < 6.5 ~ "pre-T2D",
      median_value >= 6.5 ~ "T2D")
  )

df <- left_join(df, df2, by = 'monocle_clusters')

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

# Loop over each unique monocle_cluster
for (i in unique(df$monocle_clusters)) {
  # Count the unique orig.ident in each cluster and append it to the list
  numbers <- append(numbers, length(unique(df$orig.ident[df$monocle_clusters == i])))
}

# Create a dataframe with the results
numbers_df <- data.frame(Clusters = unique(df$monocle_clusters), 
                         n_Samples = unlist(numbers))


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
cds_subset_GHRL <- cds[rownames(cds)=='GHRL',]

########################################################################################

##Figure 2A ####
A2 <- plot_cells(cds, show_trajectory_graph = F, 
                 color_cells_by = 'cluster',
                 label_cell_groups=T, cell_size = 1,
                 group_cells_by = 'cluster', 
                 group_label_size = 12)+
  ggpubr::clean_theme()
A2
##Figure 2B ####
B2 <- df %>% drop_na(Age) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  summarise(
    median_value = median(Age, na.rm = TRUE),
    mad_value = mad(Age, na.rm = TRUE)
  )  %>%
  ggplot(aes(x = monocle_clusters, y = median_value, fill = monocle_clusters)) +
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
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 18),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 7),
    legend.title  = element_text(face = 'bold', size = 14)
  )+
  xlab(NULL)

##Figure 2C ####
C2 <- df %>% drop_na(BMI) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  summarise(
    median_value = median(BMI, na.rm = TRUE),
    mad_value = mad(BMI, na.rm = TRUE)
  )  %>%
  ggplot(aes(x = monocle_clusters, y = median_value, fill = monocle_clusters)) +
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
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 18),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 7),
    legend.title  = element_text(face = 'bold', size = 14)
  )

##Figure 2A-2C ####
A2/(B2/C2)

##Figure 2D ####

D2 <- plot_cells(cds, show_trajectory_graph = F, 
                 color_cells_by = 'Study',
                 label_cell_groups=FALSE, cell_size = 1,
                 group_cells_by = 'cluster', 
                 group_label_size = 4)+
  ggpubr::clean_theme()+NoLegend()+
  scale_color_manual(values = palette[-c(12,13)])

E2 <- df %>%
  group_by( monocle_clusters,Study) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>% 
  group_by(monocle_clusters) %>% 
  mutate(per = paste0(round(100 *n/sum(n),1),'%')) %>%
  ggplot(aes(x = monocle_clusters, y = n, fill = Study)) +
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

##Figure 2F ####
F2 <- plot_cells(cds, show_trajectory_graph = F, 
                 color_cells_by = 'Condition',
                 label_cell_groups=FALSE, cell_size = 1,
                 group_cells_by = 'cluster', 
                 group_label_size = 4)+
  ggpubr::clean_theme()+NoLegend()

##Figure 2G ####
G2 <- df %>%
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

## Figures 2D-G ####
(D2/E2)|(F2/G2)

#### Figure 2H ####
H2_1<- plot_genes_violin(cds_subset_INS, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(), 
        axis.line.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = 'bold'),
        axis.text = element_text(size = 14,face = 'bold'),
        axis.title = element_text(face = 'bold', size = 14))+
  xlab(NULL) + ylab('Log10(expression)')


H2_2<- plot_genes_violin(cds_subset_GCG, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        plot.title = element_text(face = 'bold'),
        axis.text = element_text(size = 14,face = 'bold'),
        axis.title = element_text(face = 'bold', size = 14))+
  xlab('Cluster') + 
  ylab('Log10(expression)')


H2_3<- plot_genes_violin(cds_subset_PPY, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        plot.title = element_text(face = 'bold', size = 14),
        axis.text = element_text(size = 14,face = 'bold'),
        axis.title = element_text(face = 'bold', size = 14)) + 
  xlab('Cluster') + 
  ylab('Log10(expression)')


H2_4<- plot_genes_violin(cds_subset_SST, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(), axis.line.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 14,face = 'bold'),
        plot.title = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold', size = 14))+
  xlab(NULL) + ylab('Log10(expression)')

##Figure 2H-K
(H2_1/H2_2|H2_4/H2_3)

##Figure 2I

umap_cds <- reducedDims(cds)$UMAP

data.integrated[['umap_MNN']] <-
  CreateDimReducObject(embeddings = umap_cds,
                       key = "umap_MNN_")

DefaultAssay(data.integrated) <- 'SCT'

I2_1 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('INS', 'GCG'),
  pt.size = 0.5,
  order = T,
  blend = T
)

I2_2 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('INS', 'SST'),
  pt.size = 0.5,
  order = T,
  blend = T
)

I2_3 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('INS', 'PPY'),
  pt.size = 0.5,
  order = T,
  blend = T
)

(I2_1/I2_2/I2_3)

################################################################################
################################################################################

##Figure 2J
J2 <- df %>% drop_na(Hba1c) %>% 
  ggplot(., aes(monocle_clusters,Hba1c, fill = monocle_clusters))+
  geom_boxplot()+
  ggpubr::theme_pubr(border = FALSE, legend = 'none')+ 
  scale_x_discrete(limits = desired_order)+
  labs(x = "Cluster", y = "HbA1c (%)")+
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 16),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 10),
    legend.title  = element_text(face = 'bold', size = 14)
  )

##Figure 2K
K2 <- df %>% drop_na(Hba1c) %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
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
  ggplot(aes(x = monocle_clusters, 
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
    size = 5, 
    vjust = 2,
    fontface = 'bold',
    position = position_dodge(0.8)
  ) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median HbA1c (%)")+ 
  guides(fill = guide_legend(title = "HbA1c class"))+
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 16),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 14),
    legend.title  = element_text(face = 'bold', size = 14)
  )

##Figure 2L
L2 <-
  plot_cells(cds,show_trajectory_graph = F,
             color_cells_by = 'Hba1c_class2',
             label_cell_groups = F,
             group_cells_by = 'cluster',
             cell_size = 1
  ) + ggpubr::clean_theme() + NoLegend()+
  scale_color_manual(values = c('#F8766D','#FEAF16FF','#619CFF'))+
  guides(color = guide_legend(title = "HbA1c class"))

(J2/K2)
L2

######################################################################################
################################# Slingshot ##########################################
######################################################################################

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
colData(sim)$Hba1c_class = df$Hba1c_class
colData(sim)$Hba1c_class2 = df$Hba1c_class2


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

##Figure 2M
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

M2 <- trajectory_plot+
  ggpubr::clean_theme()

# print(trajectory_plot+
#        ggpubr::theme_pubr(legend = 'right')+
#         theme(
#           axis.title.y = element_text(face = 'bold', size = 16),
#           axis.text     = element_text(face = 'bold', size = 16),
#           axis.title.x = element_text(face = 'bold', size = 16),
#           legend.text   = element_text(face = 'bold', size = 14),
#           legend.title  = element_text(face = 'bold', size = 14)
#         ))
M2

##Figure 2N
colData(cds)$Lineage1 <- sim@colData$slingPseudotime_1
colData(cds)$Lineage2 <- sim@colData$slingPseudotime_2
colData(cds)$Lineage3 <- sim@colData$slingPseudotime_3

N2_1 <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = 'Lineage1',
  cell_size = 1,
  group_label_size = 10
)+ggpubr::clean_theme()+
  theme(
    legend.text   = element_text(face = 'bold', size = 14),
    legend.title  = element_text(face = 'bold', size = 14)
  )

N2_2 <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = 'Lineage2',
  cell_size = 1,
  group_label_size = 10
)+ggpubr::clean_theme()+
  theme(
    legend.text   = element_text(face = 'bold', size = 14),
    legend.title  = element_text(face = 'bold', size = 14)
  )

N2_3 <- plot_cells(
  cds,
  show_trajectory_graph = FALSE,
  color_cells_by = 'Lineage3',
  cell_size = 1,
  group_label_size = 10
)+ggpubr::clean_theme()+
  theme(
    legend.text   = element_text(face = 'bold', size = 14),
    legend.title  = element_text(face = 'bold', size = 14)
  )

M2 + N2_1 + N2_2 + N2_3

#######################################################################################
####################### Lineage subsets generation #################################
#######################################################################################

lin1_obj <- subset(data.integrated, cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage1)])
dim(lin1_obj)
saveRDS(lin1_obj, '20250601_Beta_T2D_lineage1.rds')
lin1_obj <- SetIdent(lin1_obj, value = 'assigned_label')
DefaultAssay(lin1_obj) <- 'RNA'
lin1_obj %>% 
  scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
  write.table(., 'Lineage1_sling_perc_expr.tsv',row.names = T,
              col.names = NA, sep = '\t')


lin2_obj <- subset(data.integrated, 
                   cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage2)])
dim(lin2_obj)
saveRDS(lin2_obj, '20250601_Beta_T2D_lineage2.rds')
lin2_obj <- SetIdent(lin2_obj, value = 'assigned_label')
DefaultAssay(lin2_obj) <- 'RNA'
lin2_obj %>% 
  scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
  write.table(., 'Lineage2_sling_perc_expr.tsv',
              row.names = T,col.names = NA, sep = '\t')


lin3_obj <- subset(data.integrated, 
                   cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage3)])
dim(lin3_obj)
saveRDS(lin3_obj, '20250601_Beta_T2D_lineage3.rds')
lin3_obj <- SetIdent(lin3_obj, value = 'assigned_label')
DefaultAssay(lin3_obj) <- 'RNA'
lin3_obj %>% 
  scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
  write.table(., 'Lineage3_sling_perc_expr.tsv',
              row.names = T,col.names = NA, sep = '\t')

################################################################################
################################################################################

# Extract expression data
exp <- data.integrated@assays$RNA@data
exp <- t(exp)

sce <- readRDS('sce_all.rds')

# Extract pseudotime values and cell weights
pseudotime <- as.data.frame(colData(sce)$crv)
# colnames(pseudotime)[1:3] <- c("Lineage1", "Lineage2", "Lineage3")
# Unified pseudotime assignment
pseudotime$UnifiedPseudotime <- rep(NA, nrow(pseudotime))

for (jj in 1:3) {
  lineage_column <- paste0("cellWeights.Lineage", jj)
  pseudotime_column <- paste0("pseudotime.Lineage", jj)
  pseudotime$UnifiedPseudotime[pseudotime[, lineage_column] == 1] <- pseudotime[pseudotime[, lineage_column] == 1, pseudotime_column]
}

# Normalize pseudotime
pseudotime$AdjustedPseudotime <- pseudotime$UnifiedPseudotime

for (lineage in 1:3) {
  lineage_mask <- pseudotime[, paste0("cellWeights.Lineage", lineage)] == 1
  min_pseudotime <- min(pseudotime$UnifiedPseudotime[lineage_mask], na.rm = TRUE)
  max_pseudotime <- max(pseudotime$UnifiedPseudotime[lineage_mask], na.rm = TRUE)
  pseudotime$AdjustedPseudotime[lineage_mask] <- 
    (pseudotime$UnifiedPseudotime[lineage_mask] - min_pseudotime) / (max_pseudotime - min_pseudotime)
}

# Assign lineage membership
pseudotime$Membership <- apply(pseudotime[, c("cellWeights.Lineage1", "cellWeights.Lineage2", "cellWeights.Lineage3")], 1, function(weights) {
  lineages <- c("Lineage1", "Lineage2", "Lineage3")
  paste0(lineages[weights > 0], collapse = ", ")
})

# Assign dominant lineage based on the maximum weight
pseudotime$DominantLineage <- apply(pseudotime[, c("cellWeights.Lineage1", "cellWeights.Lineage2", "cellWeights.Lineage3")], 1, function(weights) {
  lineages <- c("Lineage1", "Lineage2", "Lineage3")
  if (all(weights == 0)) return(NA)  # Handle cells with no lineage assignment
  lineages[which.max(weights)]  # Assign the lineage with the highest weight
})


######### Plotting

# Define custom colors for points and smoothers
point_colors <- c(
  "Lineage1" = "#89CFF0",  # Light Blue for Lineage1
  "Lineage2" = "#FF9999",  # Soft Red for Lineage2
  "Lineage3" = "#FFB347"   # Lighter Orange for Lineage3
)

smoother_colors <-
  c("Lineage1" = "blue",
    "Lineage2" = "red",
    "Lineage3" = "#d89000")



# Initialize plot list
plt_list <- list()

genes <- 'INS'

# Loop through the gene(s) of interest
for (gene in genes) {  # Replace with other genes if needed
  print(gene)
  # Prepare plt_data with the updated pseudotime and membership
  plt_data <- data.frame(
    AdjustedPseudotime = pseudotime$AdjustedPseudotime,
    Membership = pseudotime$Membership,
    DominantLineage = pseudotime$DominantLineage,
    gene = exp[,gene ]  # Extract gene expression values for the current gene
  )
  
  # Generate plot
  p <- ggplot(data = plt_data, aes(x = AdjustedPseudotime, y = gene)) +
    # Add points, colored by dominant lineage
    geom_point(aes(color = DominantLineage), size = 0.6, alpha = 0.6) +
    
    # Add smoothers for each lineage
    geom_smooth(data = subset(plt_data, grepl("Lineage1", Membership)),
                aes(group = "Lineage1"), size = 2, 
                color = smoother_colors["Lineage1"],
                se = FALSE, method = 'lm', formula = y ~ splines::ns(x, df = 2)) +
    geom_smooth(data = subset(plt_data, grepl("Lineage2", Membership)),
                aes(group = "Lineage2"), size = 2, 
                color = smoother_colors["Lineage2"],
                se = FALSE, method = 'lm', formula = y ~ splines::ns(x, df = 2)) +
    geom_smooth(data = subset(plt_data, grepl("Lineage3", Membership)),
                aes(group = "Lineage3"), size = 2, 
                color = smoother_colors["Lineage3"],
                se = FALSE, method = 'lm', formula = y ~ splines::ns(x, df = 2)) +
    # ylim(0,l[i]) +
    # ylim(0,2) +
    # Customize colors and appearance
    scale_color_manual(name = "Lineages", values = point_colors) +
    labs(x = "Pseudotime", y = gene) +
    theme_pubr(legend = 'right') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(face = 'bold', size = 16),
          axis.text     = element_text(face = 'bold', size = 18),
          axis.title.x = element_text(face = 'bold', size = 16),
          legend.text   = element_text(face = 'bold', size = 14),
          legend.title  = element_text(face = 'bold', size = 14)) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  # Save the plot
  plt_list[[gene]] <- p
}

plt_list

###############################################################################
# sessionInfo()
###############################################################################
