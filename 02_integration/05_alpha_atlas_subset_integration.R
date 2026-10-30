################################################################################
# Alpha-cell object: subset from atlas, SCT integration by study and Monocle3
# clustering (no trajectory inference; the Slingshot block is retained commented
# out as in the original). Produces alpha.integrated.rds and cds_alpha.rds,
# read by Figure2_S3/FigureS3_clean.R.
# Original code; only file paths and the lines marked #! were changed.
################################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)  #! restored: %>%, mutate_if, group_by, left_join, drop_na and ggplot are used below
  library(patchwork)
  library(reticulate)
  library(monocle3)
  library(slingshot)
})

set.seed(123)

reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')

palette <- c("#5A5156FF", "#E4E1E3FF", "#F6222EFF" ,
             "#FE00FAFF", "#16FF32FF" ,"#3283FEFF",
             "#FEAF16FF", "#B00068FF", "#1CFFCEFF",
             "#90AD1CFF", "#2ED9FFFF" ,"#DEA0DFFF",
             "#AA0EDFFF", "#F8A19FFF", "#325A9BFF" ,
             "#C4451CFF", "#1C8356FF" ,"#85660DFF",
             "#B10DA1FF", "#FBE426FF", "#1CBE4FFF")

cores = 2
resolution = 0.00001

################################################################################
################################################################################

## Step 0: subset alpha cells (ND and T2D) from the classified atlas
## (this block was commented out in the original run; the output Alpha_T2D_raw.rds is read below)
data <- readRDS('scRNA_integrated_harmony_classified.rds')
dim(data)

data <- subset(data, subset = Condition %in% c('ND', 'T2D'))
Idents(data) <- 'assigned_label'  #! select the atlas cell-type annotation
data <- subset(data, idents = 'Alpha')
dim(data)

saveRDS(data, 'Alpha_T2D_raw.rds')

################################################################################
######################### scTransform integration ##############################
################################################################################

data <-
  readRDS(
    'Alpha_T2D_raw.rds'
  )
dim(data)

# data <- subset(data, subset = nCount_RNA < 1500000 &
#                  nFeature_RNA < 7500)

# VlnPlot(data,features = c('nCount_RNA','nFeature_RNA','logCounts') )

DefaultAssay(data) <- 'RNA'  #! SCTransform uses raw RNA measurements (input has 18,829 genes = RNA assay)
data.list <- SplitObject(data, split.by = "study")
rm(data)
gc()

lapply(data.list, dim)
# $GSE97655
# [1] 18829    32

## Exclude datasets with insufficient cells ####
data.list <- data.list[sapply(data.list, function(x) dim(x)[2] > 200)]
gc()
"excluding datasets with fewer than 200 cells  / "

# perform SCTransform normalization
data.list <- lapply(X = data.list,
                    FUN = SCTransform,
                    vars.to.regress = c("nCount_RNA",
                                        "nFeature_RNA",
                                        "percent.mt"),
                    return.only.var.genes = F)
gc()

# select integration features
features <- SelectIntegrationFeatures(data.list)

data.list <- PrepSCTIntegration(object.list = data.list,
                                anchor.features = features)
gc()

anchors <-
  FindIntegrationAnchors(
    object.list = data.list,
    k.filter = NA,
    anchor.features = features,
    dims = 1:30
  )

rm(data.list)
gc()

# saveRDS(anchors, 'anchors.rds')

data.integrated <- IntegrateData(anchorset = anchors, dims = 1:30,
                                 normalization.method = 'SCT')

DefaultAssay(data.integrated)

data.integrated <- RunPCA(data.integrated, verbose = FALSE)

saveRDS(data.integrated,'alpha.integrated.rds')

###################################################################################
########################### Monocle preprocess ####################################
###################################################################################

data.integrated <-
  readRDS('alpha.integrated.rds')
dim(data.integrated)

print('######################## extract data for cds ######################################')

expr <- as(as.matrix(data.integrated@assays$SCT@counts), 'sparseMatrix')
gene_metadata <- data.frame(gene_short_name = row.names(expr),
                            row.names = row.names(expr))
cell_metadata <- as.data.frame(data.integrated@meta.data)

cds <- new_cell_data_set(expr,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)


colnames(pData(cds))[10] <- 'sample_name2'
fData(cds)$gene_short_name <- rownames(fData(cds))

print('################################ assign paritions #########################################')
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

meta <- pData(cds)
pData(cds) <- meta

rm(expr)
gc()

reducedDims(cds)$PCA <- data.integrated@reductions$pca@cell.embeddings
reducedDims(cds)$PCA <- reducedDims(cds)$PCA[,1:10]

meta <- as.data.frame(colData(cds))
meta <- mutate_if(meta, is.character, as.factor)

cds <- reduce_dimension(cds, cores = cores, umap.fast_sgd = F,
                        reduction_method = 'UMAP')
cds <- cluster_cells(cds, verbose = T, resolution = resolution, 
                     cluster_method = 'leiden' )

colData(cds)$monocle_clusters <- cds@clusters@listData$UMAP$clusters
data.integrated@meta.data$monocle_clusters <- pData(cds)$monocle_clusters 

df <- data.integrated@meta.data %>% as.data.frame()
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

###### adding HbA1cmeta data to cds object ######
colData(cds)$mean_value <- df$mean_value
colData(cds)$median_value <- df$median_value
colData(cds)$Hba1c_class <- df$Hba1c_class
colData(cds)$Hba1c_class2 <- df$Hba1c_class2

rownames(df) <- rownames(data.integrated@meta.data)

saveRDS(cds, 'cds_alpha.rds')   #! added: cds_alpha.rds is read by FigureS3_clean.R but was not saved here

###### adding HbA1cmeta data to seurat object ######
data.integrated@meta.data <- df

# write.table(as.data.frame(colData(cds)),col.names = NA,row.names = T,
#             sep = '\t','mncl_coldata.tsv')

# df <- read.table('mncl_coldata.tsv',header = T,row.names = 1)

df$monocle_clusters <- as.character(df$monocle_clusters)

desired_order <- df %>%
  group_by( monocle_clusters,Condition) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>%
  group_by(monocle_clusters) %>%
  mutate(per = round(100 *n/sum(n),1)) %>% 
  filter(Condition == 'T2D') %>%
  arrange(per) %>% pull(monocle_clusters)
desired_order

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

# ggplot(numbers_df, aes(x = factor(Clusters, levels = desired_order), y = n_Samples)) + 
#   geom_bar(stat = 'identity') +
#   ggpubr::theme_pubr(border = FALSE, legend = 'top', base_size = 10) +
#   xlab("Monocle Clusters") + 
#   ylab("Number of Unique Samples per Cluster") +
#   scale_x_discrete(limits = desired_order)

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
                 group_label_size = 10)+
  ggpubr::clean_theme()

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
                 color_cells_by = 'study',
                 label_cell_groups=FALSE, cell_size = 1,
                 group_cells_by = 'cluster', 
                 group_label_size = 4)+
  ggpubr::clean_theme()+NoLegend()+
  scale_color_manual(values = palette[-c(12,13,14)])

E2 <- df %>%
  group_by( monocle_clusters,study) %>%
  dplyr::count() %>%
  mutate(monocle_clusters = as.character(monocle_clusters)) %>% 
  group_by(monocle_clusters) %>% 
  mutate(per = paste0(round(100 *n/sum(n),1),'%')) %>%
  ggplot(aes(x = monocle_clusters, y = n, fill = study)) +
  geom_bar(stat = 'identity', position = "fill", color = 'black', width = 0.8) +
  ggpubr::theme_pubr(border = FALSE, legend = 'right',base_size = 10) + 
  scale_x_discrete(limits = desired_order) + 
  # geom_text(aes( label = per),position = position_fill(vjust = 0.5)) +
  ylab('Proportion') +xlab(NULL)+
  scale_fill_manual(values = palette[-c(12,13,14)])+
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 18),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 12),
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
  ggpubr::theme_pubr(border = FALSE, legend = 'right',base_size = 10) + 
  scale_x_discrete(limits = desired_order)+
  ylab('Proportion') +xlab('Cluster') + 
  geom_text(aes( label = per),size = 5, fontface = 'bold',
            position = position_fill(vjust = 0.5))+
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 18),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 12),
    legend.title  = element_text(face = 'bold', size = 14)
  )

## Figures 2D-G ####
# (D2/E2)|(F2/G2)
E2/G2

#### Figure 2H ####
H2_1<- plot_genes_violin(cds_subset_INS, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(), axis.line.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(face = 'bold'))+
  xlab(NULL)


H2_2<- plot_genes_violin(cds_subset_GCG, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        plot.title = element_text(face = 'bold', size = 14),
        axis.text = element_text(size = 16,face = 'bold'),
        axis.title = element_text(face = 'bold', size = 14))+
  xlab('Cluster') + 
  ylab('Log10(expression)')


H2_3<- plot_genes_violin(cds_subset_PPY, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(),
        plot.title = element_text(face = 'bold', size = 14),
        axis.text = element_text(size = 16,face = 'bold'),
        axis.title = element_text(face = 'bold', size = 14))+
  xlab('Cluster') + 
  ylab('Log10(expression)')


H2_4<- plot_genes_violin(cds_subset_SST, group_cells_by='monocle_clusters')+ 
  scale_x_discrete(limits = desired_order)+ggpubr::theme_pubr()+
  theme(strip.background = element_blank(), axis.line.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(face = 'bold'))+
  xlab(NULL)
##Figure 2H-K
# (H2_1/H2_2|H2_4/H2_3)
H2_2/H2_3

##Figure 2I

umap_cds <- reducedDims(cds)$UMAP

data.integrated[['umap_MNN']] <-
  CreateDimReducObject(embeddings = umap_cds,
                       key = "umap_MNN_")

DefaultAssay(data.integrated) <- 'SCT'

I2_1 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('GCG', 'INS'),
  pt.size = 0.5,
  order = T,
  blend = T
)
I2_2 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('GCG', 'SST'),
  pt.size = 0.5,
  order = T,
  blend = T
)
I2_3 <- FeaturePlot(
  data.integrated,
  reduction = 'umap_MNN',
  features = c('GCG', 'PPY'),
  pt.size = 0.5,
  order = T,
  blend = T
)
(I2_1/I2_2/I2_3)

##Figure 2J
J2 <- df %>% drop_na(Hba1c) %>% 
  ggplot(., aes(monocle_clusters,Hba1c, fill = monocle_clusters))+
  geom_boxplot()+
  ggpubr::theme_pubr(border = FALSE, legend = 'none')+ 
  scale_x_discrete(limits = desired_order)+
  labs(x = "Cluster", y = "HbA1c (%)") +
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
  ggplot(aes(x = monocle_clusters, y = median_value, fill = Hba1c_class)) +
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
  geom_text(
    aes(label = paste0(round(median_value, 1), '%')),
    size = 7,
    vjust = 2,
    fontface = 'bold',
    position = position_dodge(0.8)
  ) +
  ggpubr::theme_pubr(border = FALSE, legend = 'top',base_size = 10) + 
  scale_x_discrete(limits = desired_order) +
  labs(x = "Cluster", y = "Median HbA1c(%)")+ 
  guides(fill = guide_legend(title = "HbA1c class"))+
  theme(
    axis.title.y = element_text(face = 'bold', size = 16),
    axis.text     = element_text(face = 'bold', size = 16),
    axis.title.x = element_text(face = 'bold', size = 16),
    legend.text   = element_text(face = 'bold', size = 10),
    legend.title  = element_text(face = 'bold', size = 14)
  )

##Figure 2L
L2 <-
  plot_cells(cds,show_trajectory_graph = F,
             color_cells_by = 'Hba1c_class2',
             label_cell_groups = F,cell_size = 1,
             group_cells_by = 'cluster'
  ) + ggpubr::clean_theme() + NoLegend()+
  guides(color = guide_legend(title = "HbA1c class"))

(J2/K2)/L2


######################################################################################
################################# Slingshot ##########################################
######################################################################################
# 
# sim <-
#   SingleCellExperiment(assays =
#                          List(counts =
#                                 as.matrix(data.integrated@assays$SCT@counts)))
# 
# umap_cds <- reducedDims(cds)$UMAP
# reducedDims(sim) = SimpleList(UMAP = umap_cds)
# 
# # metadata
# colData(sim)$study = df$study
# colData(sim)$monocle_clusters =  df$monocle_clusters
# colData(sim)$Condition = df$Condition
# colData(sim)$Hba1c_class = df$Hba1c_class
# colData(sim)$Hba1c_class2 = df$Hba1c_class2
# 
# 
# sim <- slingshot(sim, 
#                  clusterLabels = 'monocle_clusters',  
#                  reducedDim = 'UMAP',  
#                  start.clus= "3",  
#                  end.clus = NULL   
# )     
# 
# lineages <- slingshot::slingCurves(SlingshotDataSet(sim))
# 
# # Assuming 'lineages' is already obtained from slingshot::slingCurves
# lineage_data <- lapply(names(lineages), function(lineage_name) {
#   curve_df <- as.data.frame(lineages[[lineage_name]]$s)
#   colnames(curve_df) <- c("Dim-1", "Dim-2")  # Rename columns based on your structure
#   curve_df$lineage <- lineage_name
#   curve_df$path_id <- paste(lineage_name, seq_len(nrow(curve_df)), sep="_")
#   return(curve_df)
# })
# 
# # Combine all lineage data frames into one
# combined_lineages <- do.call(rbind, lineage_data)
# 
# 
# ##Figure 2M
# umap_plot <- plot_cells(
#   cds,
#   show_trajectory_graph = FALSE,
#   color_cells_by = 'cluster',
#   label_groups_by_cluster = F,
#   cell_size = 1,label_cell_groups = T,
#   group_cells_by = 'cluster',
#   group_label_size = 10
# )
# 
# # Create a trajectory plot
# trajectory_plot <- umap_plot + geom_path(
#   data = combined_lineages,
#   aes(
#     x = `Dim-1`,
#     y = `Dim-2`,
#     color = lineage,
#     group = lineage
#     
#   ),
#   size = 2
# )+ scale_color_manual(values = c('1'="#F8766D",'2'= "#C49A00",'3'= "#53B400",
#                                  '4'= "#00C094",'5'= "#00B6EB" ,'6'="#A58AFF" ,
#                                  '7'="#FB61D7",'Lineage2'='red',
#                                  'Lineage1'='blue','Lineage3'="#D89000", 
#                                  'Lineage4'="#39B600", 'Lineage5'="#9590FF"))
# 
# M2 <- trajectory_plot+
#   ggpubr::clean_theme()
# # print(trajectory_plot+
# #        ggpubr::theme_pubr(legend = 'right'))
# 
# 
# ##Figure 2N
# colData(cds)$Lineage1 <- sim@colData$slingPseudotime_1
# colData(cds)$Lineage2 <- sim@colData$slingPseudotime_2
# colData(cds)$Lineage3 <- sim@colData$slingPseudotime_3
# 
# N2_1 <- plot_cells(
#   cds,
#   show_trajectory_graph = FALSE,
#   color_cells_by = 'Lineage1',
#   cell_size = 1,
#   group_label_size = 10
# )+ggpubr::clean_theme()
# 
# N2_2 <- plot_cells(
#   cds,
#   show_trajectory_graph = FALSE,
#   color_cells_by = 'Lineage2',
#   cell_size = 1,
#   group_label_size = 10
# )+ggpubr::clean_theme()
# 
# N2_3 <- plot_cells(
#   cds,
#   show_trajectory_graph = FALSE,
#   color_cells_by = 'Lineage3',
#   cell_size = 1,
#   group_label_size = 10
# )+ggpubr::clean_theme()
# 
# M2 + N2_1 + N2_2 + N2_3
# 
# 
# #######################################################################################
# ####################### Lineage pseudobulk generation #################################
# #######################################################################################
# 
# ## Calculate percentage expressed in lineage 1 cells
# 
# data.integrated@meta.data$Lineage_1 <- sim@colData$slingPseudotime_1
# data.integrated@meta.data$Lineage_2 <- sim@colData$slingPseudotime_2
# data.integrated@meta.data$Lineage_3 <- sim@colData$slingPseudotime_3
# 
# 
# 
# lin1_obj <- subset(data.integrated, cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage1)])
# saveRDS(lin1_obj, 'alpha_T2D_lineage1.rds')
# lin1_obj <- SetIdent(lin1_obj, value = 'assigned_label')
# DefaultAssay(lin1_obj) <- 'RNA'
# lin1_obj %>% 
#   scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
#   write.table(., 'Lineage1_sling_perc_expr.tsv',row.names = T,
#               col.names = NA, sep = '\t')
# 
# 
# 
# lin2_obj <- subset(data.integrated, 
#                    cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage2)])
# saveRDS(lin2_obj, 'alpha_T2D_lineage2.rds')
# lin2_obj <- SetIdent(lin2_obj, value = 'assigned_label')
# DefaultAssay(lin2_obj) <- 'RNA'
# lin2_obj %>% 
#   scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
#   write.table(., 'Lineage2_sling_perc_expr.tsv',
#               row.names = T,col.names = NA, sep = '\t')
# 
# 
# 
# lin3_obj <- subset(data.integrated, 
#                    cells = rownames(colData(cds))[!is.na(colData(cds)$Lineage3)])
# saveRDS(lin3_obj, 'alpha_T2D_lineage3.rds')
# lin3_obj <- SetIdent(lin3_obj, value = 'assigned_label')
# DefaultAssay(lin3_obj) <- 'RNA'
# lin3_obj %>% 
#   scCustomize::Percent_Expressing(., features = rownames(.)) %>% 
#   write.table(., 'Lineage3_sling_perc_expr.tsv',
#               row.names = T,col.names = NA, sep = '\t')


#######################################################################################
#######################################################################################
#######################################################################################

sessionInfo()





