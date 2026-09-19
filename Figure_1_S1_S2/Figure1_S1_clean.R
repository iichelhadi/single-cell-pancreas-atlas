###############################################################################
# Project: Pancreas scRNA atlas — figures 1 and S1
# File: Figure1_S1.R
# Author: Elhadi Iich
# Date: 2026-09-07
#
# Description:
#   Reproducible plotting script that generates main and supplemental figures
#   from an integrated Seurat object.
#
# Input:
#   - scRNA_integrated_harmony_classified.rds
#   - DEGs_allmarkers_assigned_label.tsv
#   - DEGs_allmarkers_res_0.2_ordered_pcs.tsv
#
# Dependencies:
#   Seurat, patchwork, tidyverse, doParallel, ggsci, scCustomize, ggpubr, Cairo
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(doParallel)
  library(ggsci)
  library(scCustomize)
  library(ggpubr)
  library(Cairo)
})
options(bitmapType = 'cairo')

palette <- c("#5A5156FF", "#E4E1E3FF", "#F6222EFF" ,"#FE00FAFF", "#16FF32FF" ,"#3283FEFF",
             "#FEAF16FF", "#B00068FF", "#1CFFCEFF" ,"#90AD1CFF", "#2ED9FFFF" ,"black",
             "#AA0DFEFF", "#F8A19FFF", "#325A9BFF" ,"#C4451CFF", "#1C8356FF" ,"#85660DFF",
             "#B10DA1FF", "#FBE426FF", "#1CBE4FFF")

################################################################################
# Load data
################################################################################

merged_seurat <- readRDS('scRNA_integrated_harmony_classified.rds')

################################################################################
# Figure 1A
################################################################################

DimPlot(
  merged_seurat,
  label = F,
  raster = T,cols = palette,
  reduction = 'umap_harmony'
)

################################################################################
# Figure 1B
################################################################################

myfeatures <- c('INS','GCG','SST','PPY',
                'GHRL','HLA-DRA','PRSS1',
                'PLVAP','KRT7','COL1A1','MKI67')
plot_list <- list()
for (i in myfeatures) {
  plot_list[[i]] <- Plot_Density_Custom(seurat_object = merged_seurat, features = i)+
    NoLegend() + NoAxes() + NoGrid()
}

# pdf('20241030_features_harmony_20_20.pdf', width = 10, height = 20)
cowplot::plot_grid(ncol = 2, plotlist = plot_list)

################################################################################
# Figure 1C
################################################################################

merged_seurat@meta.data %>% as.data.frame() %>% 
  group_by(assigned_label) %>% dplyr::count() %>% 
  ggplot(., aes(x = assigned_label, y = n, fill = assigned_label)) + 
  geom_bar(stat = 'identity',color = 'black',width = 0.8)+ theme_pubr(border = F) +
  ylab('Count')+xlab('')+ Seurat::NoLegend()+
  geom_text(aes(label=n),size= 5, fontface = 'bold',
            position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(
    values = c(
      "Alpha" = "#5A5156FF",
      "Ductal" = "#E4E1E3FF",
      "Acinar" = "#F6222EFF",
      "Beta" = "#FE00FAFF",
      "Pancreatic Stellate" = "#16FF32FF",
      "Endothelial" = "#3283FEFF",
      "Delta" = "#FEAF16FF",
      "Pancreatic Polypeptide" = "#B00068FF",
      "MHC Class II" = "#1CFFCEFF",
      'Epsilon' = "#90AD1CFF",
      'Mast' = "#1C8356FF",
      'Bihormonal' = "#C4451CFF",
      'Cycling' = 'black' )) +
  theme(legend.position = 'none',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        axis.text.x = element_text(face = 'bold', size = 16, angle = 45, 
                                   hjust = 1,vjust = 1),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1))

################################################################################
# Figure 1D
################################################################################

markers <- read_tsv('DEGs_allmarkers_assigned_label.tsv')

top_markers <- markers %>% group_by(cluster) %>% top_n(3, avg_log2FC) %>% pull(gene)

top_markers_sel <- c(top_markers,'DLK1')


palette <- c("#E4E1E3FF","#F6222EFF","#16FF32FF",
             "#5A5156FF","#1CFFCEFF","#3283FEFF",
             "#FE00FAFF",'black',"#C4451CFF",
             "#FEAF16FF", "#B00068FF" ,"#90AD1CFF")


# CairoPNG('S1D.png', width = 5000, height = 1500,res = 300)
Clustered_DotPlot(
  seurat_object = merged_seurat,
  features = top_markers_sel,
  k = 12,
  show_parent_dend_line = F,row_label_size = 14,
  exp_color_min = -2,
  colors_use_idents = palette,
  exp_color_max = 2,
  flip = T,
  x_lab_rotate = 90,
  colors_use_exp = CustomPalette(low = "blue", 
                                 high = "red", 
                                 mid = 'white', k = 50))
dev.off()

################################################################################
# Figure S1B
################################################################################

######################### Figure S1B
df <- as.data.frame(merged_seurat@meta.data)

p1 <- ggplot(df,aes(Study, percent.mt, fill = Study)) +
  geom_violin(scale = "width") + 
  geom_boxplot(width = 0.2,fill = 'white')+
  scale_fill_manual(values = palette) +
  theme_pubr() + NoLegend() +
  theme(legend.position = 'none',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1))+
  xlab(NULL)

p2 <- ggplot(df,aes(Study, nFeature_RNA, fill = Study)) +
  geom_violin(scale = "width") + 
  geom_boxplot(width = 0.2,fill = 'white')+
  scale_fill_manual(values = palette) + 
  theme_pubr() +  NoLegend()+
  theme(legend.position = 'none',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1))+
  xlab(NULL)

p3 <- ggplot(df,aes(Study, nCount_RNA, fill = Study)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.2,fill = 'white')+
  scale_fill_manual(values = palette) +
  theme_pubr() + NoLegend()+
  theme(legend.position = 'none',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        axis.text.x = element_text(angle = 60, vjust = 1,hjust = 1),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1))+
  xlab('Study')

(p1/p2/p3)
dev.off()

################################################################################
# Figure S1D
################################################################################

################## Figure S1D

# pdf('S1D.pmerged_seurat',width = 9)
DimPlot(merged_seurat,
        reduction = 'umap_harmony',
        group.by = 'RNA_snn_res.0.2',
        raster = T,cols = palette,
        label.size = 9,
        label = T)
dev.off()

#########################################################################
# Figure S1E
#########################################################################

################## Figure S1E
merged_seurat$RNA_snn_res.0.2 <- factor(merged_seurat$RNA_snn_res.0.2, 
                                     levels = as.character(1:12))

VlnPlot(merged_seurat, features = 'log10GenesPerUMI',
        cols = palette, pt.size = 0, group.by = 'RNA_snn_res.0.2') + 
  scale_x_discrete(limits = as.character(c(1:12)))+
  xlab('Cluster') +
  theme(legend.position = 'none',
        plot.title = element_text(face = 'bold', size = 24),
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        axis.text.x = element_text(face = 'bold', size = 16, 
                                   angle = 0,hjust = 0.5,vjust = 1),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1))

#############################################################################
# Figure S1F
#############################################################################

numbers <- list()
for (i in 1:length(unique(merged_seurat$RNA_snn_res.0.2))){
  numbers <- unlist(append(numbers,print(length(table(merged_seurat$Study[merged_seurat$RNA_snn_res.0.2 == i])))))
}

numbers_df <- data.frame(Clusters = 1:length(unique(merged_seurat$RNA_snn_res.0.2)), n_Studies = numbers)
ggplot(numbers_df, aes(x = factor(Clusters),y = n_Studies, fill = factor(Clusters)))+
  geom_bar(stat="identity",color = 'black',width = 0.8)+
  scale_fill_manual(values=palette)+theme_pubr(border = F)+ylim(0,20)+
  ylab('Number of studies') +xlab('Cluster')+
  theme(legend.position = 'none',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1)) -> S1F

# pdf('S1E.pdf',width = 6,height = 4)
S1F
dev.off()

##################################################################
# Figures S1G, S1H, S1K and stacked panels
##################################################################

######## Figure S1G

DimPlot(merged_seurat,
        reduction = 'umap_harmony',
        group.by = 'Study',
        raster = T,
        cols = palette,
        label.size = 9,
        label = F) + theme(plot.title = element_text(size = 30))
dev.off()

########### Figure S1H
DimPlot(merged_seurat,
        reduction = 'umap_harmony',
        group.by = 'Platform',
        raster = T,cols = palette,
        label.size = 9,
        label = F)+ theme(plot.title = element_text(size = 30))
dev.off()

########### Figure S1K

DimPlot(merged_seurat,
        reduction = 'umap_harmony',
        group.by = 'Condition',
        raster = T,cols = palette[-4],
        label.size = 9,
        label = F)+ theme(plot.title = element_text(size = 30))
dev.off()


######### Figure S1H
merged_seurat@meta.data %>% as.data.frame() %>%
  group_by(RNA_snn_res.0.2,Study)%>%
  dplyr::count() %>%
  ggplot(., aes(x = factor(RNA_snn_res.0.2), y = n, fill = Study)) +
  geom_bar(stat="identity",position="fill",color = 'black',width = 0.8)+
  theme_pubr(border = F,legend = 'right') +
  ylab('Proportion') + scale_fill_manual(values=palette) + 
  xlab('Cluster')+
  theme(legend.position = 'right',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1)) -> S1H

######### Figure S1J
merged_seurat@meta.data %>% as.data.frame() %>%
  group_by(RNA_snn_res.0.2,Platform)%>%
  dplyr::count() %>%
  ggplot(., aes(x = factor(RNA_snn_res.0.2), y = n, fill = Platform)) +
  geom_bar(stat="identity",position="fill",color = 'black',width = 0.8)+ 
  theme_pubr(legend = 'right')+
  ylab('Proportion')+scale_fill_manual(values=palette)+ 
  xlab('Cluster') +
  theme(legend.position = 'right',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1)) -> S1J

######### Figure S1H
merged_seurat@meta.data %>% as.data.frame() %>%
  group_by(RNA_snn_res.0.2,Condition)%>%
  dplyr::count() %>%
  ggplot(., aes(x = factor(RNA_snn_res.0.2), y = n, fill = Condition)) +
  geom_bar(stat="identity",position="fill",color = 'black',width = 0.8)+ 
  theme_pubr(legend = 'right')+
  ylab('Proportion')+scale_fill_manual(values=palette[-4])+ 
  xlab('Cluster')+
  theme(legend.position = 'right',
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x =  element_text(face = 'bold',size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        legend.title = element_text(face = 'bold', size = 14),
        legend.text = element_text(face = 'bold', size = 14,
                                   hjust = 1,vjust = 1)) -> S1L


# pdf('S1G_S1K.pdf', height = 12, width = 6)
(S1H + xlab(NULL) + 
    guides(x = "none")) / 
  (S1J + xlab(NULL) + 
     guides(x = "none")) / (S1L)
dev.off()



#############

palette2 <- c("#5A5156FF","#90AD1CFF", "#2ED9FFFF",
              'black',"#E4E1E3FF", "#F6222EFF" ,
              "#FE00FAFF", "#16FF32FF" ,"#3283FEFF",
              "#FEAF16FF", "#B00068FF", "#1CFFCEFF",
              "#DEA0DFFF",
              "#AA0EDFFF", "#F8A19FFF", "#325A9BFF" ,
              "#C4451CFF", "#1C8356FF" ,"#85660DFF",
              "#B10DA1FF", "#FBE426FF", "#1CBE4FFF")

VlnPlot_scCustom(
  merged_seurat,
  features = c('nCount_RNA'),
  group.by = 'RNA_snn_res.0.2',
  pt.size = 0,
  colors_use = palette2
) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  scale_x_discrete(limits = factor(1:12))


###############################################################################
# SOURCE DATA EXPORT — Figure 1 and Figure S1
###############################################################################

source_dir <- "../../SourceData/Figure1_S1"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Fig 1A, Fig 1B, Fig S1D, S1G, S1H, S1K
#
# One shared cell-level table:
# UMAP coordinates + annotations + marker expression
###############################################################################

umap_df <- as.data.frame(
  Seurat::Embeddings(merged_seurat, reduction = "umap_harmony")
)

colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell_id <- rownames(umap_df)

meta_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    assigned_label,
    RNA_snn_res.0.2,
    study,
    Platform,
    Condition
  )

expression_source <- Seurat::FetchData(
  merged_seurat,
  vars = myfeatures
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id")

fig1_umap_source <- umap_df %>%
  dplyr::left_join(meta_source, by = "cell_id") %>%
  dplyr::left_join(expression_source, by = "cell_id")

readr::write_tsv(
  fig1_umap_source,
  file.path(
    source_dir,
    "Fig1A_B_S1D_G_H_K_cell_data.tsv.gz"
  )
)


###############################################################################
# Fig 1C — cell-type counts
###############################################################################

fig1c_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::count(
    assigned_label,
    name = "cell_count"
  ) %>%
  dplyr::arrange(
    dplyr::desc(cell_count)
  )

readr::write_tsv(
  fig1c_source,
  file.path(
    source_dir,
    "Fig1C_celltype_counts.tsv"
  )
)


###############################################################################
# Fig 1D — dot plot
#
# Export the actual quantities represented by the dots:
#   pct.exp          = percentage of cells expressing gene
#   avg.exp          = average expression
#   avg.exp.scaled   = scaled average expression used for colour
###############################################################################

# Remove duplicated genes for source-data extraction
top_markers_source <- unique(as.character(top_markers_sel))

# Check what was duplicated
top_markers_sel[duplicated(top_markers_sel)]

fig1d_dot <- Seurat::DotPlot(
  merged_seurat,
  features = top_markers_source,
  group.by = "assigned_label"
)$data %>%
  as.data.frame() %>%
  dplyr::rename(
    cell_type = id,
    gene = features.plot,
    pct_expressing = pct.exp,
    average_expression = avg.exp,
    scaled_average_expression = avg.exp.scaled
  ) %>%
  dplyr::select(
    cell_type,
    gene,
    pct_expressing,
    average_expression,
    scaled_average_expression
  )

readr::write_tsv(
  fig1d_dot,
  file.path(source_dir, "Fig1D_dotplot.tsv")
)


###############################################################################
# Fig S1B and S1E
#
# Cell-level QC distributions.
# Also covers the final nCount_RNA violin at the end of the script.
###############################################################################

figs1_qc_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(
    cell_id,
    study,
    RNA_snn_res.0.2,
    percent.mt,
    nFeature_RNA,
    nCount_RNA,
    log10GenesPerUMI
  )

readr::write_tsv(
  figs1_qc_source,
  file.path(
    source_dir,
    "FigS1B_E_QC_distributions.tsv.gz"
  )
)


###############################################################################
# Fig S1F — number of studies represented in each cluster
###############################################################################

figs1f_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::group_by(RNA_snn_res.0.2) %>%
  dplyr::summarise(
    n_studies = dplyr::n_distinct(study),
    .groups = "drop"
  ) %>%
  dplyr::rename(
    cluster = RNA_snn_res.0.2
  )

readr::write_tsv(
  figs1f_source,
  file.path(
    source_dir,
    "FigS1F_study_counts.tsv"
  )
)


###############################################################################
# Fig S1 stacked composition panels
#
# Study / Platform / Condition composition by cluster.
###############################################################################

study_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::count(
    RNA_snn_res.0.2,
    study,
    name = "cell_count"
  ) %>%
  dplyr::group_by(RNA_snn_res.0.2) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count),
    variable = "Study",
    category = as.character(study)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = RNA_snn_res.0.2,
    variable,
    category,
    cell_count,
    proportion
  )

platform_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::count(
    RNA_snn_res.0.2,
    Platform,
    name = "cell_count"
  ) %>%
  dplyr::group_by(RNA_snn_res.0.2) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count),
    variable = "Platform",
    category = as.character(Platform)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = RNA_snn_res.0.2,
    variable,
    category,
    cell_count,
    proportion
  )

condition_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::count(
    RNA_snn_res.0.2,
    Condition,
    name = "cell_count"
  ) %>%
  dplyr::group_by(RNA_snn_res.0.2) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count),
    variable = "Condition",
    category = as.character(Condition)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster = RNA_snn_res.0.2,
    variable,
    category,
    cell_count,
    proportion
  )

figs1_composition_source <- dplyr::bind_rows(
  study_source,
  platform_source,
  condition_source
)

readr::write_tsv(
  figs1_composition_source,
  file.path(
    source_dir,
    "FigS1_stacked_composition.tsv"
  )
)


###############################################################################
# sessionInfo()
###############################################################################

