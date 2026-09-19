###############################################################################
# Project: Pancreas scRNA atlas — figures S2
# File: Figures_S2.R
# Author: Elhadi Iich
# Date: 2026-09-07
#
# Description:
#   Generates supplemental figure panels S2A to S2F.
#   Original figure-generation code with source-data export added for the 
#   revised manuscript.
#
# Inputs:
#   - DEGs_allmarkers_res_0.2_ordered_pcs.tsv
#   - scRNA_integrated_harmony_classified.rds
#
# Dependencies:
#   Seurat, patchwork, tidyverse, doParallel, ggsci, scCustomize, ggpubr,
#   scRNAtoolVis, SCpubr
###############################################################################

# Packages
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(doParallel)
  library(ggsci)
  library(scCustomize)
  library(ggpubr)
  library(scRNAtoolVis)
})

options(bitmapType = 'cairo')

# Palettes
palette <- c("#E4E1E3FF", "#5A5156FF" , "#F6222EFF" ,
             "#FE00FAFF", "#16FF32FF" ,"#3283FEFF",
             "#FEAF16FF", "#B00068FF", "#1CFFCEFF" ,
             "#90AD1CFF", "#2ED9FFFF" ,"#DEA0FDFF",
             "#AA0DFEFF", "#F8A19FFF", "#325A9BFF" ,
             "#C4451CFF", "#1C8356FF" ,"#85660DFF",
             "#B10DA1FF", "#FBE426FF", "#1CBE4FFF")

palette2 <- c("#5A5156FF","#E4E1E3FF",
              "#F6222EFF","#FE00FAFF","#16FF32FF",
              "#3283FEFF","#FEAF16FF","#B00068FF",
              "#1CFFCEFF","#90AD1CFF","#2ED9FFFF")

###############################################################################
# Figure S2A
###############################################################################
cluster_markers <-
  read_tsv('DEGs_allmarkers_res_0.2_ordered_pcs.tsv')

cluster_markers <- cluster_markers %>% Add_Pct_Diff() 
str(cluster_markers)
# cluster_markers[cluster_markers$cluster == 10,] %>% View()

# pdf('S2A.pdf', width =10)
jjVolcano(diffData = cluster_markers,
          tile.col = palette2,
          fontface = 'bold',
          legend.position = c(0.9,0.9),
          cluster.order = seq(12),
          log2FC.cutoff = 0.25) + 
  xlab('Cluster') +
  theme(axis.title = element_text(face = 'bold'),
        axis.text =  element_text(size = 14, 
                                  face = 'bold'),
        legend.text = element_text(size = 14, 
                                   face = 'bold'))
dev.off()

###################################################################################
# Figure S2C
###################################################################################
merged_seurat <-
  readRDS('scRNA_integrated_harmony_classified.rds')

merged_seurat@meta.data$scpred_no_rejection <-
  gsub(merged_seurat@meta.data$scpred_no_rejection,
       pattern = ' Cell', replacement = '')

merged_seurat@meta.data %>% as.data.frame() %>%
  group_by(RNA_snn_res.0.2,scpred_no_rejection)%>% 
  dplyr::count() %>% 
  ggplot(., aes(x = factor(RNA_snn_res.0.2), y = n, fill = scpred_no_rejection)) + 
  geom_bar(stat="identity",position="fill", color = 'black')+ 
  theme_pubr(border = F,legend = 'right')+
  ylab('Proportion')+  xlab('Cluster')+
  guides(fill=guide_legend(title="Cell type prediction"))+
  scale_fill_manual(values=c("Alpha" = "#5A5156FF",
                             "Ductal" = "#E4E1E3FF",
                             "Acinar"= "#F6222EFF",
                             "Beta" = "#FE00FAFF",
                             "Pancreatic Stellate" = "#16FF32FF",
                             "Endothelial" = "#3283FEFF",
                             "Delta" = "#FEAF16FF",
                             "Pancreatic Polypeptide" = "#B00068FF",
                             "MHC Class II" ="#1CFFCEFF",
                             'Epsilon' = "#C4451CFF",
                             'Mast' = "#1C8356FF"))+
  theme(axis.text = element_text(face = 'bold', size = 14),
        axis.title = element_text(face = 'bold', size = 16),
        legend.title = element_text(face = 'bold', size = 16),
        legend.text = element_text(face = 'bold', size = 14))-> S2C

# pdf('S2C.pdf',height = 4)
S2C
dev.off()

###############################################################################
# Figure S2D
###############################################################################
pdf('S2D.pdf',height = 4)
VlnPlot_scCustom(merged_seurat, 
                 features = 'GHRL', 
                 group.by = 'RNA_snn_res.0.2') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

##############################################################
# Figure S2E
##############################################################
pdf('S2E.pdf')
DotPlot_scCustom(merged_seurat, 
                 colors_use = c('blue','red'),
                 group.by = 'RNA_snn_res.0.2',
                 features = c('GCG','INS',
                              'SST','PPY','GHRL',
                              'HSPA1A','HSPA1B','HSPH1'),
                 x_lab_rotate = 45)


SCpubr::do_DotPlot(sample = merged_seurat,
                   use_viridis = TRUE,
                   viridis.direction = 1, 
                   viridis.palette = "F",
                   group.by = 'RNA_snn_res.0.2', 
                   features = c('GCG','INS',
                                'SST','PPY','GHRL',
                                'HSPA1A','HSPA1B','HSPH1'))
dev.off()

###############################################################################
# Figure S2F
###############################################################################
pdf('S2F.pdf')
DotPlot_scCustom(seurat_object = merged_seurat, 
                 colors_use = c('blue','red'),
                 group.by = 'RNA_snn_res.0.2',
                 features = c('UBE2C','TOP2A','CDK1',
                              'BIRC5','PBK','CDKN3',
                              'MKI67','CDC20','CCNB2',
                              'CDCA3'),
                 x_lab_rotate = 45)

SCpubr::do_DotPlot(sample = merged_seurat,
                   use_viridis = TRUE,
                   viridis.direction = 1, 
                   viridis.palette = "F",
                   group.by = 'RNA_snn_res.0.2', 
                   features = c('UBE2C','TOP2A','CDK1',
                                'BIRC5','PBK','CDKN3',
                                'MKI67','CDC20','CCNB2',
                                'CDCA3'))

dev.off()


###############################################################################
# SOURCE DATA EXPORT — Figure S2
###############################################################################

source_dir <- "../SourceData/FigureS2"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Figure S2A — cluster marker volcano plot
#
# Export the complete marker table used by jjVolcano, including the
# percentage-difference columns added by Add_Pct_Diff().
###############################################################################

readr::write_tsv(
  cluster_markers,
  file.path(
    source_dir,
    "FigS2A_volcano_data.tsv.gz"
  )
)


###############################################################################
# Figure S2C — scPred cell-type composition by cluster
###############################################################################

figs2c_source <- merged_seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::count(
    RNA_snn_res.0.2,
    scpred_no_rejection,
    name = "cell_count"
  ) %>%
  dplyr::group_by(RNA_snn_res.0.2) %>%
  dplyr::mutate(
    proportion = cell_count / sum(cell_count)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    cluster = RNA_snn_res.0.2,
    predicted_cell_type = scpred_no_rejection
  )

readr::write_tsv(
  figs2c_source,
  file.path(
    source_dir,
    "FigS2C_scpred_composition.tsv"
  )
)


###############################################################################
# Figures S2D, S2E and S2F
#
# Cell-level expression values underlying:
#   S2D — GHRL violin plot
#   S2E — endocrine/stress marker dot plot
#   S2F — cycling marker dot plot
###############################################################################

s2_genes <- unique(c(
  # S2D
  "GHRL",
  
  # S2E
  "GCG", "INS", "SST", "PPY", "GHRL",
  "HSPA1A", "HSPA1B", "HSPH1",
  
  # S2F
  "UBE2C", "TOP2A", "CDK1", "BIRC5", "PBK",
  "CDKN3", "MKI67", "CDC20", "CCNB2", "CDCA3"
))

figs2_expression <- Seurat::FetchData(
  merged_seurat,
  vars = c("RNA_snn_res.0.2", s2_genes)
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::rename(
    cluster = RNA_snn_res.0.2
  )

readr::write_tsv(
  figs2_expression,
  file.path(
    source_dir,
    "FigS2D_E_F_expression.tsv.gz"
  )
)


###############################################################################
# sessionInfo()
###############################################################################

