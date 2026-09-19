################################################################################
# Project: scRNA beta cell atlas - Figure 5H SCENIC figures
# File: Figure5H_SCENIC_figures.R
# Author: Elhadi Iich
#
# Description:
#   Generates the revised Figure 5H heatmap, supplementary donor-level
#   regulon boxplots, and SCENIC robustness table from saved analysis outputs.
################################################################################


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


fig5h_donor_data <- read_tsv("Figure5H_regulon_donor_data.tsv",
                             show_col_types = FALSE)

pairwise_results <- read_tsv("Figure5H_pairwise_mixed_model.tsv",
                             show_col_types = FALSE)

robustness_summary <- read_tsv("Figure5H_regulon_robustness_summary.tsv",
                               show_col_types = FALSE)

sig70_threshold_check <- read_tsv("Figure5H_sig70_regulons_across_thresholds.tsv",
                                  show_col_types = FALSE)


################################################################################
### Supplementary figure: donor-level regulon activity
################################################################################

plot_data <- fig5h_donor_data %>%
  dplyr::mutate(
    regulon = gsub("\\(\\+\\)", "", regulon),
    LineageGroup = factor(LineageGroup,
                          levels = c("Lineage2_1",
                                     "Shared",
                                     "Lineage2_2"))
  )

sig_pairs <- pairwise_results %>%
  dplyr::filter(p_adj < 0.05) %>%
  dplyr::mutate(
    regulon = gsub("\\(\\+\\)", "", regulon),
    group1 = sub(" - .*", "", contrast),
    group2 = sub(".* - ", "", contrast),
    p_label = dplyr::case_when(p_adj < 0.001 ~ "***",
                               p_adj < 0.01  ~ "**",
                               p_adj < 0.05  ~ "*",
                               TRUE ~ "ns")
  )

# calculate annotation heights separately for each regulon
y_positions <- plot_data %>%
  dplyr::group_by(regulon) %>%
  dplyr::summarise(
    ymax = max(mean_auc, na.rm = TRUE),
    ymin = min(mean_auc, na.rm = TRUE),
    yrange = ymax - ymin,
    .groups = "drop"
  )

sig_pairs <- sig_pairs %>%
  dplyr::left_join(y_positions,
                   by = "regulon") %>%
  dplyr::group_by(regulon) %>%
  dplyr::arrange(p_adj, .by_group = TRUE) %>%
  dplyr::mutate(y.position = ymax + (0.08 + 0.10 * (dplyr::row_number() - 1)) * yrange) %>%
  dplyr::ungroup()

p_fig5h_stats <- ggplot(plot_data,
                        aes(x = LineageGroup,
                            y = mean_auc,
                            fill = LineageGroup)) +
  geom_boxplot(outlier.shape = NA,
               width = 0.65) +
  geom_jitter(width = 0.12,
              alpha = 0.45,
              size = 0.8) +
  facet_wrap(~ regulon,
             scales = "free_y",
             ncol = 4) +
  ggpubr::stat_pvalue_manual(
    sig_pairs,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.35,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  labs(x = NULL,
       y = "Mean regulon AUC per donor") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    legend.position = "none"
  )

p_fig5h_stats


################################################################################
### Revised Figure 5H - heatmap of robust donor-level regulon activity
################################################################################

GROUP_ORDER <- c("Lineage2_1",
                 "Shared",
                 "Lineage2_2")

################################################################################
### 1. Define final robust regulons
################################################################################

robust_regulons <- robustness_summary %>%
  dplyr::filter(n_thresholds_FDR05 >= 2,
                same_peak_all,!any_singular) %>%
  dplyr::pull(regulon)

cat("\nRobust regulons used in Figure 5H:\n")
print(robust_regulons)

################################################################################
### 2. Calculate mean donor-level AUC within each lineage group
################################################################################

heatmap_df <- fig5h_donor_data %>%
  dplyr::filter(regulon %in% robust_regulons) %>%
  dplyr::mutate(LineageGroup = factor(LineageGroup,
                                      levels = GROUP_ORDER)) %>%
  dplyr::group_by(regulon,
                  LineageGroup) %>%
  dplyr::summarise(mean_auc = mean(mean_auc, na.rm = TRUE),
                   .groups = "drop") %>%
  tidyr::pivot_wider(names_from = LineageGroup,
                     values_from = mean_auc)

################################################################################
### 3. Convert to regulon x lineage matrix
################################################################################

mat_auc <- heatmap_df %>%
  tibble::column_to_rownames("regulon") %>%
  as.matrix()

mat_auc <- mat_auc[,
                   GROUP_ORDER,
                   drop = FALSE]

print(mat_auc)

################################################################################
### 4. Row-scale each regulon
#
# This shows RELATIVE activity of each regulon across the three lineage groups.
# Do not compare absolute colour intensity between different regulons.
################################################################################

mat_z <- t(scale(t(mat_auc),
                 center = TRUE,
                 scale = TRUE))

################################################################################
### 5. Order regulons by reproducible peak lineage
################################################################################

robust_order <- robustness_summary %>%
  dplyr::filter(regulon %in% robust_regulons) %>%
  dplyr::mutate(peak_group = factor(peak_groups,
                                    levels = c("Lineage2_1",
                                               "Shared",
                                               "Lineage2_2"))) %>%
  dplyr::arrange(peak_group,
                 regulon) %>%
  dplyr::pull(regulon)

mat_z <- mat_z[robust_order,
               ,
               drop = FALSE]

################################################################################
### 6. Clean labels
################################################################################

rownames(mat_z) <- gsub("\\(\\+\\)",
                        "",
                        rownames(mat_z))

colnames(mat_z) <- c("Lineage2_1",
                     "Shared",
                     "Lineage2_2")

################################################################################
### 7. Colour scale
################################################################################

max_abs <- max(abs(mat_z),
               na.rm = TRUE)

col_fun <- circlize::colorRamp2(c(-max_abs,
                                  0,
                                  max_abs),
                                c("#2166AC",
                                  "white",
                                  "#B2182B"))

################################################################################
### 8. Draw focused Figure 5H heatmap
################################################################################

ht_fig5h <- ComplexHeatmap::Heatmap(
  mat_z,
  
  name = "Relative\nregulon activity",
  
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_side = "left",
  
  row_names_gp = grid::gpar(fontsize = 11,
                            fontface = "bold"),
  
  column_names_gp = grid::gpar(fontsize = 10,
                               fontface = "bold"),
  
  column_names_rot = 0,
  column_names_centered = TRUE,
  
  rect_gp = grid::gpar(col = "white",
                       lwd = 1),
  
  border = TRUE,
  
  heatmap_legend_param = list(
    title = "Relative\nregulon activity",
    title_gp = grid::gpar(fontsize = 10,
                          fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9)
  )
)

ComplexHeatmap::draw(ht_fig5h,
                     heatmap_legend_side = "right")

################################################################################
### 9. Save PDF and PNG
################################################################################
#
# pdf(
#   "Figure5H_SCENIC_heatmap.pdf",
#   width = 5,
#   height = 6
# )
#
# ComplexHeatmap::draw(
#   ht_fig5h,
#   heatmap_legend_side = "right"
# )
#
# dev.off()
#


################################################################################
### 10. Save values underlying heatmap
################################################################################

heatmap_source_data <- as.data.frame(mat_z) %>%
  tibble::rownames_to_column("regulon")

readr::write_tsv(heatmap_source_data,
                 "Figure5H_heatmap_source_data.tsv")

cat("\nFigure 5H complete.\n")



