################################################################################
# Project: scRNA beta-cell atlas — C7 donor-level HbA1c analysis
# File: Fig2I_S3F_HbA1c_donor_analysis.R
# Author: Elhadi Iich
#
# Purpose:
#   1. Reconstruct harmonized donor IDs.
#   2. Generate Fig 2I: donor-level HbA1c distributions by beta-cell cluster.
#   3. Quantify donor-level beta-cell cluster occupancy.
#   4. Test continuous associations between HbA1c and cluster occupancy.
#   5. Generate Fig S3F: donor HbA1c versus cluster occupancy.
#   6. Assess sensitivity to HbA1c thresholds of 5.7%, 6.0%, and 6.5%.
#   7. Export source-data tables and analysis summaries.
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggbeeswarm)
  library(ggh4x)
})

set.seed(123)

################################################################################
### 1. Paths and plotting settings
################################################################################

beta_path <- "./revision/Fig2_S3/beta.integrated.rds"

ann81608_path <-
  "./revision/R1-9__R2-15/revision_out/GSE81608.csv"

ann153855_path <-
  "./revision/R1-9__R2-15/revision_out/GSE153855.csv"

outdir <- "Figure2_S3"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

desired_order <- c("3", "5", "1", "4", "6", "2", "7")

# Cluster palette used for panel labels/fills.
cluster_cols <- scales::hue_pal()(7)
names(cluster_cols) <- as.character(seq_len(7))

################################################################################
### 2. Load beta-cell object and metadata
################################################################################

beta.integrated <- readRDS(beta_path)

metadata <- beta.integrated@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::transmute(
    cell_id,
    orig.ident = as.character(orig.ident),
    sample_name = as.character(sample_name),
    Study = as.character(study),
    cluster = as.character(monocle_clusters),
    Age = as.numeric(Age),
    BMI = as.numeric(BMI),
    Hba1c = as.numeric(Hba1c),
    Gender = as.character(Gender),
    Ethnicity = as.character(Ethnicity),
    Condition = as.character(Condition)
  )


################################################################################
### 3. Load explicit donor mappings
################################################################################

ann81608 <- readr::read_csv(
  ann81608_path,
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    sample_name = as.character(Run),
    annotated_donor = as.character(donor_id)
  ) %>%
  dplyr::filter(
    !is.na(sample_name),
    sample_name != "",
    !is.na(annotated_donor),
    annotated_donor != ""
  ) %>%
  dplyr::distinct()

ann153855 <- readr::read_csv(
  ann153855_path,
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    sample_name = as.character(Run),
    annotated_biosample = as.character(BioSample)
  ) %>%
  dplyr::filter(
    !is.na(sample_name),
    sample_name != "",
    !is.na(annotated_biosample),
    annotated_biosample != ""
  ) %>%
  dplyr::distinct()

metadata <- metadata %>%
  dplyr::left_join(ann81608, by = "sample_name") %>%
  dplyr::left_join(ann153855, by = "sample_name")

################################################################################
### 4. Reconstruct harmonized donor IDs
################################################################################

fallback_studies <- c(
  "GSE98887",
  "GSE124742",
  "EMTAB5061",
  "GSE86469",
  "GSE81547"
)

signature_value <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x
}

metadata <- metadata %>%
  dplyr::mutate(
    donor_signature = paste(
      signature_value(Age),
      signature_value(BMI),
      signature_value(Hba1c),
      signature_value(Gender),
      signature_value(Ethnicity),
      signature_value(Condition),
      sep = "|"
    ),
    donor_id = dplyr::case_when(
      orig.ident == "GSE81608" ~
        paste0("GSE81608::", annotated_donor),
      orig.ident == "GSE153855" ~
        paste0("GSE153855::", annotated_biosample),
      Study %in% fallback_studies ~
        paste0(Study, "::metadata::", donor_signature),
      TRUE ~
        paste0(Study, "::", orig.ident)
    )
  )

################################################################################
### 5. Verify donor-level HbA1c consistency
################################################################################

metadata_inconsistency <- metadata %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(
    n_HbA1c = dplyr::n_distinct(Hba1c, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_HbA1c > 1)

if (nrow(metadata_inconsistency) > 0) {
  print(metadata_inconsistency)
  stop("Some reconstructed donors have inconsistent HbA1c values.")
}

################################################################################
### 6. Donor-level HbA1c distribution among donors represented in each cluster
################################################################################

hba1c_donor_cluster <- metadata %>%
  dplyr::filter(
    !is.na(Hba1c),
    !is.na(cluster)
  ) %>%
  dplyr::group_by(donor_id, cluster) %>%
  dplyr::summarise(
    Study = dplyr::first(Study),
    Hba1c = dplyr::first(Hba1c),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

cluster_hba1c_summary <- hba1c_donor_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_donors = dplyr::n(),
    donor_mean = mean(Hba1c, na.rm = TRUE),
    donor_median = median(Hba1c, na.rm = TRUE),
    donor_IQR = IQR(Hba1c, na.rm = TRUE),
    .groups = "drop"
  )

################################################################################
### 7. Build complete donor x cluster occupancy matrix
################################################################################

all_donors <- metadata %>%
  dplyr::filter(!is.na(Hba1c)) %>%
  dplyr::distinct(donor_id, Hba1c)

all_clusters <- sort(unique(metadata$cluster))

donor_cluster_occ <- metadata %>%
  dplyr::filter(
    !is.na(Hba1c),
    !is.na(cluster)
  ) %>%
  dplyr::count(donor_id, cluster, name = "n_cells") %>%
  tidyr::complete(
    donor_id = all_donors$donor_id,
    cluster = all_clusters,
    fill = list(n_cells = 0)
  ) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(
    total_beta_cells = sum(n_cells),
    occupancy = n_cells / total_beta_cells
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(all_donors, by = "donor_id")

################################################################################
### 8. Continuous association: HbA1c vs donor-level cluster occupancy
################################################################################

occupancy_cor <- donor_cluster_occ %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_donors = dplyr::n(),
    rho = suppressWarnings(
      cor(Hba1c, occupancy, method = "spearman")
    ),
    p_value = suppressWarnings(
      cor.test(
        Hba1c,
        occupancy,
        method = "spearman",
        exact = FALSE
      )$p.value
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_adj)

print(occupancy_cor)

################################################################################
### 9. Sensitivity to alternative HbA1c thresholds
################################################################################

thresholds <- c(5.7, 6.0, 6.5)

threshold_results <- lapply(
  thresholds,
  function(thr) {
    donor_cluster_occ %>%
      dplyr::mutate(
        glycaemic_group = ifelse(Hba1c >= thr, "Above", "Below")
      ) %>%
      dplyr::group_by(cluster) %>%
      dplyr::group_modify(
        ~{
          wt <- wilcox.test(
            occupancy ~ glycaemic_group,
            data = .x,
            exact = FALSE
          )
          
          tibble::tibble(
            threshold = thr,
            n_below = sum(.x$Hba1c < thr),
            n_above = sum(.x$Hba1c >= thr),
            median_below = median(
              .x$occupancy[.x$Hba1c < thr],
              na.rm = TRUE
            ),
            median_above = median(
              .x$occupancy[.x$Hba1c >= thr],
              na.rm = TRUE
            ),
            delta_occupancy = median_above - median_below,
            p_value = wt$p.value
          )
        }
      ) %>%
      dplyr::ungroup()
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(threshold) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      delta_occupancy > 0 ~ "Higher above threshold",
      delta_occupancy < 0 ~ "Lower above threshold",
      TRUE ~ "No change"
    ),
    significant = p_adj < 0.05
  )

print(threshold_results, n = Inf)

################################################################################
### 10. Fig 2I — donor-level HbA1c distributions by beta-cell cluster
################################################################################

n_labels <- hba1c_donor_cluster %>%
  dplyr::count(cluster, name = "n_donors") %>%
  dplyr::mutate(
    cluster = factor(cluster, levels = desired_order),
    label = paste0("n=", n_donors)
  )

Fig2I <- ggplot(
  hba1c_donor_cluster,
  aes(
    x = factor(cluster, levels = desired_order),
    y = Hba1c
  )
) +
  geom_boxplot(
    aes(fill = factor(cluster, levels = desired_order)),
    outlier.shape = NA,
    width = 0.55,
    alpha = 1
  ) +
  ggbeeswarm::geom_quasirandom(
    width = 0.18,
    size = 1.15,
    alpha = 1
  ) +
  geom_hline(
    yintercept = 5.7,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 6.5,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_text(
    data = n_labels,
    aes(x = cluster, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.3,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = cluster_cols) +
  labs(
    x = "Cluster",
    y = "Donor HbA1c (%)"
  ) +
  coord_cartesian(ylim = c(4, 12), clip = "off") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14)
  )

print(Fig2I)

################################################################################
### 11. Fig S3F — donor HbA1c versus beta-cell cluster occupancy
################################################################################

occ_plot <- donor_cluster_occ %>%
  dplyr::left_join(
    occupancy_cor %>%
      dplyr::select(cluster, rho, p_adj),
    by = "cluster"
  ) %>%
  dplyr::mutate(
    cluster = factor(cluster, levels = desired_order)
  )

occ_labels <- occupancy_cor %>%
  dplyr::mutate(
    cluster = factor(cluster, levels = desired_order),
    lab = sprintf(
      "rho = %.2f\nBH p = %.3f",
      rho,
      p_adj
    )
  )

strip_cols <- cluster_cols[desired_order]

FigS3F <- ggplot(
  occ_plot,
  aes(x = Hba1c, y = occupancy)
) +
  geom_vline(
    xintercept = c(5.7, 6.0, 6.5),
    linetype = "dashed",
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    size = 1.2,
    alpha = 0.6
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 3),
    se = FALSE,
    linewidth = 0.6,
    colour = "#0072B2"
  ) +
  geom_text(
    data = occ_labels,
    aes(x = Inf, y = Inf, label = lab),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.2,
    fontface = "bold"
  ) +
  ggh4x::facet_wrap2(
    ~ cluster,
    nrow = 2,
    strip = ggh4x::strip_themed(
      background_x = lapply(
        strip_cols,
        function(x) element_rect(fill = x, colour = "black")
      ),
      text_x = lapply(
        strip_cols,
        function(x) element_text(face = "bold", size = 14)
      )
    )
  ) +
  labs(
    x = "Donor HbA1c (%)",
    y = "Proportion of donor\nbeta cells in cluster"
  ) +
  ggpubr::theme_pubr() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14)
  )

print(FigS3F)

################################################################################
### 12. Export source-data tables
################################################################################

readr::write_tsv(
  hba1c_donor_cluster,
  file.path(
    outdir,
    "Fig2I_HbA1c_donor_distribution_source_data.tsv"
  )
)

readr::write_tsv(
  cluster_hba1c_summary,
  file.path(
    outdir,
    "Fig2I_HbA1c_cluster_summary.tsv"
  )
)

readr::write_tsv(
  donor_cluster_occ,
  file.path(
    outdir,
    "FigS3F_donor_cluster_occupancy_source_data.tsv"
  )
)

readr::write_tsv(
  occupancy_cor,
  file.path(
    outdir,
    "FigS3F_HbA1c_cluster_occupancy_associations.tsv"
  )
)

readr::write_tsv(
  threshold_results,
  file.path(
    outdir,
    "FigS3F_HbA1c_threshold_occupancy_sensitivity.tsv"
  )
)

################################################################################
### 13. Save Figs
################################################################################

ggsave(
  file.path(outdir, "Fig2I_HbA1c_donor_distribution.pdf"),
  Fig2I,
  width = 6.5,
  height = 4.5
)

ggsave(
  file.path(outdir, "FigS3F_HbA1c_cluster_occupancy.pdf"),
  FigS3F,
  width = 8,
  height = 6
)

################################################################################
# sessionInfo()
################################################################################
