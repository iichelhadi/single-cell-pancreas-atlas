################################################################################
# Project: scRNA beta-cell atlas - revised Lineage 2 module-trajectory analysis
# File: Fig4_IJK_tradeSeq_modules_clean.R
# Author: Elhadi Iich
# Date: 2026-09-08
#
# Purpose:
#   Generate revised Figs 4I–4K using tradeSeq-fitted expression trajectories.
#
# Rationale:
#   This workflow replaces the previous custom module-analysis block that:
#     1. selected pre/post cells by positional indices rather than matched cell IDs;
#     2. smoothed expression over cell order rather than directly using fitted
#        pseudotime trajectories; and
#     3. selected variable genes after row-wise z-scoring.
#
#   The revised workflow:
#     1. selects genes associated with pseudotime using tradeSeq associationTest;
#     2. extracts fitted GAM trajectories using tradeSeq predictSmooth;
#     3. maps the Slingshot-derived bifurcation point onto fitted pseudotime grids;
#     4. separates shared pre-branch, Lineage2_1 post-branch, and Lineage2_2
#        post-branch trajectories;
#     5. clusters standardized fitted trajectories separately for each segment;
#     6. selects K by average silhouette width;
#     7. uses K = 2 for each segment, as supported by the silhouette analysis;
#     8. performs GO enrichment relative to the pseudotime-associated gene universe.
#
# Inputs:
#   - sce_lin2.rds
#   - slingshot_sim_lin2.rds
#
# Outputs:
#   - Fig 4I–K
#   - pseudotime association results
#   - K-selection results
#   - module assignments
#   - module trajectory source data
#   - GO enrichment tables
#
# Manuscript lineage naming:
#   tradeSeq / Slingshot "Lineage1" = manuscript Lineage2_1 (toward C2)
#   tradeSeq / Slingshot "Lineage2" = manuscript Lineage2_2 (toward C5)
################################################################################

suppressPackageStartupMessages({
  library(tradeSeq)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(cluster)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

set.seed(123)

################################################################################
### 1. Paths and output directory
################################################################################

sce_path <- "sce_lin2.rds"
sim_path <- "slingshot_sim_lin2.rds"

outdir <- "Fig4_IJK_tradeSeq_outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

################################################################################
### 2. Load fitted tradeSeq and Slingshot objects
################################################################################

sce <- readRDS(sce_path)
sim <- readRDS(sim_path)

cat(
  "tradeSeq genes:", nrow(sce), "\n",
  "tradeSeq cells:", ncol(sce), "\n",
  "tradeSeq knots:", tradeSeq::nknots(sce), "\n"
)

################################################################################
### 3. Recalculate the validated Lineage 2 bifurcation
################################################################################

pt_sim <- slingshot::slingPseudotime(sim)
wt_sim <- slingshot::slingCurveWeights(sim)


pt1_all <- pt_sim[, "Lineage1"]
pt2_all <- pt_sim[, "Lineage2"]
w1 <- wt_sim[, "Lineage1"]
w2 <- wt_sim[, "Lineage2"]

names(pt1_all) <- rownames(pt_sim)
names(pt2_all) <- rownames(pt_sim)
names(w1) <- rownames(wt_sim)
names(w2) <- rownames(wt_sim)

shared_idx <- which(
  abs(w1 - w2) < 0.2 &
    !is.na(pt1_all) &
    !is.na(pt2_all)
)


bifurcation_pt <- max(
  pmin(
    pt1_all[shared_idx],
    pt2_all[shared_idx]
  ),
  na.rm = TRUE
)

cat("Estimated bifurcation pseudotime:", bifurcation_pt, "\n")

################################################################################
### 4. Identify genes associated with pseudotime
################################################################################

assoc_res <- tradeSeq::associationTest(
  sce,
  global = TRUE,
  lineages = FALSE,
  l2fc = 1
)

assoc_df <- as.data.frame(assoc_res) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    p_adj = p.adjust(pvalue, method = "BH")
  )

genes_use <- assoc_df %>%
  dplyr::filter(
    !is.na(p_adj),
    p_adj < 0.05
  ) %>%
  dplyr::pull(gene)

cat("Pseudotime-associated genes:", length(genes_use), "\n")


readr::write_tsv(
  assoc_df,
  file.path(outdir, "Fig4I-K_tradeSeq_associationTest.tsv")
)

################################################################################
### 5. Extract fitted tradeSeq trajectories
################################################################################

nPoints <- 100

yhat <- tradeSeq::predictSmooth(
  sce,
  gene = genes_use,
  nPoints = nPoints,
  tidy = FALSE
)


fit_l1 <- yhat[, seq_len(nPoints), drop = FALSE]
fit_l2 <- yhat[, nPoints + seq_len(nPoints), drop = FALSE]

################################################################################
### 6. Map fitted trajectory grids onto Slingshot pseudotime
################################################################################

pt1_range <- range(pt1_all, na.rm = TRUE)
pt2_range <- range(pt2_all, na.rm = TRUE)

grid_l1 <- seq(pt1_range[1], pt1_range[2], length.out = nPoints)
grid_l2 <- seq(pt2_range[1], pt2_range[2], length.out = nPoints)

pre_idx_l1 <- which(grid_l1 <= bifurcation_pt)
pre_idx_l2 <- which(grid_l2 <= bifurcation_pt)

post_idx_l1 <- which(grid_l1 > bifurcation_pt)
post_idx_l2 <- which(grid_l2 > bifurcation_pt)

cat(
  "Pre points L2_1:", length(pre_idx_l1), "\n",
  "Pre points L2_2:", length(pre_idx_l2), "\n",
  "Post points L2_1:", length(post_idx_l1), "\n",
  "Post points L2_2:", length(post_idx_l2), "\n"
)

################################################################################
### 7. Construct shared pre-branch and lineage-specific post-branch matrices
################################################################################

common_pre_grid <- seq(
  min(c(
    grid_l1[pre_idx_l1],
    grid_l2[pre_idx_l2]
  )),
  bifurcation_pt,
  length.out = 50
)

interp_matrix <- function(mat, x_old, x_new) {
  
  out <- t(
    apply(
      mat,
      1,
      function(y) {
        stats::approx(
          x = x_old,
          y = y,
          xout = x_new,
          rule = 2
        )$y
      }
    )
  )
  
  rownames(out) <- rownames(mat)
  out
}

fit_pre_l1 <- interp_matrix(
  fit_l1[, pre_idx_l1, drop = FALSE],
  grid_l1[pre_idx_l1],
  common_pre_grid
)

fit_pre_l2 <- interp_matrix(
  fit_l2[, pre_idx_l2, drop = FALSE],
  grid_l2[pre_idx_l2],
  common_pre_grid
)

fit_pre <- (fit_pre_l1 + fit_pre_l2) / 2

fit_post_l1 <- fit_l1[, post_idx_l1, drop = FALSE]
fit_post_l2 <- fit_l2[, post_idx_l2, drop = FALSE]

################################################################################
### 8. Row-standardize fitted trajectories
################################################################################

scale_rows <- function(mat) {
  
  z <- t(scale(t(mat)))
  
  keep <- apply(
    z,
    1,
    function(x) all(is.finite(x))
  )
  
  z[keep, , drop = FALSE]
}

z_pre <- scale_rows(fit_pre)
z_l1 <- scale_rows(fit_post_l1)
z_l2 <- scale_rows(fit_post_l2)

cat(
  "Pre matrix:", nrow(z_pre), "genes x", ncol(z_pre), "points\n",
  "L2_1 post matrix:", nrow(z_l1), "genes x", ncol(z_l1), "points\n",
  "L2_2 post matrix:", nrow(z_l2), "genes x", ncol(z_l2), "points\n"
)

################################################################################
### 9. Select module number using average silhouette width
################################################################################

choose_k <- function(mat, k_range = 2:6) {
  
  pca <- stats::prcomp(
    mat,
    center = FALSE,
    scale. = FALSE
  )
  
  npc <- min(10, ncol(pca$x))
  x <- pca$x[, seq_len(npc), drop = FALSE]
  d <- stats::dist(x)
  
  out <- lapply(
    k_range,
    function(k) {
      
      set.seed(123)
      
      km <- stats::kmeans(
        x,
        centers = k,
        nstart = 50,
        iter.max = 200,
        algorithm = "Lloyd"
      )
      
      sil <- cluster::silhouette(km$cluster, d)
      
      data.frame(
        k = k,
        average_silhouette = mean(sil[, "sil_width"])
      )
    }
  )
  
  dplyr::bind_rows(out)
}

k_pre <- choose_k(z_pre)
k_l1 <- choose_k(z_l1)
k_l2 <- choose_k(z_l2)

cat("\nPre-branch K selection:\n")
print(k_pre)

cat("\nLineage2_1 K selection:\n")
print(k_l1)

cat("\nLineage2_2 K selection:\n")
print(k_l2)

best_k_pre <- k_pre$k[which.max(k_pre$average_silhouette)]
best_k_l1 <- k_l1$k[which.max(k_l1$average_silhouette)]
best_k_l2 <- k_l2$k[which.max(k_l2$average_silhouette)]

if (!all(c(best_k_pre, best_k_l1, best_k_l2) == 2)) {
  warning(
    "K = 2 is not the highest-silhouette solution for all segments. ",
    "Review K-selection before using the fixed K = 2 Fig."
  )
}

k_results <- dplyr::bind_rows(
  dplyr::mutate(k_pre, segment = "Pre-branch"),
  dplyr::mutate(k_l1, segment = "Lineage2_1_post"),
  dplyr::mutate(k_l2, segment = "Lineage2_2_post")
) %>%
  dplyr::select(segment, dplyr::everything())

readr::write_tsv(
  k_results,
  file.path(outdir, "Fig4I-K_K_silhouette.tsv")
)

################################################################################
### 10. Final clustering: K = 2
################################################################################

cluster_segment <- function(mat, k = 2) {
  
  pca <- stats::prcomp(
    mat,
    center = FALSE,
    scale. = FALSE
  )
  
  npc <- min(10, ncol(pca$x))
  x <- pca$x[, seq_len(npc), drop = FALSE]
  
  set.seed(123)
  
  km <- stats::kmeans(
    x,
    centers = k,
    nstart = 100,
    iter.max = 200,
    algorithm = "Lloyd"
  )
  
  clusters <- km$cluster
  names(clusters) <- rownames(mat)
  
  list(
    cluster = clusters,
    pca = pca,
    kmeans = km
  )
}

cl_pre_raw <- cluster_segment(z_pre, k = 2)
cl_l1_raw <- cluster_segment(z_l1, k = 2)
cl_l2_raw <- cluster_segment(z_l2, k = 2)

################################################################################
### 11. Canonicalize module numbering
################################################################################

canonicalize_k2 <- function(mat, clusters, segment = c("pre", "post")) {
  
  segment <- match.arg(segment)
  cluster_ids <- sort(unique(clusters))
  
  deltas <- sapply(
    cluster_ids,
    function(cl) {
      
      genes <- names(clusters)[clusters == cl]
      
      mean(
        mat[genes, ncol(mat), drop = TRUE],
        na.rm = TRUE
      ) -
        mean(
          mat[genes, 1, drop = TRUE],
          na.rm = TRUE
        )
    }
  )
  
  increasing_cluster <- cluster_ids[which.max(deltas)]
  other_cluster <- setdiff(cluster_ids, increasing_cluster)
  
  canonical <- rep(NA_integer_, length(clusters))
  names(canonical) <- names(clusters)
  
  if (segment == "pre") {
    canonical[clusters == increasing_cluster] <- 2L
    canonical[clusters == other_cluster] <- 1L
  } else {
    canonical[clusters == increasing_cluster] <- 1L
    canonical[clusters == other_cluster] <- 2L
  }
  
  mapping <- tibble::tibble(
    raw_cluster = cluster_ids,
    delta_end_minus_start = as.numeric(deltas),
    canonical_module = if (segment == "pre") {
      ifelse(cluster_ids == increasing_cluster, 2L, 1L)
    } else {
      ifelse(cluster_ids == increasing_cluster, 1L, 2L)
    }
  )
  
  list(
    cluster = canonical,
    mapping = mapping
  )
}

canon_pre <- canonicalize_k2(z_pre, cl_pre_raw$cluster, segment = "pre")
canon_l1 <- canonicalize_k2(z_l1, cl_l1_raw$cluster, segment = "post")
canon_l2 <- canonicalize_k2(z_l2, cl_l2_raw$cluster, segment = "post")

cl_pre <- canon_pre$cluster
cl_l1 <- canon_l1$cluster
cl_l2 <- canon_l2$cluster

cat("\nCanonical module sizes:\n")

cat("Pre-branch:\n")
print(table(cl_pre))

cat("Lineage2_1:\n")
print(table(cl_l1))

cat("Lineage2_2:\n")
print(table(cl_l2))

canonical_mapping <- dplyr::bind_rows(
  dplyr::mutate(canon_pre$mapping, segment = "Pre-branch"),
  dplyr::mutate(canon_l1$mapping, segment = "Lineage2_1_post"),
  dplyr::mutate(canon_l2$mapping, segment = "Lineage2_2_post")
) %>%
  dplyr::select(segment, dplyr::everything())

readr::write_tsv(
  canonical_mapping,
  file.path(outdir, "Fig4I-K_cluster_label_mapping.tsv")
)

################################################################################
### 12. Final gene-module assignments
################################################################################

modules_pre <- tibble::tibble(
  gene = names(cl_pre),
  module = paste0("Pre_M", cl_pre)
)

modules_l1 <- tibble::tibble(
  gene = names(cl_l1),
  module = paste0("L2_1_M", cl_l1)
)

modules_l2 <- tibble::tibble(
  gene = names(cl_l2),
  module = paste0("L2_2_M", cl_l2)
)

readr::write_tsv(
  modules_pre,
  file.path(outdir, "Fig4I_tradeSeq_prebranch_modules.tsv")
)

readr::write_tsv(
  modules_l1,
  file.path(outdir, "Fig4J_tradeSeq_Lineage2_1_modules.tsv")
)

readr::write_tsv(
  modules_l2,
  file.path(outdir, "Fig4K_tradeSeq_Lineage2_2_modules.tsv")
)

################################################################################
### 13. Build module-trajectory source data
################################################################################

make_module_plot_data <- function(
    mat,
    clusters,
    pseudotime,
    prefix
) {
  
  mods <- sort(unique(clusters))
  
  lapply(
    mods,
    function(m) {
      
      genes_m <- names(clusters)[clusters == m]
      
      tibble::tibble(
        pseudotime = pseudotime,
        expression = colMeans(
          mat[genes_m, , drop = FALSE],
          na.rm = TRUE
        ),
        module = paste0(prefix, "_M", m),
        n_genes = length(genes_m)
      )
    }
  ) %>%
    dplyr::bind_rows()
}

plot_pre_df <- make_module_plot_data(
  z_pre,
  cl_pre,
  common_pre_grid,
  "Pre"
)

plot_l1_df <- make_module_plot_data(
  z_l1,
  cl_l1,
  grid_l1[post_idx_l1],
  "L2_1"
)

plot_l2_df <- make_module_plot_data(
  z_l2,
  cl_l2,
  grid_l2[post_idx_l2],
  "L2_2"
)

readr::write_tsv(
  plot_pre_df,
  file.path(outdir, "Fig4I_module_trajectory_source_data.tsv")
)

readr::write_tsv(
  plot_l1_df,
  file.path(outdir, "Fig4J_module_trajectory_source_data.tsv")
)

readr::write_tsv(
  plot_l2_df,
  file.path(outdir, "Fig4K_module_trajectory_source_data.tsv")
)

################################################################################
### 14. GO enrichment
################################################################################

run_module_go <- function(
    module_table,
    universe_genes,
    ont = "BP"
) {
  
  module_list <- split(
    module_table$gene,
    module_table$module
  )
  
  lapply(
    module_list,
    function(genes) {
      
      clusterProfiler::enrichGO(
        gene = genes,
        universe = universe_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500
      )
    }
  )
}

go_pre <- run_module_go(modules_pre, genes_use)
go_l1 <- run_module_go(modules_l1, genes_use)
go_l2 <- run_module_go(modules_l2, genes_use)

################################################################################
### 15. Print and export GO results
################################################################################

print_module_go <- function(
    go_list,
    label,
    n = 20
) {
  
  cat("\n###", label, "###\n")
  
  for (nm in names(go_list)) {
    
    cat("\n", nm, "\n")
    ego <- go_list[[nm]]
    
    if (is.null(ego) || nrow(ego@result) == 0) {
      cat("No GO terms returned\n")
      next
    }
    
    tab <- ego@result %>%
      as.data.frame() %>%
      dplyr::arrange(p.adjust, pvalue) %>%
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust
      ) %>%
      head(n)
    
    print(tab, row.names = FALSE)
  }
}

print_module_go(go_pre, "Pre-branch")
print_module_go(go_l1, "Lineage2_1")
print_module_go(go_l2, "Lineage2_2")

write_go_tables <- function(go_list, output_dir) {
  
  for (nm in names(go_list)) {
    
    ego <- go_list[[nm]]
    
    if (is.null(ego)) {
      next
    }
    
    readr::write_tsv(
      as.data.frame(ego),
      file.path(
        output_dir,
        paste0(nm, "_GO_BP.tsv")
      )
    )
  }
}

write_go_tables(go_pre, outdir)
write_go_tables(go_l1, outdir)
write_go_tables(go_l2, outdir)

################################################################################
### 16. Functional labels used in Figs 4I–4K
#
# * = descriptive functional annotation based on leading GO terms;
#     no GO term reached BH-adjusted P < 0.05 for that module.
################################################################################

module_annotations <- tibble::tribble(
  ~module,    ~label,                                       ~go_significant,
  "Pre_M1",   "Ciliary Organization +\nProteostasis",        TRUE,
  "Pre_M2",   "Chromatin +\nRNA Regulation*",                FALSE,
  "L2_1_M1",  "Inflammatory +\nDefense Response",            TRUE,
  "L2_1_M2",  "Chromatin +\nRNA Processing*",                FALSE,
  "L2_2_M1",  "Wound Healing +\nCell Migration",             TRUE,
  "L2_2_M2",  "Nucleotide Metabolism +\nVesicle Transport*", FALSE
)

readr::write_tsv(
  module_annotations,
  file.path(outdir, "Fig4I-K_module_annotations.tsv")
)

################################################################################
### 17. Plot Figs 4I–4K
################################################################################

asterisk_note <-
  "* Descriptive functional annotation based on leading GO terms; no GO term reached BH-adjusted P < 0.05."

plot_module_panel <- function(
    plot_df,
    module_ids,
    title,
    ylab_text = NULL,
    significant_color = "#E31A1C"
) {
  
  ann <- module_annotations %>%
    dplyr::filter(module %in% module_ids)
  
  label_vec <- stats::setNames(ann$label, ann$module)
  
  color_vec <- stats::setNames(
    ifelse(
      ann$go_significant,
      significant_color,
      "#999999"
    ),
    ann$module
  )
  
  ggplot(
    plot_df,
    aes(
      x = pseudotime,
      y = expression,
      color = module
    )
  ) +
    geom_line(linewidth = 1.3) +
    geom_vline(
      xintercept = bifurcation_pt,
      linetype = "dotted",
      linewidth = 0.6,
      color = "blue"
    ) +
    scale_color_manual(
      values = color_vec,
      breaks = module_ids,
      labels = label_vec
    ) +
    labs(
      title = title,
      x = "Pseudotime",
      y = ylab_text,
      color = "Module"
    ) +
    ggpubr::theme_pubr(legend = 'right') +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      axis.title = element_text(
        face = "bold",
        size = 14
      ),
      axis.text = element_text(
        face = "bold",
        size = 13
      ),
      legend.title = element_text(
        face = "bold",
        size = 13
      ),
      legend.text = element_text(
        face = "bold",
        size = 12
      )
    )
}

library(geomtextpath)

p_pre <- plot_module_panel(
  plot_pre_df,
  module_ids = c("Pre_M1", "Pre_M2"),
  title = "Pre-bifurcation",
  ylab_text = "Mean scaled fitted expression",
  significant_color = "black"
)+
  coord_cartesian(
    xlim = c(0, 14)) +
  geom_vline(
    xintercept = bifurcation_pt,
    linetype = "dotted",
    color = "blue",
    linewidth = 1.5
  ) +
  geom_textvline(
    label = "Pre bifurcation",
    xintercept = bifurcation_pt,
    vjust = 1.2,
    size = 5,
    fontface = "bold"
  )

p_l1 <- plot_module_panel(
  plot_l1_df,
  module_ids = c("L2_1_M1", "L2_1_M2"),
  title = "Lineage2_1",
  ylab_text = NULL,
  significant_color = "#F62424"
)+
  geom_vline(
    xintercept = bifurcation_pt,
    linetype = "dotted",
    color = "blue",
    linewidth = 1.5
  ) +
  coord_cartesian(
    xlim = c(13.2, 18.3))+
  geom_textvline(
    label = "Post bifurcation",
    xintercept = bifurcation_pt,
    vjust = -0.3,
    size = 5,
    fontface = "bold"
  )

p_l2 <- plot_module_panel(
  plot_l2_df,
  module_ids = c("L2_2_M1", "L2_2_M2"),
  title = "Lineage2_2",
  ylab_text = NULL,
  significant_color = "#ffb300"
)+
  geom_vline(
    xintercept = bifurcation_pt,
    linetype = "dotted",
    color = "blue",
    linewidth = 1.5
  ) +
  coord_cartesian(
    xlim = c(13.2, 18.3))+
  geom_textvline(
    label = "Post bifurcation",
    xintercept = bifurcation_pt,
    vjust = -0.3,
    size = 5,
    fontface = "bold"
  )

Fig4_IJK <- (
  p_pre |
    p_l1 |
    p_l2
) +
  patchwork::plot_annotation(
    caption = asterisk_note
  )

print(Fig4_IJK)


################################################################################
### 18. Export analysis summary and session information
################################################################################

analysis_summary <- tibble::tibble(
  metric = c(
    "Bifurcation pseudotime",
    "Pseudotime-associated genes",
    "Pre fitted points L2_1",
    "Pre fitted points L2_2",
    "Post fitted points L2_1",
    "Post fitted points L2_2",
    "Best K pre-branch",
    "Best K Lineage2_1 post",
    "Best K Lineage2_2 post"
  ),
  value = c(
    as.character(bifurcation_pt),
    as.character(length(genes_use)),
    as.character(length(pre_idx_l1)),
    as.character(length(pre_idx_l2)),
    as.character(length(post_idx_l1)),
    as.character(length(post_idx_l2)),
    as.character(best_k_pre),
    as.character(best_k_l1),
    as.character(best_k_l2)
  )
)

readr::write_tsv(
  analysis_summary,
  file.path(outdir, "Fig4I-K_analysis_summary.tsv")
)

###############################################################################
# sessionInfo()
###############################################################################


