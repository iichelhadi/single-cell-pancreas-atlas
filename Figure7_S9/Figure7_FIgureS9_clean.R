###############################################################################
# Project: Pancreas scRNA atlas — Static MI Network and Hypergraph Modules
# File: Figure7_FIgureS9_clean.R
# Author: Elhadi Iich
# Date: 2025-08-29
#
# Description:
#   Constructs mutual information (MI) networks for L2_C5 and L3_C3 clusters,
#   applies percolation and community detection, post-processes with 
#   neighborhood overlap, and builds hypergraphs. Outputs include:
#     - Consensus MI networks (Leiden + label propagation modules)
#     - Hypergraph incidence matrices
#     - Cluster/community summaries (size, density, hub genes)
#     - GO enrichment, bar plots, dot plots
#     - Network and chord diagrams for SCENIC–MI TF–target integration
#     - Module-score trajectories along pseudotime
#
# Inputs:
#   - lin2.beta.reintegrated.rds
#   - lin3.beta.reintegrated.rds
#   - SCENIC regulonTargetsInfo RDS files (L2_C5, L3_C3)
#
# Dependencies:
#   Seurat, infotheo, igraph, parallel, ggplot2, tidyverse, ggraph, ggbreak,
#   clusterProfiler, org.Hs.eg.db, patchwork, flextable, VennDiagram, hypergraph,
#   ComplexHeatmap, circlize, ggpubr
#
# Notes:
# MI-hypergraph analysis packaged into the R package miEdgeR
# https://github.com/iichelhadi/miEdgeR
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(infotheo)
  library(igraph)
  library(parallel)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(patchwork)
  library(tidyverse)
  library(ggraph)
  library(ggbreak)
  library(dplyr)
  library(flextable)
  library(VennDiagram)
  library(hypergraph)
  library(circlize)
  library(ggpubr)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
num_cores <- detectCores()

# Define lineage end clusters
lineage_clusters <- list(
  L2_C5 = list(obj = readRDS("lin2.beta.reintegrated.rds"), cluster = "5"),
  L3_C3 = list(obj = readRDS("lin3.beta.reintegrated.rds"), cluster = "3")
)

# Preload Seurat objects
seurat_objects <- lapply(lineage_clusters, function(x) x$obj)
names(seurat_objects) <- names(lineage_clusters)

# rm(lineage_clusters)

# Utility Functions
filter_housekeeping <- function(gene_list) {
  gene_list[!grepl("^(RPL|RPS|GAPDH|ACTB)", gene_list)]
}

choose_adaptive_nbins <- function(expr_data, min_bins = 5, max_bins = 30, min_cells_per_bin = 75) {
  if (is.null(dim(expr_data))) n_total <- length(expr_data) else n_total <- ncol(expr_data)
  non_na_vals <- expr_data[!is.na(expr_data)]
  n <- length(non_na_vals)
  if (n < 2 || all(diff(sort(non_na_vals)) == 0)) return(min_bins)
  sd_val <- sd(expr_data, na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) return(min_bins)
  h <- 3.5 * sd_val / (n_total^(1/3))
  range_data <- diff(range(expr_data, na.rm = TRUE))
  if (is.na(h) || h == 0) return(min_bins)
  nbins <- ceiling(range_data / h)
  if (is.na(nbins) || nbins <= 0 || is.infinite(nbins)) return(min_bins)
  avg_cells_per_bin <- n_total / nbins
  while (avg_cells_per_bin < min_cells_per_bin && nbins > min_bins) {
    nbins <- nbins - 1
    avg_cells_per_bin <- n_total / nbins
  }
  max(min_bins, min(max_bins, nbins))
}

compute_mi_matrix_parallel <- function(expr_subset, n_cores = num_cores) {
  sd_vals <- apply(expr_subset, 1, function(x) sd(x, na.rm = TRUE))
  if (any(is.na(sd_vals) | sd_vals == 0)) {
    cat("Warning: Found genes with NA or zero standard deviation:\n")
    print(names(which(is.na(sd_vals) | sd_vals == 0)))
  }
  nbins_vec <- apply(expr_subset, 1, choose_adaptive_nbins)
  nbins <- round(median(nbins_vec, na.rm = TRUE))
  cat("Using adaptive nbins =", nbins, "\n")
  disc_expr <- infotheo::discretize(t(expr_subset), nbins = nbins)
  n_genes <- nrow(expr_subset)
  gene_pairs <- combn(n_genes, 2, simplify = FALSE)
  mi_values <- mclapply(gene_pairs, function(pair) {
    mutinformation(disc_expr[, pair[1]], disc_expr[, pair[2]])
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  mi_mat <- matrix(0, n_genes, n_genes)
  idx <- combn(n_genes, 2)
  mi_mat[idx] <- mi_mat[t(idx)] <- unlist(mi_values)
  diag(mi_mat) <- 0
  rownames(mi_mat) <- colnames(mi_mat) <- rownames(expr_subset)
  mi_mat
}

get_expr_data <- function(seurat_obj, cluster_id, top_n_genes = 5000, max_cells = 1500, min_expr_pct = 0.05) {
  cluster_cells <- colnames(seurat_obj)[seurat_obj$monocle_clusters2 == cluster_id]
  if (length(cluster_cells) > max_cells) {
    set.seed(123)
    cluster_cells <- sample(cluster_cells, max_cells)
  }
  counts_data <- GetAssayData(seurat_obj[, cluster_cells], slot = "counts", assay = "SCT")
  min_cells <- ceiling(min_expr_pct * ncol(counts_data))
  gene_expr_counts <- rowSums(counts_data > 0)
  genes_to_keep <- names(gene_expr_counts[gene_expr_counts >= min_cells])
  message("Cluster ", cluster_id, ": ", length(genes_to_keep), " genes kept after 5% filter")
  rm(counts_data); gc(verbose = FALSE)
  norm_data <- GetAssayData(seurat_obj[, cluster_cells], slot = "data", assay = "SCT")
  norm_data <- norm_data[genes_to_keep, , drop = FALSE]
  gene_vars <- apply(norm_data, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_n_genes, length(gene_vars))]
  filtered_genes <- filter_housekeeping(top_genes)
  message("Cluster ", cluster_id, ": ", length(filtered_genes), " genes after housekeeping filter")
  as.matrix(norm_data[filtered_genes, , drop = FALSE])
}

build_mi_graph <- function(mi_matrix, threshold_method = "percentile", percentile = 0.95) {
  mi_vals <- mi_matrix[upper.tri(mi_matrix)]
  threshold <- quantile(mi_vals, percentile, na.rm = TRUE)
  adj_matrix <- matrix(0, nrow = nrow(mi_matrix), ncol = ncol(mi_matrix))
  adj_matrix[mi_matrix >= threshold] <- mi_matrix[mi_matrix >= threshold]
  diag(adj_matrix) <- 0
  rownames(adj_matrix) <- colnames(adj_matrix) <- rownames(mi_matrix)
  graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
}

compute_mi_network <- function(seurat_obj, cluster_id, top_n_genes = 5000, max_cells = 1500, n_cores = num_cores) {
  expr_subset <- get_expr_data(seurat_obj, cluster_id, top_n_genes, max_cells)
  mi_mat <- compute_mi_matrix_parallel(expr_subset, n_cores)
  g <- build_mi_graph(mi_mat, threshold_method = "percentile", percentile = 0.95)
  list(graph = g, mi_matrix = mi_mat, genes = rownames(expr_subset))
}

################################################################################
################################################################################
################################################################################
# Compute MI networks
results <- lapply(names(seurat_objects), function(cluster_name) {
  tryCatch({
    cat("Starting", cluster_name, "\n")
    timing <- system.time({
      result <- compute_mi_network(
        seurat_objects[[cluster_name]],
        lineage_clusters[[cluster_name]]$cluster,
        n_cores = num_cores
      )
    })
    print(timing)
    result
  }, error = function(e) {
    message("Error in MI computation for ", cluster_name, ": ", e$message)
    NULL
  })
})
names(results) <- names(seurat_objects)
################################################################################
################################################################################
################################################################################

# Existing L2_C5 Percolation and Community Detection
g_full <- results[["L2_C5"]]$graph
el_full <- igraph::as_data_frame(g_full, what = "edges") %>% dplyr::select(from, to, weight)
top_thr <- quantile(el_full$weight, 0.95)
el_thresh <- el_full %>% filter(weight >= top_thr)
g_thresh <- graph_from_data_frame(el_thresh, directed = FALSE)

################################################################################
#                           graph percolation                                  #
################################################################################
n_iterations <- 100
consensus_edge_list <- list()
for (i in seq_len(n_iterations)) {
  set.seed(321 + i)
  edges_to_keep <- sample(E(g_thresh), size = round(length(E(g_thresh)) * 0.6))
  g_perc_iter <- subgraph.edges(g_thresh, edges_to_keep)
  consensus_edge_list[[i]] <- igraph::as_data_frame(g_perc_iter, what = "edges")
}

combined_edges <- bind_rows(consensus_edge_list) %>%
  group_by(from, to) %>%
  summarise(freq = n(), mean_weight = mean(weight), .groups = "drop")
freq_threshold <- n_iterations * 0.5
consensus_edges <- combined_edges %>% filter(freq >= freq_threshold)
g_consensus_l2_c5 <- graph_from_data_frame(consensus_edges, directed = FALSE)
if("mean_weight" %in% list.edge.attributes(g_consensus_l2_c5)) {
  E(g_consensus_l2_c5)$weight <- E(g_consensus_l2_c5)$mean_weight
}
################################################################################

set.seed(123)
lev_l2_c5 <- cluster_leiden(g_consensus_l2_c5, resolution_parameter = 0.08)
ig_mods_l2_c5 <- split(V(g_consensus_l2_c5)$name, membership(lev_l2_c5))
ig_mods_l2_c5 <- Filter(function(x) length(x) >= 10, ig_mods_l2_c5)
names(ig_mods_l2_c5) <- paste0("IG", seq_along(ig_mods_l2_c5))

################################################################################
#                              label propagation                               #
################################################################################
# Label propagation involves genes adopting community labels based on their 
# neighbors’ labels, weighted by connection strength.
lc_comm_l2_c5 <- cluster_label_prop(g_consensus_l2_c5, weights = E(g_consensus_l2_c5)$weight)
lc_groups_l2_c5 <- split(V(g_consensus_l2_c5)$name, membership(lc_comm_l2_c5))
lc_groups_l2_c5 <- Filter(function(x) length(x) >= 10, lc_groups_l2_c5) #filtered to keep only those with at least 10 genes

################################################################################
#                    Post-Processing with Neighborhood Overlap                 #
################################################################################

adj_matrix_l2_c5 <- as_adjacency_matrix(g_consensus_l2_c5, attr = "weight", sparse = FALSE)
overlap_threshold <- 0.1
for (node in V(g_consensus_l2_c5)$name) {
  idx <- which(V(g_consensus_l2_c5)$name == node)
  neighbors <- V(g_consensus_l2_c5)$name[adj_matrix_l2_c5[idx, ] > 0]
  if (length(neighbors) == 0) next
  for (i in seq_along(lc_groups_l2_c5)) {
    comm_genes <- lc_groups_l2_c5[[i]]
    ov <- length(intersect(neighbors, comm_genes)) / length(neighbors)
    if (ov >= overlap_threshold && !(node %in% comm_genes)) {
      lc_groups_l2_c5[[i]] <- c(comm_genes, node)
    }
  }
}
names(lc_groups_l2_c5) <- paste0("Module ", seq_along(lc_groups_l2_c5))

################################################################################
membership_vector_l2_c5 <- rep(0L, vcount(g_consensus_l2_c5))
names(membership_vector_l2_c5) <- V(g_consensus_l2_c5)$name
for (i in seq_along(lc_groups_l2_c5)) {
  genes <- lc_groups_l2_c5[[i]]
  idx <- match(genes, names(membership_vector_l2_c5))
  valid_idx <- idx[!is.na(idx)]
  membership_vector_l2_c5[valid_idx][membership_vector_l2_c5[valid_idx] == 0L] <- i
}

edges_l2_c5 <- as_edgelist(g_consensus_l2_c5, names = TRUE)
weights_l2_c5 <- E(g_consensus_l2_c5)$weight
m_l2_c5 <- sum(weights_l2_c5)
modularity_manual_l2_c5 <- sum(sapply(unique(membership_vector_l2_c5[membership_vector_l2_c5 > 0]), function(k) {
  ek <- sum(weights_l2_c5[membership_vector_l2_c5[edges_l2_c5[,1]] == k & membership_vector_l2_c5[edges_l2_c5[,2]] == k], na.rm = TRUE)
  dk <- sum(weights_l2_c5[apply(edges_l2_c5, 1, function(e) membership_vector_l2_c5[e[1]] == k | membership_vector_l2_c5[e[2]] == k)], na.rm = TRUE)
  (ek / m_l2_c5) - (dk / (2 * m_l2_c5))^2
}), na.rm = TRUE)

overlap_coeff_l2_c5 <- if (length(lc_groups_l2_c5) > 1) {
  mean(unlist(sapply(1:(length(lc_groups_l2_c5)-1), function(i) {
    sapply((i+1):length(lc_groups_l2_c5), function(j) {
      length(intersect(lc_groups_l2_c5[[i]], lc_groups_l2_c5[[j]])) / 
        min(length(lc_groups_l2_c5[[i]]), length(lc_groups_l2_c5[[j]]))
    })
  })))
} else {
  0
}

################################################################################
# New Percolation and Community Detection for L3_C3
g_full_l3_c3 <- results[["L3_C3"]]$graph
el_full_l3_c3 <- igraph::as_data_frame(g_full_l3_c3, what = "edges") %>% 
  dplyr::select(from, to, weight)
top_thr_l3_c3 <- quantile(el_full_l3_c3$weight, 0.95)
el_thresh_l3_c3 <- el_full_l3_c3 %>% filter(weight >= top_thr_l3_c3)
g_thresh_l3_c3 <- graph_from_data_frame(el_thresh_l3_c3, directed = FALSE)

n_iterations <- 100
consensus_edge_list_l3_c3 <- list()
for (i in seq_len(n_iterations)) {
  set.seed(321 + i)
  edges_to_keep_l3_c3 <- sample(E(g_thresh_l3_c3), size = round(length(E(g_thresh_l3_c3)) * 0.6))
  g_perc_iter_l3_c3 <- subgraph.edges(g_thresh_l3_c3, edges_to_keep_l3_c3)
  consensus_edge_list_l3_c3[[i]] <- igraph::as_data_frame(g_perc_iter_l3_c3, what = "edges")
}

combined_edges_l3_c3 <- bind_rows(consensus_edge_list_l3_c3) %>%
  group_by(from, to) %>%
  summarise(freq = n(), mean_weight = mean(weight), .groups = "drop")
freq_threshold_l3_c3 <- n_iterations * 0.5
consensus_edges_l3_c3 <- combined_edges_l3_c3 %>% filter(freq >= freq_threshold_l3_c3)
g_consensus_l3_c3 <- graph_from_data_frame(consensus_edges_l3_c3, directed = FALSE)
if("mean_weight" %in% list.edge.attributes(g_consensus_l3_c3)) {
  E(g_consensus_l3_c3)$weight <- E(g_consensus_l3_c3)$mean_weight
}

set.seed(123)
lev_l3_c3 <- cluster_leiden(g_consensus_l3_c3, resolution_parameter = 0.1)
ig_mods_l3_c3 <- split(V(g_consensus_l3_c3)$name, membership(lev_l3_c3))
ig_mods_l3_c3 <- Filter(function(x) length(x) >= 10, ig_mods_l3_c3)
names(ig_mods_l3_c3) <- paste0("IG", seq_along(ig_mods_l3_c3))
names(ig_mods_l3_c3)


################################################################################
#                              label propagation                               #
################################################################################
# Label propagation involves genes adopting community labels based on their 
# neighbors’ labels, weighted by connection strength.
lc_comm_l3_c3 <- cluster_label_prop(g_consensus_l3_c3, weights = E(g_consensus_l3_c3)$weight)
lc_groups_l3_c3 <- split(V(g_consensus_l3_c3)$name, membership(lc_comm_l3_c3))
lc_groups_l3_c3 <- Filter(function(x) length(x) >= 10, lc_groups_l3_c3)
names(lc_groups_l3_c3)

adj_matrix_l3_c3 <- as_adjacency_matrix(g_consensus_l3_c3, attr = "weight", sparse = FALSE)
overlap_threshold <- 0.1
for (node in V(g_consensus_l3_c3)$name) {
  idx <- which(V(g_consensus_l3_c3)$name == node)
  neighbors <- V(g_consensus_l3_c3)$name[adj_matrix_l3_c3[idx, ] > 0]
  if (length(neighbors) == 0) next
  for (i in seq_along(lc_groups_l3_c3)) {
    comm_genes <- lc_groups_l3_c3[[i]]
    ov <- length(intersect(neighbors, comm_genes)) / length(neighbors)
    if (ov >= overlap_threshold && !(node %in% comm_genes)) {
      lc_groups_l3_c3[[i]] <- c(comm_genes, node)
    }
  }
}
names(lc_groups_l3_c3) <- paste0("Module ", seq_along(lc_groups_l3_c3))

membership_vector_l3_c3 <- rep(0L, vcount(g_consensus_l3_c3))
names(membership_vector_l3_c3) <- V(g_consensus_l3_c3)$name
for (i in seq_along(lc_groups_l3_c3)) {
  genes <- lc_groups_l3_c3[[i]]
  idx <- match(genes, names(membership_vector_l3_c3))
  valid_idx <- idx[!is.na(idx)]
  membership_vector_l3_c3[valid_idx][membership_vector_l3_c3[valid_idx] == 0L] <- i
}

edges_l3_c3 <- as_edgelist(g_consensus_l3_c3, names = TRUE)
weights_l3_c3 <- E(g_consensus_l3_c3)$weight
m_l3_c3 <- sum(weights_l3_c3)
modularity_manual_l3_c3 <- sum(sapply(unique(membership_vector_l3_c3[membership_vector_l3_c3 > 0]), function(k) {
  ek <- sum(weights_l3_c3[membership_vector_l3_c3[edges_l3_c3[,1]] == k & membership_vector_l3_c3[edges_l3_c3[,2]] == k], na.rm = TRUE)
  dk <- sum(weights_l3_c3[apply(edges_l3_c3, 1, function(e) membership_vector_l3_c3[e[1]] == k | membership_vector_l3_c3[e[2]] == k)], na.rm = TRUE)
  (ek / m_l3_c3) - (dk / (2 * m_l3_c3))^2
}), na.rm = TRUE)


overlap_coeff_l3_c3 <- if (length(lc_groups_l3_c3) > 1) {
  mean(unlist(sapply(1:(length(lc_groups_l3_c3)-1), function(i) {
    sapply((i+1):length(lc_groups_l3_c3), function(j) {
      length(intersect(lc_groups_l3_c3[[i]], lc_groups_l3_c3[[j]])) / 
        min(length(lc_groups_l3_c3[[i]]), length(lc_groups_l3_c3[[j]]))
    })
  })))
} else {
  0  # No overlap possible with only one module
}

################################################################################
#                             Hypergraph creation                              #
################################################################################
# Hypergraph Creation
make_hypergraph <- function(group_list, all_genes) {
  all_genes <- unique(all_genes)
  hyperedges <- lapply(group_list, function(group) {
    valid_genes <- intersect(group, all_genes)
    if (length(valid_genes) == 0) return(NULL)
    Hyperedge(valid_genes)
  })
  hyperedges <- Filter(Negate(is.null), hyperedges)
  hg <- if (length(hyperedges) > 0) Hypergraph(all_genes, hyperedges) else Hypergraph(all_genes, list())
  
  incidence_matrix <- matrix(0, nrow = length(all_genes), ncol = length(group_list))
  rownames(incidence_matrix) <- all_genes
  colnames(incidence_matrix) <- names(group_list)
  for (i in seq_along(group_list)) {
    valid_genes <- intersect(group_list[[i]], all_genes)
    if (length(valid_genes) > 0) {
      incidence_matrix[valid_genes, i] <- 1
    }
  }
  
  list(
    hypergraph = hg,
    incidence_matrix = incidence_matrix,
    genes = all_genes,
    groups = group_list
  )
}

hypergraph_results <- lapply(names(results), function(cluster_name) {
  cat("Working on", cluster_name, "\n")
  res <- results[[cluster_name]]
  gene_groups <- if (cluster_name == "L2_C5") lc_groups_l2_c5 else lc_groups_l3_c3
  names(gene_groups) <- paste0("Module ", seq_along(gene_groups))
  filtered_graph <- if (cluster_name == "L2_C5") g_consensus_l2_c5 else g_consensus_l3_c3
  hypergraph <- make_hypergraph(gene_groups, res$genes)
  list(hypergraph = hypergraph, genes = res$genes, filtered_graph = filtered_graph)
})
names(hypergraph_results) <- names(results)

# Update hypergraph_results for L2_C5 and L3_C3
hypergraph_results[["L2_C5"]]$filtered_graph <- g_consensus_l2_c5
hypergraph_results[["L2_C5"]]$hypergraph <- make_hypergraph(group_list = lc_groups_l2_c5, results[["L2_C5"]]$genes)
hypergraph_results[["L3_C3"]]$filtered_graph <- g_consensus_l3_c3
hypergraph_results[["L3_C3"]]$hypergraph <- make_hypergraph(lc_groups_l3_c3, results[["L3_C3"]]$genes)

#########################################################################

# Cluster Summaries
for (cluster_name in names(hypergraph_results)) {
  if (!is.null(hypergraph_results[[cluster_name]])) {
    hg <- hypergraph_results[[cluster_name]]$hypergraph
    groups <- hg$groups
    num_groups <- length(groups)
    group_sizes <- sapply(groups, length)
    total_genes <- sum(group_sizes)
    cat("Cluster:", cluster_name, "\n")
    cat("  Number of groups:", num_groups, "\n")
    cat("  Group sizes:", paste(group_sizes, collapse = " "), "\n")
    cat("  Total genes in groups:", total_genes, "\n")
    small_graph <- hypergraph_results[[cluster_name]]$filtered_graph
    for (i in seq_along(groups)) {
      subgraph <- induced_subgraph(small_graph, V(small_graph)[name %in% groups[[i]]])
      cat("  Group", i, "size:", length(groups[[i]]), "edges:", gsize(subgraph), "\n")
    }
    avg_size <- mean(group_sizes)
    full_graph <- results[[cluster_name]]$graph
    density <- gsize(full_graph) / (vcount(full_graph) * (vcount(full_graph) - 1) / 2)
    cat("  Avg group size:", round(avg_size, 1), "Density:", round(density, 3), "\n\n")
  }
}

# Overlap Analysis
l2_c5_group1 <- hypergraph_results[["L2_C5"]]$hypergraph$groups[["Module 1"]]
l2_c5_group2 <- hypergraph_results[["L2_C5"]]$hypergraph$groups[["Module 2"]]
l3_c3_group1 <- hypergraph_results[["L3_C3"]]$hypergraph$groups[["Module 1"]]
overlap <- intersect(l2_c5_group1, l3_c3_group1)
unique_l2_c5 <- setdiff(l2_c5_group1, l3_c3_group1)
unique_l3_c3 <- setdiff(l3_c3_group1, l2_c5_group1)
common <- intersect(l2_c5_group1, l3_c3_group1)

# GO Enrichment and Plot Data
plot_data <- data.frame(
  Cluster = character(),
  Group = character(),
  Size = numeric(),
  Density = numeric(),
  Edges = numeric(),
  stringsAsFactors = FALSE
)
enrich_results <- list()

for (cluster_name in names(hypergraph_results)) {
  if (!is.null(hypergraph_results[[cluster_name]])) {
    groups <- hypergraph_results[[cluster_name]]$hypergraph$groups
    small_graph <- hypergraph_results[[cluster_name]]$filtered_graph
    for (i in seq_along(groups)) {
      genes <- groups[[i]]
      filtered_genes <- filter_housekeeping(genes)
      num_nodes <- length(filtered_genes)
      subgraph <- induced_subgraph(small_graph, V(small_graph)[name %in% filtered_genes])
      num_edges <- gsize(subgraph)
      density <- if(num_nodes > 1) num_edges / (num_nodes * (num_nodes - 1) / 2) else NA
      plot_data <- rbind(plot_data, data.frame(Cluster = cluster_name, Group = paste("Group", i), Size = num_nodes, Density = density, Edges = num_edges, stringsAsFactors = FALSE))
      enrich_result <- enrichGO(gene = filtered_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2) %>% clusterProfiler::simplify()
      if (!is.null(enrich_result) && nrow(enrich_result) > 0) {
        enrich_results[[paste(cluster_name, "Group", i, sep = "_")]] <- enrich_result
      }
    }
  }
}

plot_data

############# Bar Plot ###################### Figure 7B
plot_data$Density <- round(plot_data$Density, 2)
plot_data$Group <- str_replace_all(plot_data$Group, pattern = 'Group', replacement = 'Module')
plot1 <- ggplot(subset(plot_data, Cluster %in% c('L2_C5')), aes(x = Group, y = Size, fill = Density)) +
  geom_bar(stat = "identity") +
  scale_y_break(c(40,600), scales = 1) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Density") +
  labs(title = "L2_C5", x = "", y = "Number of Genes") +
  ggpubr::theme_pubr(legend = 'right') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = 'bold'))
plot2 <- ggplot(subset(plot_data, Cluster %in% c('L3_C3')), aes(x = Group, y = Size, fill = Density)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Density") +
  labs(title = "L3_C3", x = "", y = "Number of Genes") +
  ggpubr::theme_pubr(legend = 'right') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = 'bold'))
combined_bar_plot <- plot1 | plot2 + plot_layout(guides = "collect")
combined_bar_plot


filter_graph <- function(graph) {
  deg <- igraph::degree(graph)
  graph <- induced_subgraph(graph, V(graph)[deg >= median(deg)])
  if (vcount(graph) > 2000) {
    set.seed(123)
    sampled_names <- sample(V(graph)$name, 2000)
    graph <- induced_subgraph(graph, V(graph)$name %in% sampled_names)
  }
  graph
}

simplify_graph <- function(g, max_nodes = 50) {
  if (vcount(g) > max_nodes) {
    deg_val <- igraph::degree(g)
    deg <- as.numeric(unlist(deg_val))
    names(deg) <- V(g)$name
    keep <- names(sort(deg, decreasing = TRUE)[1:max_nodes])
    g <- induced_subgraph(g, V(g)[name %in% keep])
  }
  g
}

################################## Figure S9C ##################################
# Network Plots
network_plots <- list()
for (cluster_name in names(hypergraph_results)) {
  if (!is.null(hypergraph_results[[cluster_name]])) {
    groups <- hypergraph_results[[cluster_name]]$hypergraph$groups
    small_graph <- filter_graph(results[[cluster_name]]$graph)
    color <- switch(as.character(cluster_name), "L2_C5" = "#00B6EB", "L3_C3" = "#53B400", "grey")
    for (i in seq_along(groups)) {
      genes <- groups[[i]]
      subgraph <- induced_subgraph(small_graph, V(small_graph)[name %in% genes])
      subgraph <- simplify_graph(subgraph, max_nodes = 50)
      if (vcount(subgraph) == 0) next
      if (gsize(subgraph) > 0) {
        edge_weights <- E(subgraph)$weight
        threshold <- quantile(edge_weights, probs = 0.75)
        subgraph <- delete_edges(subgraph, E(subgraph)[weight < threshold])
      }
      avg_deg <- mean(igraph::degree(subgraph))
      low_deg_threshold <- max(1, avg_deg * 0.25)
      low_deg_nodes <- V(subgraph)[igraph::degree(subgraph) < low_deg_threshold]
      if (length(low_deg_nodes) > 0) {
        subgraph <- delete_vertices(subgraph, low_deg_nodes)
      }
      if (vcount(subgraph) == 0) next
      deg_filtered <- igraph::degree(subgraph)
      V(subgraph)$size <- if (length(deg_filtered) > 0 && max(deg_filtered) > min(deg_filtered)) {
        (deg_filtered - min(deg_filtered)) / (max(deg_filtered) - min(deg_filtered)) * 12 + 3
      } else {
        5
      }
      top_nodes <- names(sort(deg_filtered, decreasing = TRUE)[1:min(10, vcount(subgraph))])
      V(subgraph)$label <- ifelse(V(subgraph)$name %in% top_nodes, V(subgraph)$name, NA)
      p <- ggraph(subgraph, layout = "fr") +
        geom_edge_link(aes(edge_alpha = weight, edge_width = weight), color = "grey60", show.legend = FALSE) +
        geom_node_point(aes(size = size), shape = 21, fill = color, color = "black", show.legend = FALSE) +
        geom_node_text(aes(label = label), fontface = 'bold', repel = TRUE, size = 6, color = "black") +
        scale_size_continuous(range = c(3, 15), guide = "none") +
        scale_edge_alpha(range = c(0.2, 1)) +
        theme_void() +
        labs(title = paste(cluster_name, "Group", i, "Network (", vcount(subgraph), "genes,", gsize(subgraph), "edges)"))
      network_plots[[paste(cluster_name, "Group", i, sep = "_")]] <- p
    }
  }
}

combined_network_plot <- wrap_plots(network_plots, ncol = 3)
combined_network_plot

################################### Figure S9B #################################
# Hub Table
gene_data_list <- lapply(names(hypergraph_results), function(cluster_name) {
  groups <- hypergraph_results[[cluster_name]]$hypergraph$groups
  graph_cluster <- results[[cluster_name]]$graph
  df_list <- lapply(names(groups), function(group_name) {
    genes <- groups[[group_name]]
    subgraph <- induced_subgraph(graph_cluster, V(graph_cluster)[name %in% genes])
    if (vcount(subgraph) == 0) return(NULL)
    deg <- igraph::degree(subgraph)
    data.frame(Cluster = cluster_name, Community = group_name, Gene = names(deg), Score = deg, stringsAsFactors = FALSE)
  })
  do.call(rbind, df_list)
})
gene_data <- do.call(rbind, gene_data_list)
hub_table <- gene_data %>%
  filter(Cluster %in% c('L3_C3', 'L2_C5')) %>%
  group_by(Cluster, Community) %>%
  arrange(desc(Score)) %>%
  slice_head(n = 10) %>%
  summarise(Top_Hub_Genes = paste(Gene, collapse = ", "), .groups = "drop")
hub_table

ft <- flextable(hub_table) %>%
  autofit() %>%
  set_header_labels(Cluster = "Cluster", Community = "Community", Top_Hub_Genes = "Top Hub Genes") %>%
  bg(part = "header", bg = "gray") %>%
  merge_v(j = "Cluster") %>%
  align(j = "Cluster", align = "center", part = "body") %>%
  theme_vanilla()

ft

################################ Figure 7C ######################################
# Venn Diagram
venn.plot <- venn.diagram(
  x = list(L2_C5 = l2_c5_group1, L3_C3 = l3_c3_group1),
  disable.logging = T,
  filename = NULL,
  output = TRUE,
  fill = c("#00B6EB", "#53B400"),
  alpha = 0.5,
  label.col = "black",
  cex = 2.5,
  fontfamily = "sans",
  cat.cex = 2,
  cat.fontfamily = "sans",
  cat.col = c("#00B6EB", "#53B400"),
  cat.pos = c(-40, 40),
  cat.dist = c(0.05, 0.05)
)
grid.draw(venn.plot)
dev.off()

################################################################################
############################### Figure 7D ######################################
# L2_C5 Specific Analysis
groups <- hypergraph_results[["L2_C5"]]$hypergraph$groups
gene_memberships <- table(unlist(groups))
overlapping_genes <- names(gene_memberships[gene_memberships > 1])
if (length(overlapping_genes) >= 3) {
  enr <- enrichGO(gene = overlapping_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  if (!is.null(enr) && nrow(enr@result) > 0) {
    enr_simp <- if (nrow(enr@result) > 1) clusterProfiler::simplify(enr) else enr
  }
}

edges <- data.frame()
for (i in seq_along(groups)) {
  edges <- rbind(edges, data.frame(Gene = groups[[i]], Community = paste("Module", i, sep = "")))
}
g_viz <- graph_from_data_frame(edges, directed = FALSE)
V(g_viz)$type <- V(g_viz)$name %in% paste0("Module", seq_along(groups))
V(g_viz)$color <- ifelse(V(g_viz)$name %in% overlapping_genes, "red", ifelse(V(g_viz)$type, "green", "blue"))
hypergraph_plot <- ggraph(g_viz, layout = "fr") +
  geom_edge_link(color = "grey") +
  geom_node_point(aes(color = color, size = ifelse(type, 5, 3))) +
  geom_node_text(aes(label = ifelse(type | name %in% c(overlapping_genes, 'PRSS1', 'PRSS2'), name, "")), repel = TRUE, size = 3) +
  scale_color_identity(guide = "legend", labels = c("Gene", "Community", "Overlapping Gene"), breaks = c("blue", "green", "red")) +
  scale_size_identity() +
  theme_void() +
  labs(title = "L2_C5 Hypergraph with Overlapping Communities")

hypergraph_plot

# Regulon Analysis
regulonTargetsInfo_l2_c5 <- readRDS('L2_C5_regulonTargetsInfo.Rds')
regulon_overlap <- regulonTargetsInfo_l2_c5[regulonTargetsInfo_l2_c5$gene %in% overlapping_genes, ]
tf_summary_l2c5 <- regulon_overlap %>%
  group_by(TF) %>%
  summarise(Targets = paste(gene, collapse = ", "), 
            Num_Targets = n(), 
            Avg_CoexWeight = mean(Genie3Weight), 
            HighConfAnnot = any(highConfAnnot == "TRUE")) %>%
  arrange(desc(Num_Targets), desc(Avg_CoexWeight))

################################################################################
############################### Figure 7F ######################################

exocrine_markers <- c("CTRB2", "REG1A", "REG1B",
                      "SPINK1","PRSS1","PRSS2")

# Get SCENIC-only candidate edges (marker-filtered)
chord_data_l2_c5 <- regulonTargetsInfo_l2_c5 %>%
  filter(gene %in% exocrine_markers, NES >= 1.65) %>%
  dplyr::rename(TargetGene = gene) %>%
  mutate(Cluster = "L2_C5") %>%
  dplyr::select(TF, TargetGene, NES)

dim(chord_data_l2_c5)

# Get MI edges for L2_C5
mi_edges_l2_c5 <- as.data.frame(as_edgelist(results[["L2_C5"]]$graph, 
                                            names = TRUE))
colnames(mi_edges_l2_c5) <- c("Gene1", "Gene2")
mi_edges_l2_c5$MI_score <- E(results[["L2_C5"]]$graph)$weight

dim(mi_edges_l2_c5)

# Filter MI network for non-housekeeping genes
mi_genes_l2_c5 <- unique(c(mi_edges_l2_c5$Gene1, mi_edges_l2_c5$Gene2))
mi_genes_l2_c5 <- filter_housekeeping(mi_genes_l2_c5)
dim(mi_edges_l2_c5)

# Create MI pairs as strings (original)
mi_pairs_raw <- unique(
  apply(mi_edges_l2_c5, 1, function(row) paste(sort(row), collapse = "_"))
)

# Clean MI pairs by extracting gene names (assuming format: "<numeric>_<gene1>_<gene2>")
clean_mi_pairs <- sapply(mi_pairs_raw, function(x) {
  parts <- unlist(strsplit(x, "_"))
  if(length(parts) >= 3){
    # Use the 2nd and 3rd elements as the gene names
    sorted_pair <- paste(sort(parts[2:3]), collapse = "_")
  } else {
    sorted_pair <- x
  }
  return(sorted_pair)
})
clean_mi_pairs <- unique(clean_mi_pairs)
head(clean_mi_pairs)

# Annotate each SCENIC edge with MI evidence using cleaned MI pairs.
# For each SCENIC edge, create a pair string "TF_TargetGene" (sorted alphabetically)
chord_data_l2_c5 <- chord_data_l2_c5 %>%
  mutate(pair = paste0(pmin(TF, TargetGene), "_", pmax(TF, TargetGene)),
         Source = ifelse(pair %in% clean_mi_pairs, 
                         "SCENIC + MI", "SCENIC only"))

# Debug: Check intersection
print("MI-supported pairs (cleaned):")
print(intersect(chord_data_l2_c5$pair, clean_mi_pairs))
print(table(chord_data_l2_c5$Source))

# Define link colors based on evidence
# Green for SCENIC + MI; Orange for SCENIC only.
link_colors <- ifelse(chord_data_l2_c5$Source == "SCENIC + MI", "#1b9e77", "#d95f02")
# Assign names in the same order as rows in chord_data_l2_c5.
link_colors <- setNames(link_colors, paste(chord_data_l2_c5$TF, chord_data_l2_c5$TargetGene, sep = "_"))

# Define sector colors (TFs in tomato; targets in lightblue)
sectors <- unique(c(chord_data_l2_c5$TF, chord_data_l2_c5$TargetGene))
tf_colors <- rep("tomato", length(unique(chord_data_l2_c5$TF)))
names(tf_colors) <- unique(chord_data_l2_c5$TF)
target_colors <- rep("lightblue", length(setdiff(sectors, names(tf_colors))))
names(target_colors) <- setdiff(sectors, names(tf_colors))
grid.col <- c(tf_colors, target_colors)

# Clear circos settings and plot the chord diagram
circos.clear()
circos.par(start.degree = 90, gap.degree = 1)

chordDiagram(
  chord_data_l2_c5[, c("TF", "TargetGene", "NES")],
  grid.col = grid.col,
  col = link_colors,
  transparency = 0.3,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.2),
  directional = 1,
  direction.type = "arrows",
  diffHeight = -0.02,
  link.arr.type = "big.arrow"
)

# Add sector labels (highlight marker genes and TFs)
highlight_genes <- exocrine_markers
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  label <- CELL_META$sector.index
  if (label %in% highlight_genes | label %in% unique(chord_data_l2_c5$TF)) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] - 1, label,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.9)
  }
}, bg.border = NA)

# title("SCENIC–MI TF–Target Interactions (L2_C5)")

# Add a manual legend
legend("topright",
       legend = c("SCENIC only", "SCENIC + MI"),
       col = c("#d95f02", "#1b9e77"),
       lwd = 2,
       title = expression(bold("Edge Source")),
       bty = "n",
       cex = 0.8)
dev.off()

################################################################################
############################### Figure 7G ######################################

## --------------------------------------------------------------------------
## pseudotime columns & cluster-ID map
## --------------------------------------------------------------------------
seurat_objects$L2_C5$pseudotime <- seurat_objects$L2_C5$Lineage2_2
cluster_id_map <- c(L2_C5 = "5")

###############################################################################
## Sanity checks & visualizations
###############################################################################
## Module-score trajectories along pseudotime ------------------------------
# define modules
mods <- list(
  L2_C5_JUND = c('JUND'),
  L2_C5_M2_prss = c('PRSS1','PRSS2'),
  L2_C5_M2_bridge = c('SPINK1','REG1A','REG1B','CTRB2')
)
# compute scores
for (modName in names(mods)) {
  # extract cluster key
  parts       <- strsplit(modName, "_")[[1]]
  clusterName <- paste(parts[1:2], collapse = "_")
  so <- seurat_objects[[clusterName]]
  DefaultAssay(so) <- "SCT"
  so <- AddModuleScore(so, features = list(mods[[modName]]), name = modName)
  seurat_objects[[clusterName]] <- so
}


traj_df <- bind_rows(lapply(names(mods), function(modName) {
  parts       <- strsplit(modName, "_")[[1]]
  clusterName <- paste(parts[1:2], collapse = "_")
  so          <- seurat_objects[[clusterName]]
  
  # only cells in that monocle cluster and with non-NA pseudotime
  keep_cells <- colnames(so)[
    so$monocle_clusters2 == cluster_id_map[clusterName] &
      !is.na(so$pseudotime)
  ]
  score_col <- paste0(modName, "1")
  
  data.frame(
    Cell       = keep_cells,
    Cluster    = clusterName,
    Module     = modName,
    Pseudotime = so$pseudotime[keep_cells],
    Score      = so@meta.data[keep_cells, score_col],
    stringsAsFactors = FALSE
  )
}))

module_labels <- c(
  "L2_C5_JUND" = "JUND Signature",
  "L2_C5_M2_prss" = "PRSS1/2 Signature",
  "L2_C5_M2_bridge" = "SPINK1|REG1A/B|CTRB2 Signature"
)

module_color <- list(
  "L2_C5_JUND" = "#FEAF16FF",
  "L2_C5_M2_prss" = "#3283FEFF",
  "L2_C5_M2_bridge" = "#F8766D"
)

p_modtraj <- ggplot(traj_df, aes(x = Pseudotime, y = Score, 
                                 color = Module, fill  = Module)) +
  geom_smooth(se = TRUE, span = 0.5) +
  facet_wrap(~ Cluster, scales = "free_x") +
  labs(title = "Module-score trajectories along pseudotime",
       x = "Pseudotime", y = "Module score") +
  theme_pubr(legend = 'right') + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 'bold', size = 18))+
  scale_color_manual(values = module_color, labels = module_labels) + # Define color legend labels
  scale_fill_manual(values = module_color, labels = module_labels)   # Define fill legend labels

print(p_modtraj)

################################################################################
# sessionInfo()
################################################################################
