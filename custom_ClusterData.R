### custom_ClusterData.R

library(data.table)

filter.std <- function(eset = NULL, min.std = 0, visu = FALSE) {
  if (!requireNamespace("matrixStats", quietly = TRUE)) 
    install.packages("matrixStats")
  
  # 10x faster SD calculation
  row_sds <- matrixStats::rowSds(eset, na.rm = TRUE)  
  
  index <- row_sds > min.std
  index[is.na(index)] <- FALSE  # Exclude NAs instead of keeping
  
  cat(paste(sum(!index), "genes excluded.\n"))
  
  if (visu) plot(sort(row_sds), xlab = "Ordered genes", ylab = "Sd")
  
  eset[index, , drop = FALSE]
}

clusterData_custom <- function(exp = NULL, scaleData = FALSE, 
                               min.std = 0, seed = 123) {
  if (!requireNamespace("fastcluster", quietly = TRUE)) 
    install.packages("fastcluster")
  if (!requireNamespace("Rfast", quietly = TRUE)) 
    install.packages("Rfast")
  
  # Filter rows efficiently
  exp <- filter.std(exp, min.std = min.std)
  
  # Ensure input is a matrix and optionally scale data
  hclust_matrix <- if (scaleData) t(scale(t(exp))) else as.matrix(exp)
  
  # Fast distance computation using Rfast
  dist_matrix <- Rfast::Dist(hclust_matrix, method = "euclidean")
  
  # Convert dist_matrix to a valid distance object
  dist_object <- as.dist(dist_matrix)
  
  # Perform hierarchical clustering using fastcluster
  hc <- fastcluster::hclust(dist_object)
  
  # Get ordered indices
  ordered_indices <- hc$order
  
  # Prepare wide format data
  wide.r <- data.table::data.table(
    gene = rownames(hclust_matrix)[ordered_indices],
    cluster = 1,  # Single cluster due to row_km = 1
    as.data.table(hclust_matrix[ordered_indices, ])
  )
  
  # Convert to long format using data.table
  df <- data.table::melt(
    wide.r,
    id.vars = c("cluster", "gene"),
    variable.name = "cell_type",
    value.name = "norm_value"
  )
  
  # Add cluster name
  df[, cluster_name := sprintf("cluster 1 (%d)", nrow(wide.r))]
  
  return(list(
    wide.res = wide.r,
    long.res = df,
    type = "hierarchical",
    geneMode = "none",
    geneType = "none"
  ))
}
