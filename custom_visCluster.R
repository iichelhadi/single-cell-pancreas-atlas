# custom_visCluster.R

# Improved visCluster_custom with support for row_names_side and returning Heatmap object
visCluster_custom <- function(object = NULL,
                              ht.col.list = list(
                                col_range = c(-2, 0, 2),
                                col_color = c("#53B400", "#00FFFF", "#A58AFF")
                              ),
                              border = TRUE,
                              top_annotation = NULL,
                              row_names = NULL,
                              row_label_size = 10,
                              row_names_side = c("right", "left")[1],
                              raster_device = "png",
                              return_ht = FALSE) {
  # Load required packages
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) install.packages("ComplexHeatmap")
  if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
  if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")
  
  ComplexHeatmap::ht_opt(message = FALSE)
  
  # Precompute the color mapping function
  col_fun <- circlize::colorRamp2(ht.col.list$col_range, ht.col.list$col_color)
  
  # Extract and convert data to a dense, numeric matrix
  data <- object$wide.res
  mat <- data.matrix(data[, -(1:2)])
  rownames(mat) <- data$gene
  
  # Create heatmap object
  htf <- ComplexHeatmap::Heatmap(
    mat,
    name             = "Z-score",
    cluster_columns  = FALSE,
    show_row_names =  FALSE,
    # show_row_names   = !is.null(row_names),
    row_names_side   = row_names_side,
    show_row_dend    = FALSE,
    show_column_names= FALSE,
    border           = border,
    col              = col_fun,
    top_annotation   = top_annotation,
    use_raster       = TRUE,
    raster_device    = raster_device
  )
  
  # Add gene markers if requested
  if (!is.null(row_names)) {
    gene_indices <- match(row_names, rownames(mat))
    valid_indices <- !is.na(gene_indices)
    if (any(valid_indices)) {
      htf <- htf + ComplexHeatmap::rowAnnotation(
        mark = ComplexHeatmap::anno_mark(
          at = gene_indices[valid_indices],
          labels = row_names[valid_indices],
          labels_gp = grid::gpar(fontsize = row_label_size)
        )
      )
    }
  }
  
  # Return or draw
  if (return_ht) {
    return(htf)
  } else {
    ComplexHeatmap::draw(htf, merge_legend = FALSE)
  }
}
