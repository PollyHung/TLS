#' Run K-Nearest Neighbors for TLS Identification
#'
#' Identifies tumor-associated lymphoid structures (TLS) in spatial transcriptomics data
#' using an adaptive K-nearest neighbors approach that automatically adjusts to each sample's
#' spatial characteristics.
#'
#' @param seurat A \code{Seurat} object containing spatial transcriptomics data with
#'        coordinates stored in the \code{images} slot.
#' @param exp_threshold A numeric value (default: 0.98) indicating the quantile threshold
#'        for identifying high-expression spots (e.g., 0.98 = top 2% of spots).
#' @param min_spots An integer (default: 3) representing the minimum number of adjacent spots
#'        required to form a valid TLS cluster.
#' @param max_distance An optional numeric value defining the absolute maximum distance between
#'        adjacent spots (in coordinate units). If NULL (default), calculates automatically.
#' @param distance_multiplier A numeric value (default: 2) scaling factor applied to the
#'        median nearest-neighbor distance when auto-calculating \code{max_distance}.
#'
#' @return A modified \code{Seurat} object with:
#' \itemize{
#'   \item New metadata column \code{TLS_identity} ("TLS" or "not TLS")
#'   \item Preserved all original data and reductions
#' }
#'
#' @details The function performs adaptive TLS detection by:
#' \enumerate{
#'   \item Calculating sample-specific distance thresholds based on spot spacing
#'   \item Identifying high-expression spots using the specified quantile threshold
#'   \item Building KNN graphs with dynamic neighbor counts (k = min(6, n_spots-1))
#'   \item Filtering neighbors by adaptive spatial constraints
#'   \item Identifying connected components as TLS candidates
#'   \item Applying size thresholds to validate TLS clusters
#' }
#'
#' @section Adaptive Distance Calculation:
#' When \code{max_distance = NULL}, the threshold is calculated as:
#' \deqn{max\_distance = distance\_multiplier \times median(nearest\_neighbor\_distances)}
#' This ensures automatic adjustment for samples with different spatial resolutions.
#'
#' @examples
#' \donttest{
#' # Automatic distance calculation (recommended for multi-sample datasets)
#' seurat <- runKNN(seurat, exp_threshold = 0.98, distance_multiplier = 2.5)
#'
#' # Manual distance threshold
#' seurat <- runKNN(seurat, exp_threshold = 0.95, max_distance = 5)
#' }
#'
#' @importFrom Seurat GetTissueCoordinates AddMetaData
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom stats median quantile
#' @export
runKNN <- function(seurat,
                   exp_threshold = 0.98,
                   min_spots = 3,
                   max_distance = NULL,  # Now optional
                   distance_multiplier = 2) {   # Maximum distance between spots

  calculate_adaptive_distance <- function(coords) {
    nn1 <- RANN::nn2(coords[1:2], k = 3)$nn.dists[,2]
    median_dist <- median(nn1)
    return(median_dist * distance_multiplier)
  }

  message("Now, for each sample...")
  samples <- unique(seurat$orig.ident)

  myList <- lapply(samples, function(x) {
    tryCatch({
      # Subsetting Seurat Object
      message(paste0("\nProcessing sample: ", x))
      seurat.sub <- subset(seurat, subset = orig.ident == x)
      seurat.sub <- AddMetaData(seurat.sub, GetTissueCoordinates(seurat.sub))
      coords <- GetTissueCoordinates(seurat.sub)

      # Calculate min distance
      if(is.null(max_distance)) {
        sample_max_dist <- calculate_adaptive_distance(coords)
        message("Auto-calculated max_distance: ", round(sample_max_dist, 2))
      } else {
        sample_max_dist <- max_distance
      }

      # Keep on
      df <- cbind(coords, seurat.sub@meta.data["TLS"])
      colnames(df) <- c("x", "y", "cell", "TLS_score")

      # Find spots with high expression
      exp_threshold <- quantile(df$TLS_score, probs = exp_threshold)
      high_exp_spots <- df$TLS_score >= exp_threshold
      n_high <- sum(high_exp_spots)
      message(paste0("Found ", n_high, " high expression spots"))

      # Early return if no high-exp spots
      if(n_high == 0) {
        message("No high expression spots - skipping neighborhood analysis")
        seurat.sub$TLS_identity <- "not TLS"
        return(seurat.sub)
      }

      # Dynamic neighbor calculation
      k <- min(6, n_high - 1)  # Ensure k <= (n-1)
      if(k < 1) {
        message("Insufficient spots for neighbor analysis")
        seurat.sub$TLS_identity <- "not TLS"
        return(seurat.sub)
      }

      # Build nearest neighbor graph
      message("Calculating nearest neighbors with k = ", k)
      nn <- RANN::nn2(df[high_exp_spots, c("x", "y")], k = k)

      # Create adjacency matrix
      adj_matrix <- matrix(0, n_high, n_high)
      for(i in 1:n_high) {
        neighbors <- nn$nn.idx[i, ]
        distances <- nn$nn.dists[i, ]
        valid_neighbors <- neighbors[distances < sample_max_dist]
        adj_matrix[i, valid_neighbors] <- 1
      }

      # Find connected components
      g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
      components <- igraph::components(g)

      # Identify valid TLS
      valid_tls <- which(components$csize >= min_spots)
      if(length(valid_tls) == 0) {
        message("No TLS clusters found meeting size threshold")
        seurat.sub$TLS_identity <- "not TLS"
        return(seurat.sub)
      }

      # Create TLS labels
      tls_labels <- rep(0, ncol(seurat.sub))
      for(i in valid_tls) {
        component_spots <- which(high_exp_spots)[components$membership == i]
        tls_labels[component_spots] <- i
      }

      seurat.sub$TLS_identity <- ifelse(tls_labels == 0, "not TLS", "TLS")
      return(seurat.sub)

    }, error = function(e) {
      message("\nError processing ", x, ": ", e$message)
      # Return object with default TLS identity
      seurat.sub$TLS_identity <- "not TLS"
      return(seurat.sub)
    })
  })
  seurat2 <- Reduce(function(x, y) merge(x, y), myList)
  return(seurat2)
}



