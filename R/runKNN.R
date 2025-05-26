#' Run K-Nearest Neighbors for TLS Identification
#'
#' Identifies tumor-associated lymphoid structures (TLS) in spatial transcriptomics data
#' using a K-nearest neighbors approach based on expression thresholds and spatial proximity.
#'
#' @param seurat A \code{Seurat} object containing spatial transcriptomics data.
#' @param exp_threshold A numeric value (default: 0.70) indicating the expression threshold
#'        for identifying high-expression spots.
#' @param min_spots An integer (default: 3) representing the minimum number of adjacent spots
#'        required to form a valid TLS cluster.
#' @param max_distance A numeric value (default: 8) defining the maximum distance between
#'        adjacent spots to be considered neighbors.
#'
#' @return A modified \code{Seurat} object with an additional metadata column \code{TLS_identity}
#'         indicating whether each spot is part of a TLS or not.
#'
#' @details This function processes each sample in the provided \code{Seurat} object to:
#' \enumerate{
#'   \item Subset data to individual samples.
#'   \item Identify high expression spots based on the specified threshold.
#'   \item Calculate dynamic neighbors using a K-nearest neighbors approach.
#'   \item Create an adjacency matrix for identified neighbors.
#'   \item Find connected components in the adjacency matrix to identify valid TLS clusters.
#'   \item Label spots as "TLS" or "not TLS" based on their cluster membership.
#' }
#'
#' @examples
#' \donttest{
#' seurat_obj <- runKNN(seurat, exp_threshold = 0.70, min_spots = 3, max_distance = 8)
#' }
#'
#' @importFrom Seurat AddMetaData
#' @importFrom dplyr %>%
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_adjacency_matrix components
#' @export
runKNN <- function(seurat,
                   exp_threshold = 0.70,
                   min_spots = 3,     # Minimum number of adjacent spots
                   max_distance = 8) {   # Maximum distance between spots

  message("Now, for each sample...")
  samples <- seurat$orig.ident %>% unique

  myList <- lapply(samples, function(x) {
    tryCatch({
      # Subsetting Seurat Object
      message(paste0("\nProcessing sample: ", x))
      seurat.sub <- subset(seurat, subset = orig.ident == x)
      seurat.sub <- AddMetaData(seurat.sub, GetTissueCoordinates(seurat.sub))
      df <- seurat.sub@meta.data

      # Find spots with high expression
      high_exp_spots <- df$TLS >= exp_threshold
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
        valid_neighbors <- neighbors[distances < max_distance]
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



