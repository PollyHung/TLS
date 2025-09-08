#' Run K-Nearest Neighbors for TLS Identification
#'
#' Identifies tumor-associated lymphoid structures (TLS) in spatial transcriptomics data
#' using an adaptive K-nearest neighbors approach that automatically adjusts to each sample's
#' spatial characteristics.
#'
#' @param seurat A \code{Seurat} object containing spatial transcriptomics data with
#'        coordinates stored in the \code{images} slot.
#' @param exp_threshold A numeric value (default: 0.98) indicating the quantile threshold
#'        for identifying high-expression spots when use_statistical_threshold = FALSE.
#' @param min_spots An integer (default: 3) representing the minimum number of adjacent spots
#'        required to form a valid TLS cluster.
#' @param max_distance An optional numeric value defining the absolute maximum distance between
#'        adjacent spots (in coordinate units). If NULL (default), calculates automatically.
#' @param distance_multiplier A numeric value (default: 2) scaling factor applied to the
#'        median nearest-neighbor distance when auto-calculating \code{max_distance}.
#' @param use_statistical_threshold Logical (default: TRUE). Use statistical significance
#'        testing instead of quantile-based thresholding.
#' @param alpha Significance level (default: 0.05) for statistical thresholding.
#' @param correction_method Multiple testing correction method ("fdr", "bonferroni").
#'
#' @return A modified \code{Seurat} object with:
#' \itemize{
#'   \item New metadata column \code{TLS_identity} ("TLS" or "not TLS")
#'   \item New metadata column \code{TLS_confidence} (0-1 confidence score)
#'   \item New metadata column \code{TLS_cluster_id} (cluster identifier for TLS spots)
#'   \item New metadata column \code{spatial_coherence} (spatial clustering strength)
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
#' @importFrom stats median quantile pnorm p.adjust
#' @export

#' Calculate statistically rigorous expression thresholds
#' 
#' @param expression_values Vector of expression values for TLS signature
#' @param alpha Significance level (default: 0.05)
#' @param method Multiple testing correction method ("fdr", "bonferroni")
#' @return List with threshold value and statistical measures
calculateStatisticalThreshold <- function(expression_values, alpha = 0.05, method = "fdr") {
  # Remove missing values
  exp_vals <- expression_values[!is.na(expression_values)]
  n <- length(exp_vals)
  
  if (n < 10) {
    warning("Small sample size (n < 10). Results may be unreliable.")
    return(list(
      threshold = quantile(exp_vals, 0.98),
      method = "quantile_fallback",
      p_values = NA,
      adjusted_p_values = NA
    ))
  }
  
  # Calculate mean and standard deviation
  mean_exp <- mean(exp_vals)
  sd_exp <- sd(exp_vals)
  
  # Calculate z-scores for each spot
  z_scores <- (exp_vals - mean_exp) / sd_exp
  
  # Calculate p-values assuming normal distribution
  p_values <- 1 - pnorm(z_scores)
  
  # Apply multiple testing correction
  adjusted_p_values <- p.adjust(p_values, method = method)
  
  # Find threshold corresponding to significance level
  significant_spots <- adjusted_p_values < alpha
  
  if (sum(significant_spots) == 0) {
    # If no spots are significant, use less stringent threshold
    threshold <- quantile(exp_vals, 1 - alpha)
    warning(paste("No statistically significant spots found. Using", 
                  round((1-alpha)*100), "th percentile as threshold."))
  } else {
    # Use minimum expression value among significant spots
    threshold <- min(exp_vals[significant_spots])
  }
  
  return(list(
    threshold = threshold,
    method = paste0("statistical_", method),
    p_values = p_values,
    adjusted_p_values = adjusted_p_values,
    n_significant = sum(significant_spots),
    alpha = alpha
  ))
}
runKNN <- function(seurat,
                   exp_threshold = 0.98,
                   min_spots = 3,
                   max_distance = NULL,
                   distance_multiplier = 2,
                   use_statistical_threshold = TRUE,
                   alpha = 0.05,
                   correction_method = "fdr") {

  # Validate input data
  input_validation <- validateSeuratInput(seurat, require_spatial = TRUE, require_tls_score = TRUE)
  if (!input_validation$valid) {
    printValidationErrors(input_validation, "input")
    stop("Input validation failed. See errors above.")
  }
  
  # Print warnings if any
  if (length(input_validation$warnings) > 0) {
    for (warning_msg in input_validation$warnings) {
      warning(warning_msg)
    }
  }
  
  # Validate parameters
  param_validation <- validateTLSParameters(
    exp_threshold = exp_threshold,
    min_spots = min_spots,
    max_distance = max_distance,
    distance_multiplier = distance_multiplier,
    alpha = alpha,
    use_statistical_threshold = use_statistical_threshold
  )
  
  if (!param_validation$valid) {
    printValidationErrors(param_validation, "parameters")
    stop("Parameter validation failed. See errors above.")
  }
  
  # Print parameter warnings if any
  if (length(param_validation$warnings) > 0) {
    for (warning_msg in param_validation$warnings) {
      warning(warning_msg)
    }
  }

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

      # Find spots with high expression using statistical or quantile thresholding
      if (use_statistical_threshold) {
        stat_result <- calculateStatisticalThreshold(df$TLS_score, alpha = alpha, method = correction_method)
        threshold_value <- stat_result$threshold
        message(paste0("Using statistical threshold: ", round(threshold_value, 4), 
                      " (method: ", stat_result$method, ", n_significant: ", stat_result$n_significant, ")"))
      } else {
        threshold_value <- quantile(df$TLS_score, probs = exp_threshold)
        message(paste0("Using quantile threshold (", exp_threshold*100, "th percentile): ", round(threshold_value, 4)))
      }
      
      high_exp_spots <- df$TLS_score >= threshold_value
      n_high <- sum(high_exp_spots)
      message(paste0("Found ", n_high, " high expression spots"))

      # Early return if no high-exp spots
      if(n_high == 0) {
        message("No high expression spots - skipping neighborhood analysis")
        seurat.sub$TLS_identity <- "not TLS"
        seurat.sub$TLS_confidence <- 0
        seurat.sub$TLS_cluster_id <- 0
        seurat.sub$spatial_coherence <- 0
        return(seurat.sub@meta.data)
      }

      # Dynamic neighbor calculation
      k <- min(6, n_high - 1)  # Ensure k <= (n-1)
      if(k < 1) {
        message("Insufficient spots for neighbor analysis")
        seurat.sub$TLS_identity <- "not TLS"
        seurat.sub$TLS_confidence <- 0
        seurat.sub$TLS_cluster_id <- 0
        seurat.sub$spatial_coherence <- 0
        return(seurat.sub@meta.data)
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
        seurat.sub$TLS_confidence <- 0
        seurat.sub$TLS_cluster_id <- 0
        seurat.sub$spatial_coherence <- 0
        return(seurat.sub@meta.data)
      }

      # Initialize confidence scores and cluster IDs
      tls_labels <- rep(0, ncol(seurat.sub))
      tls_confidence <- rep(0, ncol(seurat.sub))
      spatial_coherence <- rep(0, ncol(seurat.sub))
      
      # Calculate confidence scores for each valid TLS cluster
      for(i in valid_tls) {
        component_spots_idx <- which(components$membership == i)
        component_spots <- which(high_exp_spots)[component_spots_idx]
        tls_labels[component_spots] <- i
        
        # Calculate confidence based on multiple factors
        cluster_size <- components$csize[i]
        size_score <- min(cluster_size / 10, 1)  # Size factor (saturates at 10 spots)
        
        # Expression strength score (normalized TLS scores within cluster)
        cluster_exp_scores <- df$TLS_score[component_spots]
        exp_strength <- mean(cluster_exp_scores) / max(df$TLS_score)
        
        # Spatial coherence score based on average distance to cluster centroid
        if(length(component_spots) > 1) {
          cluster_coords <- df[component_spots, c("x", "y")]
          centroid <- c(mean(cluster_coords$x), mean(cluster_coords$y))
          distances_to_centroid <- sqrt(rowSums(sweep(cluster_coords, 2, centroid)^2))
          spatial_coherence_score <- 1 / (1 + mean(distances_to_centroid) / sample_max_dist)
        } else {
          spatial_coherence_score <- 0.5  # Single spot gets moderate coherence
        }
        
        # Combined confidence score (weighted average)
        confidence <- 0.4 * exp_strength + 0.3 * size_score + 0.3 * spatial_coherence_score
        
        tls_confidence[component_spots] <- confidence
        spatial_coherence[component_spots] <- spatial_coherence_score
      }

      seurat.sub$TLS_identity <- ifelse(tls_labels == 0, "not TLS", "TLS")
      seurat.sub$TLS_confidence <- tls_confidence
      seurat.sub$TLS_cluster_id <- tls_labels
      seurat.sub$spatial_coherence <- spatial_coherence
      
      return(seurat.sub@meta.data)

    }, error = function(e) {
      message("\nError processing ", x, ": ", e$message)
      # Return object with default values
      seurat.sub$TLS_identity <- "not TLS"
      seurat.sub$TLS_confidence <- 0
      seurat.sub$TLS_cluster_id <- 0
      seurat.sub$spatial_coherence <- 0
      return(seurat.sub@meta.data)
    })
  })
  new_meta <- do.call(rbind, myList)
  new_meta <- new_meta[Cells(seurat), ]
  seurat@meta.data <- new_meta
  return(seurat)
}



