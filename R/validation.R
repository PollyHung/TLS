#' Validate input data for TLS analysis
#' 
#' Comprehensive validation of Seurat objects before TLS detection analysis.
#' Checks data structure, spatial coordinates, and data quality.
#'
#' @param seurat A \code{Seurat} object to validate
#' @param min_cells Minimum number of cells required (default: 10)
#' @param min_features Minimum number of features required (default: 100)
#' @param require_spatial Logical whether spatial coordinates are required (default: TRUE)
#' @param require_tls_score Logical whether TLS scores are required (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{valid}: Logical indicating if object passes validation
#'   \item \code{errors}: Character vector of error messages
#'   \item \code{warnings}: Character vector of warning messages
#'   \item \code{summary}: List of data quality metrics
#' }
#'
#' @details The validation checks:
#' \itemize{
#'   \item Basic Seurat object structure and required slots
#'   \item Presence and validity of spatial coordinates
#'   \item Data completeness and quality metrics
#'   \item TLS signature scores if required
#'   \item Sample identifiers and metadata structure
#' }
#'
#' @examples
#' \donttest{
#' validation_result <- validateSeuratInput(seurat_obj)
#' if (!validation_result$valid) {
#'   stop("Validation failed: ", paste(validation_result$errors, collapse = "; "))
#' }
#' }
#'
#' @importFrom Seurat GetTissueCoordinates
#' @export
validateSeuratInput <- function(seurat,
                                min_cells = 10,
                                min_features = 100,
                                require_spatial = TRUE,
                                require_tls_score = TRUE) {
  
  errors <- character(0)
  warnings <- character(0)
  
  # Basic object validation
  if (!inherits(seurat, "Seurat")) {
    errors <- c(errors, "Input is not a Seurat object")
    return(list(valid = FALSE, errors = errors, warnings = warnings, summary = NULL))
  }
  
  # Check basic dimensions
  n_cells <- ncol(seurat)
  n_features <- nrow(seurat)
  
  if (n_cells < min_cells) {
    errors <- c(errors, paste("Insufficient number of cells:", n_cells, "<", min_cells))
  }
  
  if (n_features < min_features) {
    errors <- c(errors, paste("Insufficient number of features:", n_features, "<", min_features))
  }
  
  # Check for required assays
  if (!"Spatial" %in% names(seurat@assays) && !"RNA" %in% names(seurat@assays)) {
    errors <- c(errors, "No 'Spatial' or 'RNA' assay found in Seurat object")
  }
  
  # Validate spatial coordinates
  spatial_coords <- NULL
  if (require_spatial) {
    tryCatch({
      spatial_coords <- GetTissueCoordinates(seurat)
      
      if (is.null(spatial_coords) || nrow(spatial_coords) == 0) {
        errors <- c(errors, "No spatial coordinates found")
      } else {
        # Check coordinate validity
        if (ncol(spatial_coords) < 2) {
          errors <- c(errors, "Spatial coordinates must have at least 2 dimensions")
        }
        
        if (any(is.na(spatial_coords[, 1:2]))) {
          errors <- c(errors, "Spatial coordinates contain missing values")
        }
        
        # Check for degenerate coordinates (all points the same)
        if (nrow(spatial_coords) > 1) {
          coord_ranges <- apply(spatial_coords[, 1:2], 2, function(x) diff(range(x, na.rm = TRUE)))
          if (any(coord_ranges == 0)) {
            warnings <- c(warnings, "Some coordinate dimensions have no variation")
          }
        }
      }
    }, error = function(e) {
      errors <- c(errors, paste("Error accessing spatial coordinates:", e$message))
    })
  }
  
  # Check metadata structure
  if (is.null(seurat@meta.data) || nrow(seurat@meta.data) == 0) {
    errors <- c(errors, "No metadata found in Seurat object")
  } else {
    # Check sample identifiers
    if (!"orig.ident" %in% colnames(seurat@meta.data)) {
      warnings <- c(warnings, "No 'orig.ident' column found in metadata")
    } else {
      n_samples <- length(unique(seurat$orig.ident))
      if (n_samples == 0) {
        errors <- c(errors, "No valid sample identifiers found")
      }
    }
  }
  
  # Validate TLS scores if required
  if (require_tls_score) {
    if (!"TLS" %in% colnames(seurat@meta.data)) {
      errors <- c(errors, "No 'TLS' score column found in metadata - run preprocessing first")
    } else {
      tls_scores <- seurat@meta.data$TLS
      
      if (all(is.na(tls_scores))) {
        errors <- c(errors, "All TLS scores are missing")
      } else if (any(is.na(tls_scores))) {
        warnings <- c(warnings, paste("Some TLS scores are missing:", sum(is.na(tls_scores)), "out of", length(tls_scores)))
      }
      
      # Check score distribution
      if (length(unique(tls_scores[!is.na(tls_scores)])) == 1) {
        warnings <- c(warnings, "All TLS scores are identical - may indicate preprocessing issues")
      }
    }
  }
  
  # Data quality checks
  quality_summary <- list()
  
  # Expression data checks
  tryCatch({
    default_assay <- DefaultAssay(seurat)
    expr_data <- GetAssayData(seurat, slot = "counts", assay = default_assay)
    
    # Check for empty cells/features
    cells_with_counts <- colSums(expr_data) > 0
    features_with_counts <- rowSums(expr_data) > 0
    
    quality_summary$empty_cells <- sum(!cells_with_counts)
    quality_summary$empty_features <- sum(!features_with_counts)
    quality_summary$total_counts <- sum(expr_data)
    quality_summary$median_counts_per_cell <- median(colSums(expr_data))
    
    if (quality_summary$empty_cells > 0) {
      warnings <- c(warnings, paste(quality_summary$empty_cells, "cells have zero counts"))
    }
    
    if (quality_summary$empty_features > n_features * 0.5) {
      warnings <- c(warnings, "More than 50% of features have zero counts across all cells")
    }
    
  }, error = function(e) {
    warnings <- c(warnings, paste("Error accessing expression data:", e$message))
  })
  
  # Spatial data quality if available
  if (!is.null(spatial_coords) && nrow(spatial_coords) > 0) {
    quality_summary$spatial_range_x <- diff(range(spatial_coords[, 1], na.rm = TRUE))
    quality_summary$spatial_range_y <- diff(range(spatial_coords[, 2], na.rm = TRUE))
    
    # Check for suspicious spatial patterns
    if (nrow(spatial_coords) > 2) {
      # Calculate minimum distances
      min_dists <- apply(spatial_coords[, 1:2], 1, function(point) {
        dists <- sqrt(rowSums(sweep(spatial_coords[, 1:2], 2, point)^2))
        min(dists[dists > 0])
      })
      
      quality_summary$min_spatial_distance <- min(min_dists, na.rm = TRUE)
      quality_summary$median_spatial_distance <- median(min_dists, na.rm = TRUE)
      
      if (quality_summary$min_spatial_distance == 0) {
        warnings <- c(warnings, "Some spatial coordinates are identical")
      }
    }
  }
  
  # Overall validation result
  is_valid <- length(errors) == 0
  
  # Create summary
  summary_info <- list(
    n_cells = n_cells,
    n_features = n_features,
    n_samples = if("orig.ident" %in% colnames(seurat@meta.data)) length(unique(seurat$orig.ident)) else 1,
    has_spatial_coords = !is.null(spatial_coords) && nrow(spatial_coords) > 0,
    has_tls_scores = "TLS" %in% colnames(seurat@meta.data),
    quality_metrics = quality_summary
  )
  
  return(list(
    valid = is_valid,
    errors = errors,
    warnings = warnings,
    summary = summary_info
  ))
}

#' Validate parameters for TLS detection functions
#' 
#' @param exp_threshold Expression threshold value
#' @param min_spots Minimum spots parameter
#' @param max_distance Maximum distance parameter
#' @param distance_multiplier Distance multiplier parameter
#' @param alpha Statistical significance level
#' @return List with validation results
validateTLSParameters <- function(exp_threshold = 0.98,
                                  min_spots = 3,
                                  max_distance = NULL,
                                  distance_multiplier = 2,
                                  alpha = 0.05,
                                  use_statistical_threshold = TRUE) {
  
  errors <- character(0)
  warnings <- character(0)
  
  # Validate exp_threshold
  if (!is.numeric(exp_threshold) || length(exp_threshold) != 1) {
    errors <- c(errors, "exp_threshold must be a single numeric value")
  } else if (exp_threshold < 0 || exp_threshold > 1) {
    errors <- c(errors, "exp_threshold must be between 0 and 1")
  } else if (exp_threshold < 0.8) {
    warnings <- c(warnings, "exp_threshold < 0.8 may result in many false positives")
  }
  
  # Validate min_spots
  if (!is.numeric(min_spots) || length(min_spots) != 1) {
    errors <- c(errors, "min_spots must be a single numeric value")
  } else if (min_spots < 1 || min_spots != round(min_spots)) {
    errors <- c(errors, "min_spots must be a positive integer")
  } else if (min_spots > 20) {
    warnings <- c(warnings, "min_spots > 20 may miss smaller TLS structures")
  }
  
  # Validate max_distance
  if (!is.null(max_distance)) {
    if (!is.numeric(max_distance) || length(max_distance) != 1) {
      errors <- c(errors, "max_distance must be a single numeric value or NULL")
    } else if (max_distance <= 0) {
      errors <- c(errors, "max_distance must be positive")
    }
  }
  
  # Validate distance_multiplier
  if (!is.numeric(distance_multiplier) || length(distance_multiplier) != 1) {
    errors <- c(errors, "distance_multiplier must be a single numeric value")
  } else if (distance_multiplier <= 0) {
    errors <- c(errors, "distance_multiplier must be positive")
  } else if (distance_multiplier > 5) {
    warnings <- c(warnings, "distance_multiplier > 5 may create overly large neighborhoods")
  }
  
  # Validate alpha
  if (!is.numeric(alpha) || length(alpha) != 1) {
    errors <- c(errors, "alpha must be a single numeric value")
  } else if (alpha <= 0 || alpha >= 1) {
    errors <- c(errors, "alpha must be between 0 and 1")
  } else if (alpha > 0.2) {
    warnings <- c(warnings, "alpha > 0.2 may be too permissive for multiple testing")
  }
  
  # Validate use_statistical_threshold
  if (!is.logical(use_statistical_threshold) || length(use_statistical_threshold) != 1) {
    errors <- c(errors, "use_statistical_threshold must be a single logical value")
  }
  
  return(list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings
  ))
}

#' Print informative error messages with suggestions
#' 
#' @param validation_result Result from validateSeuratInput or validateTLSParameters
#' @param context Character string describing the context ("input" or "parameters")
#' @export
printValidationErrors <- function(validation_result, context = "input") {
  
  if (!validation_result$valid) {
    cat("Validation failed for", context, ":\n\n")
    
    if (length(validation_result$errors) > 0) {
      cat("ERRORS:\n")
      for (i in seq_along(validation_result$errors)) {
        cat(paste0("  ", i, ". ", validation_result$errors[i], "\n"))
      }
      cat("\n")
    }
    
    if (length(validation_result$warnings) > 0) {
      cat("WARNINGS:\n")
      for (i in seq_along(validation_result$warnings)) {
        cat(paste0("  ", i, ". ", validation_result$warnings[i], "\n"))
      }
      cat("\n")
    }
    
    # Provide suggestions based on common errors
    if (context == "input") {
      cat("SUGGESTIONS:\n")
      cat("  - Ensure your Seurat object was created from spatial transcriptomics data\n")
      cat("  - Run the preprocess() function before TLS detection\n")
      cat("  - Check that spatial coordinates are properly loaded\n")
      cat("  - Verify that TLS module scores have been calculated\n")
    } else if (context == "parameters") {
      cat("SUGGESTIONS:\n")
      cat("  - Use exp_threshold between 0.90-0.99 for most datasets\n")
      cat("  - Start with min_spots = 3-5 for initial exploration\n")
      cat("  - Consider use_statistical_threshold = TRUE for more robust results\n")
      cat("  - Set alpha = 0.05 for standard statistical significance\n")
    }
  } else {
    cat("Validation passed for", context, "\n")
    
    if (length(validation_result$warnings) > 0) {
      cat("\nWARNINGS:\n")
      for (i in seq_along(validation_result$warnings)) {
        cat(paste0("  ", i, ". ", validation_result$warnings[i], "\n"))
      }
    }
  }
}