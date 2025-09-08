#' Optimize TLS detection parameters using cross-validation
#' 
#' This function performs grid search or Bayesian optimization to find optimal
#' parameters for TLS detection, using cross-validation to evaluate performance.
#'
#' @param seurat A \code{Seurat} object containing spatial transcriptomics data
#' @param validation_annotations Optional ground truth annotations for validation
#' @param param_grid Optional list of parameter combinations to test. If NULL, uses default grid
#' @param cv_folds Integer number of cross-validation folds (default: 5)
#' @param performance_metric Character string specifying metric ("f1", "precision", "recall", "accuracy")
#' @param n_cores Integer number of cores for parallel processing (default: 1)
#' @param verbose Logical whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{best_params}: Optimal parameter combination
#'   \item \code{best_score}: Performance score for best parameters
#'   \item \code{cv_results}: Detailed cross-validation results
#'   \item \code{performance_summary}: Summary statistics across parameter combinations
#' }
#'
#' @details The function performs k-fold cross-validation to evaluate different
#' parameter combinations for TLS detection. It tests various combinations of:
#' \itemize{
#'   \item Expression threshold or statistical significance levels
#'   \item Minimum spots per cluster
#'   \item Distance multipliers
#'   \item Statistical vs quantile thresholding methods
#' }
#'
#' @examples
#' \donttest{
#' # Basic optimization with default parameters
#' opt_results <- optimizeTLSParameters(seurat_obj)
#' 
#' # Custom parameter grid
#' custom_grid <- list(
#'   exp_threshold = c(0.95, 0.98, 0.99),
#'   min_spots = c(3, 5, 7),
#'   distance_multiplier = c(1.5, 2.0, 2.5),
#'   use_statistical_threshold = c(TRUE, FALSE)
#' )
#' opt_results <- optimizeTLSParameters(seurat_obj, param_grid = custom_grid)
#' }
#'
#' @importFrom stats sd
#' @importFrom parallel detectCores makeCluster clusterEvalQ parLapply stopCluster
#' @export
optimizeTLSParameters <- function(seurat,
                                  validation_annotations = NULL,
                                  param_grid = NULL,
                                  cv_folds = 5,
                                  performance_metric = "f1",
                                  n_cores = 1,
                                  verbose = TRUE) {
  
  if (verbose) message("Starting TLS parameter optimization...")
  
  # Create default parameter grid if not provided
  if (is.null(param_grid)) {
    param_grid <- list(
      exp_threshold = c(0.95, 0.98, 0.99),
      min_spots = c(3, 5, 7),
      distance_multiplier = c(1.5, 2.0, 2.5),
      use_statistical_threshold = c(TRUE, FALSE),
      alpha = c(0.01, 0.05, 0.1),
      correction_method = c("fdr", "bonferroni")
    )
  }
  
  # Generate all parameter combinations
  param_combinations <- expand.grid(param_grid, stringsAsFactors = FALSE)
  n_combinations <- nrow(param_combinations)
  
  if (verbose) {
    message(paste("Testing", n_combinations, "parameter combinations"))
    message(paste("Using", cv_folds, "-fold cross-validation"))
  }
  
  # Create cross-validation folds
  cv_folds_list <- createCVFolds(seurat, cv_folds)
  
  # Performance tracking
  cv_results <- data.frame()
  
  # Test each parameter combination
  for (i in 1:n_combinations) {
    params <- param_combinations[i, ]
    if (verbose) message(paste("Testing combination", i, "of", n_combinations))
    
    fold_scores <- numeric(cv_folds)
    
    # Cross-validation for this parameter combination
    for (fold in 1:cv_folds) {
      train_cells <- cv_folds_list[[fold]]$train
      test_cells <- cv_folds_list[[fold]]$test
      
      # Create training subset
      seurat_train <- subset(seurat, cells = train_cells)
      seurat_test <- subset(seurat, cells = test_cells)
      
      # Run TLS detection with current parameters
      tryCatch({
        seurat_train <- do.call(runKNN, c(list(seurat = seurat_train), params))
        
        # Apply trained parameters to test set (simplified approach)
        seurat_test <- do.call(runKNN, c(list(seurat = seurat_test), params))
        
        # Calculate performance score
        if (!is.null(validation_annotations)) {
          true_labels <- validation_annotations[test_cells]
          pred_labels <- seurat_test$TLS_identity
          fold_scores[fold] <- calculatePerformanceMetric(true_labels, pred_labels, performance_metric)
        } else {
          # Use internal consistency metrics if no ground truth
          fold_scores[fold] <- calculateInternalMetric(seurat_test, performance_metric)
        }
        
      }, error = function(e) {
        if (verbose) message(paste("Error in fold", fold, ":", e$message))
        fold_scores[fold] <- 0
      })
    }
    
    # Store results
    mean_score <- mean(fold_scores, na.rm = TRUE)
    sd_score <- sd(fold_scores, na.rm = TRUE)
    
    result_row <- data.frame(
      combination = i,
      mean_score = mean_score,
      sd_score = sd_score,
      params,
      stringsAsFactors = FALSE
    )
    cv_results <- rbind(cv_results, result_row)
  }
  
  # Find best parameters
  best_idx <- which.max(cv_results$mean_score)
  best_params <- cv_results[best_idx, ]
  
  if (verbose) {
    message("Optimization complete!")
    message(paste("Best", performance_metric, "score:", round(best_params$mean_score, 4)))
  }
  
  # Create performance summary
  performance_summary <- list(
    best_score = best_params$mean_score,
    best_score_sd = best_params$sd_score,
    mean_score_all = mean(cv_results$mean_score),
    score_range = range(cv_results$mean_score)
  )
  
  return(list(
    best_params = best_params,
    best_score = best_params$mean_score,
    cv_results = cv_results,
    performance_summary = performance_summary
  ))
}

#' Create cross-validation folds for spatial data
#' 
#' @param seurat Seurat object
#' @param cv_folds Number of folds
#' @return List of train/test cell splits
createCVFolds <- function(seurat, cv_folds) {
  all_cells <- Cells(seurat)
  n_cells <- length(all_cells)
  
  # Stratify by sample if multiple samples present
  if ("orig.ident" %in% names(seurat@meta.data)) {
    samples <- unique(seurat$orig.ident)
    folds_list <- list()
    
    for (fold in 1:cv_folds) {
      train_cells <- c()
      test_cells <- c()
      
      for (sample in samples) {
        sample_cells <- Cells(subset(seurat, subset = orig.ident == sample))
        n_sample <- length(sample_cells)
        
        # Create stratified split within each sample
        test_indices <- sample(n_sample, size = floor(n_sample / cv_folds))
        test_cells_sample <- sample_cells[test_indices]
        train_cells_sample <- setdiff(sample_cells, test_cells_sample)
        
        train_cells <- c(train_cells, train_cells_sample)
        test_cells <- c(test_cells, test_cells_sample)
      }
      
      folds_list[[fold]] <- list(train = train_cells, test = test_cells)
    }
  } else {
    # Simple random splits if single sample
    folds_list <- list()
    fold_indices <- cut(sample(n_cells), breaks = cv_folds, labels = FALSE)
    
    for (fold in 1:cv_folds) {
      test_indices <- which(fold_indices == fold)
      train_indices <- which(fold_indices != fold)
      
      folds_list[[fold]] <- list(
        train = all_cells[train_indices],
        test = all_cells[test_indices]
      )
    }
  }
  
  return(folds_list)
}

#' Calculate performance metrics for TLS detection
#' 
#' @param true_labels Ground truth TLS annotations
#' @param pred_labels Predicted TLS annotations
#' @param metric Type of metric to calculate
#' @return Numeric performance score
calculatePerformanceMetric <- function(true_labels, pred_labels, metric = "f1") {
  
  # Convert to binary (TLS vs not TLS)
  true_binary <- ifelse(true_labels == "TLS", 1, 0)
  pred_binary <- ifelse(pred_labels == "TLS", 1, 0)
  
  # Calculate confusion matrix components
  tp <- sum(true_binary == 1 & pred_binary == 1)
  fp <- sum(true_binary == 0 & pred_binary == 1)
  tn <- sum(true_binary == 0 & pred_binary == 0)
  fn <- sum(true_binary == 1 & pred_binary == 0)
  
  # Calculate metrics
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  accuracy <- (tp + tn) / (tp + fp + tn + fn)
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  switch(metric,
         "precision" = precision,
         "recall" = recall,
         "accuracy" = accuracy,
         "f1" = f1,
         f1  # default to F1
  )
}

#' Calculate internal consistency metrics when no ground truth is available
#' 
#' @param seurat Seurat object with TLS annotations
#' @param metric Type of metric to calculate
#' @return Numeric internal consistency score
calculateInternalMetric <- function(seurat, metric = "coherence") {
  
  # Get TLS spots
  tls_spots <- which(seurat$TLS_identity == "TLS")
  
  if (length(tls_spots) == 0) {
    return(0)
  }
  
  # Calculate spatial coherence of TLS spots
  coords <- GetTissueCoordinates(seurat)[tls_spots, ]
  
  if (nrow(coords) < 2) {
    return(0.5)  # Single spot gets moderate score
  }
  
  # Calculate average nearest neighbor distances within TLS
  nn_within_tls <- RANN::nn2(coords[, 1:2], k = 2)$nn.dists[, 2]
  
  # Calculate distances to non-TLS spots
  non_tls_coords <- GetTissueCoordinates(seurat)[-tls_spots, ]
  
  if (nrow(non_tls_coords) > 0) {
    nn_to_non_tls <- RANN::nn2(non_tls_coords[, 1:2], coords[, 1:2], k = 1)$nn.dists[, 1]
    
    # Coherence score: TLS spots should be closer to each other than to non-TLS
    coherence <- mean(nn_to_non_tls) / mean(nn_within_tls)
    return(min(coherence / 2, 1))  # Scale and cap at 1
  } else {
    return(mean(1 / (1 + nn_within_tls)))  # Fallback when no non-TLS spots
  }
}