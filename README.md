# annotateTLS

**Detection and Annotation of Tertiary Lymphoid Structures in Spatial Transcriptomics Data**

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

The **annotateTLS** package provides a comprehensive, statistically rigorous framework for detecting and annotating tertiary lymphoid structures (TLS) in spatial transcriptomics data. Originally designed for 10X Genomics Visium data, the package now features enhanced statistical methods, confidence scoring, and automated parameter optimization.

### Key Features

- 🔬 **Statistical Rigor**: Replace arbitrary quantile thresholds with FDR-corrected statistical significance testing
- 🎯 **Confidence Scoring**: Multi-factor confidence assessment for each TLS prediction
- ⚙️ **Parameter Optimization**: Automated cross-validation-based parameter tuning
- ✅ **Robust Validation**: Comprehensive input validation with informative error messages  
- 🧪 **Comprehensive Testing**: Extensive test suite ensuring reliability and reproducibility
- 📊 **Enhanced Output**: Rich metadata including confidence scores and spatial coherence measures

## What's New in v0.2.0

### 🆕 Major Enhancements

#### Statistical Threshold Calculation
- **New Function**: `calculateStatisticalThreshold()`
- Replaces quantile-based thresholding with statistical significance testing
- Implements FDR correction using Benjamini-Hochberg procedure
- Provides p-values and adjusted p-values for transparency

#### Parameter Optimization Framework  
- **New Function**: `optimizeTLSParameters()`
- Automated grid search with k-fold cross-validation
- Performance metrics: precision, recall, F1-score, accuracy
- Works with or without ground truth annotations
- Reduces arbitrary parameter selection

#### Enhanced Confidence Scoring
- **Enhanced**: `runKNN()` now outputs confidence measures
- **New Columns**:
  - `TLS_confidence`: Probabilistic confidence score (0-1)
  - `TLS_cluster_id`: Unique cluster identifier
  - `spatial_coherence`: Spatial clustering strength measure
- Multi-factor scoring: expression strength + cluster size + spatial coherence

#### Robust Input Validation
- **New Functions**: `validateSeuratInput()`, `validateTLSParameters()`, `printValidationErrors()`
- Comprehensive data quality checks
- Informative error messages with actionable suggestions
- Prevents common user errors and data format issues

#### Comprehensive Testing Framework
- Unit tests for all major functions using `testthat`
- Property-based testing for spatial algorithms
- Edge case handling validation  
- Performance benchmarking capabilities

## Installation

### From Local Source
```r
# Install from local directory
devtools::install_local("path/to/annotateTLS")
```

### Dependencies
The package automatically installs required dependencies. For optional advanced features:

```r
# Required dependencies (automatically installed)
install.packages(c("Seurat", "RANN", "igraph", "harmony", "magrittr", "dplyr"))

# Optional dependencies for extended functionality  
install.packages(c("RCTD", "SPOTlight", "spatstat", "FNN", "cluster"))

# Testing dependencies (for developers)
install.packages("testthat")
```

## Data Preparation

Organize your spatial data with this directory structure:
```
raw-spatial/
└── E-MTAB-13084/
    ├── sample1/
    │   └── spatial/          # Contains filtered_feature_bc_matrix
    └── sample2/
        └── spatial/
```

## Quick Start Guide

### Basic Workflow (Recommended for New Users)

```r
library(annotateTLS)
library(Seurat)

# 1. Preprocess spatial transcriptomics data
seurat_obj <- preprocess(
  samples = c("sample1", "sample2"),
  path = "path/to/data"
)

# 2. Run TLS detection with statistical thresholding (recommended)
seurat_obj <- runKNN(
  seurat_obj,
  use_statistical_threshold = TRUE,  # Use statistical instead of quantile thresholds
  alpha = 0.05,                     # Significance level
  correction_method = "fdr",        # Multiple testing correction
  min_spots = 3                     # Minimum cluster size
)

# 3. Examine results
head(seurat_obj@meta.data[, c("TLS_identity", "TLS_confidence", "TLS_cluster_id", "spatial_coherence")])
```

### Advanced Workflow with Parameter Optimization

```r
# 1. Optimize parameters for your specific dataset
opt_results <- optimizeTLSParameters(
  seurat_obj, 
  cv_folds = 5,                    # Cross-validation folds
  performance_metric = "f1",       # Optimization metric
  n_cores = 2                      # Parallel processing
)

# 2. View optimization results
print(opt_results$best_params)
cat("Best F1 score:", round(opt_results$best_score, 3), "\n")

# 3. Apply optimized parameters
seurat_obj <- runKNN(
  seurat_obj, 
  exp_threshold = opt_results$best_params$exp_threshold,
  min_spots = opt_results$best_params$min_spots,
  distance_multiplier = opt_results$best_params$distance_multiplier,
  use_statistical_threshold = opt_results$best_params$use_statistical_threshold,
  alpha = opt_results$best_params$alpha
)
```

### Quality Control and Analysis

```r
# 1. Examine confidence score distribution
hist(seurat_obj$TLS_confidence, 
     main = "TLS Confidence Scores", 
     xlab = "Confidence Score", 
     breaks = 20)

# 2. Filter high-confidence TLS predictions
high_conf_threshold <- 0.7
high_conf_tls <- seurat_obj@meta.data$TLS_confidence > high_conf_threshold

# 3. Summary statistics
cat("Total TLS spots:", sum(seurat_obj$TLS_identity == "TLS"), "\n")
cat("High-confidence TLS spots:", sum(high_conf_tls & seurat_obj$TLS_identity == "TLS"), "\n")

# 4. Cluster analysis
tls_clusters <- table(seurat_obj$TLS_cluster_id[seurat_obj$TLS_identity == "TLS"])
cat("Number of TLS clusters:", length(tls_clusters[tls_clusters > 0]), "\n")
```

### Visualization

```r
# Spatial visualization of TLS identity
SpatialDimPlot(seurat_obj, 
               group.by = "TLS_identity", 
               pt.size.factor = 2,
               cols = c("gray80", "red"))

# Confidence score spatial plot  
SpatialFeaturePlot(seurat_obj, 
                   features = "TLS_confidence",
                   pt.size.factor = 2)

# Spatial coherence visualization
SpatialFeaturePlot(seurat_obj, 
                   features = "spatial_coherence",
                   pt.size.factor = 2)
```

## Output Structure

### Enhanced Metadata Columns

After running the enhanced `runKNN()` function, your Seurat object contains rich metadata:

| Column | Type | Description | Range/Values |
|--------|------|-------------|--------------|
| `TLS_identity` | Factor | Binary TLS classification | "TLS", "not TLS" |
| `TLS_confidence` | Numeric | Multi-factor confidence score | 0-1 (higher = more confident) |
| `TLS_cluster_id` | Integer | Unique identifier for each TLS cluster | 0 = not TLS, >0 = cluster ID |
| `spatial_coherence` | Numeric | Spatial clustering coherence measure | 0-1 (higher = more spatially coherent) |
| `TLS` | Numeric | Original TLS module score from preprocessing | Continuous |

### Statistical Output (when using statistical thresholding)

The `calculateStatisticalThreshold()` function provides:

```r
# Example output structure
$threshold
[1] 2.341

$method
[1] "statistical_fdr"

$n_significant
[1] 127

$alpha
[1] 0.05
```

### Parameter Optimization Results

The `optimizeTLSParameters()` function returns:

```r
# Example optimization output
$best_params
  combination mean_score exp_threshold min_spots distance_multiplier use_statistical_threshold alpha correction_method
1          15      0.763          0.98         3                 2.0                     TRUE  0.05               fdr

$best_score
[1] 0.763

$performance_summary
$best_score
[1] 0.763

$mean_score_all  
[1] 0.524

$score_range
[1] 0.234 0.763
```  


## Advanced Configuration

### Key Parameters for Enhanced `runKNN()`

| Parameter | Default | Description | Recommendations |
|-----------|---------|-------------|-----------------|
| `use_statistical_threshold` | `TRUE` | Use statistical vs quantile thresholding | **TRUE** for robust results |
| `alpha` | `0.05` | Significance level for statistical testing | 0.01-0.1 depending on stringency |
| `correction_method` | `"fdr"` | Multiple testing correction | `"fdr"` or `"bonferroni"` |
| `exp_threshold` | `0.98` | Quantile threshold (when statistical=FALSE) | 0.90-0.99 for different sensitivities |
| `min_spots` | `3` | Minimum cluster size | 3-7 depending on resolution |
| `distance_multiplier` | `2` | Spatial distance scaling factor | 1.5-3.0 for different neighborhood sizes |

### Validation Functions

```r
# Validate your input data before analysis
validation_result <- validateSeuratInput(seurat_obj)
if (!validation_result$valid) {
  printValidationErrors(validation_result, "input")
}

# Check parameter validity
param_validation <- validateTLSParameters(
  exp_threshold = 0.98,
  alpha = 0.05,
  min_spots = 3
)
```

### Performance Tuning

```r
# For large datasets, consider parallel processing
opt_results <- optimizeTLSParameters(
  seurat_obj,
  n_cores = parallel::detectCores() - 1,  # Use available cores
  cv_folds = 3  # Reduce for faster computation
)
```

## Dependencies and System Requirements

### Required
- **R** ≥ 4.0
- **Seurat** ≥ 4.0.0
- **RANN**, **igraph**, **stats**, **parallel**
- **magrittr**, **dplyr**, **harmony**

### Suggested (for extended functionality)
- **RCTD**, **SPOTlight** (cell type deconvolution)
- **spatstat**, **FNN**, **cluster** (advanced spatial analysis)
- **testthat** ≥ 3.0.0 (testing framework)

### System Requirements
- **Memory**: Minimum 8GB RAM (16GB+ recommended for large datasets)
- **Storage**: ~100MB for package + data-dependent space
- **OS**: Windows, macOS, or Linux with R support

## Troubleshooting

### Common Issues

**Error: "Input validation failed"**
- Ensure Seurat object has spatial coordinates
- Run `preprocess()` function first to calculate TLS module scores
- Check data format with `validateSeuratInput()`

**Warning: "No statistically significant spots found"**
- Try less stringent alpha (e.g., 0.1)
- Consider `use_statistical_threshold = FALSE`
- Check TLS module score distribution

**Low confidence scores**
- Verify preprocessing quality (TLS module scores)
- Consider parameter optimization with `optimizeTLSParameters()`
- Check spatial coordinate quality

## Citation

If you use this package in your research, please cite:

```bibtex
@misc{annotateTLS2024,
  title={annotateTLS: Statistical Detection and Annotation of Tertiary Lymphoid Structures in Spatial Transcriptomics Data},
  author={Hung, Polly and Contributors},
  year={2024},
  note={R package version 0.2.0}
}
```

## License

MIT License - see LICENSE file for details

## Support and Contributing

- **Issues**: Report bugs and feature requests via GitHub issues
- **Questions**: Contact the maintainer or open a discussion
- **Contributing**: Pull requests welcome - see CONTRIBUTING.md for guidelines

---

**Developed with ❤️ for the spatial transcriptomics community** 
