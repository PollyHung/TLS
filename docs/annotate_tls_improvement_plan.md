# annotateTLS Package Improvement Plan

## Phase 1: Statistical Rigor and Validation Framework (Weeks 1-3)

### 1.1 Implement Statistical Testing and Multiple Correction

**Files to modify:** `R/runKNN.R`

**New function:** `calculateStatisticalThreshold()`
```r
#' Calculate statistically rigorous expression thresholds
#' 
#' @param expression_matrix Matrix of expression values
#' @param alpha Significance level (default: 0.05)
#' @param method Multiple testing correction method ("fdr", "bonferroni")
#' @return List with threshold value and adjusted p-values
calculateStatisticalThreshold <- function(expression_matrix, alpha = 0.05, method = "fdr") {
  # Implement statistical testing for each gene
  # Use appropriate distribution assumptions (negative binomial for counts)
  # Apply multiple testing correction
  # Return statistically justified thresholds
}
```

**Tasks:**
- [ ] Replace quantile-based thresholding with statistical significance testing
- [ ] Implement FDR correction using Benjamini-Hochberg procedure
- [ ] Add parameter validation for statistical parameters
- [ ] Document statistical assumptions and limitations

### 1.2 Add Parameter Optimization Framework

**New file:** `R/optimizeParameters.R`

```r
#' Optimize TLS detection parameters using cross-validation
#' 
#' @param seurat Seurat object
#' @param validation_annotations Ground truth annotations (if available)
#' @param param_grid List of parameter combinations to test
#' @return Optimized parameters with performance metrics
optimizeTLSParameters <- function(seurat, validation_annotations = NULL, param_grid = NULL) {
  # Implement grid search or Bayesian optimization
  # Use cross-validation for parameter selection
  # Return optimal parameters with confidence intervals
}
```

**Tasks:**
- [ ] Create default parameter grids for optimization
- [ ] Implement k-fold cross-validation
- [ ] Add performance metrics calculation (precision, recall, F1)
- [ ] Generate parameter sensitivity analysis plots

### 1.3 Implement Confidence Scoring

**Modify:** `R/runKNN.R` - update `runKNN()` function

**New output columns:**
- `TLS_confidence`: Probabilistic confidence score (0-1)
- `TLS_maturity_score`: Continuous measure of TLS development stage
- `spatial_coherence`: Measure of spatial clustering strength

**Tasks:**
- [ ] Calculate confidence based on multiple factors (expression strength, spatial coherence, cluster size)
- [ ] Implement TLS maturity scoring based on marker gene combinations
- [ ] Add uncertainty quantification for edge cases

## Phase 2: Biological Validation and Enhancement (Weeks 4-6)

### 2.1 Cell Type Deconvolution Integration

**New file:** `R/cellTypeDeconvolution.R`

```r
#' Integrate cell type deconvolution for TLS detection
#' 
#' @param seurat Seurat object
#' @param reference_data Single-cell reference for deconvolution
#' @param method Deconvolution method ("RCTD", "SPOTlight", "Cell2location")
#' @return Seurat object with cell type proportions
addCellTypeDeconvolution <- function(seurat, reference_data, method = "RCTD") {
  # Implement cell type deconvolution
  # Focus on TLS-relevant cell types (B cells, T cells, FDCs)
  # Calculate cell type co-occurrence scores
}
```

**Tasks:**
- [ ] Integrate with RCTD or SPOTlight packages
- [ ] Define TLS-relevant cell type combinations
- [ ] Add cell type spatial organization analysis
- [ ] Validate against known TLS cell type patterns

### 2.2 Multi-scale Spatial Analysis

**New file:** `R/multiScaleSpatial.R`

```r
#' Analyze TLS features at multiple spatial scales
#' 
#' @param seurat Seurat object with coordinates
#' @param scales Vector of distance scales to analyze
#' @return Multi-scale spatial features
calculateMultiScaleFeatures <- function(seurat, scales = c(1, 2, 3, 5)) {
  # Calculate spatial features at different neighborhood sizes
  # Implement spatial autocorrelation measures
  # Add morphological shape features
}
```

**Tasks:**
- [ ] Implement Moran's I for spatial autocorrelation
- [ ] Calculate shape descriptors (compactness, elongation)
- [ ] Add multi-scale clustering analysis
- [ ] Integrate with existing KNN approach

### 2.3 Enhanced Gene Signature Analysis

**New file:** `R/geneSignatures.R`

```r
#' Comprehensive TLS gene signature analysis
#' 
#' @param seurat Seurat object
#' @param custom_signatures Optional custom gene signatures
#' @return Enhanced signature scores and classifications
calculateTLSSignatures <- function(seurat, custom_signatures = NULL) {
  # Implement multiple TLS signature methods
  # Add maturation stage-specific signatures
  # Calculate signature reliability scores
}
```

**Tasks:**
- [ ] Curate comprehensive TLS signature database
- [ ] Implement stage-specific signatures (early, intermediate, mature TLS)
- [ ] Add signature correlation analysis
- [ ] Validate against literature-defined signatures

## Phase 3: Software Engineering and Testing (Weeks 7-8)

### 3.1 Comprehensive Testing Suite

**New directory:** `tests/testthat/`

**Test files to create:**
- `test-runKNN.R`: Test core TLS detection algorithm
- `test-preprocess.R`: Test preprocessing pipeline
- `test-parameters.R`: Test parameter optimization
- `test-validation.R`: Test validation functions

**Tasks:**
- [ ] Create synthetic datasets with known TLS patterns
- [ ] Test edge cases (empty datasets, single cells, extreme parameters)
- [ ] Implement property-based testing for spatial algorithms
- [ ] Add performance benchmarking tests

### 3.2 Input Validation and Error Handling

**New file:** `R/validation.R`

```r
#' Validate input data for TLS analysis
#' 
#' @param seurat Seurat object to validate
#' @return Boolean with detailed error messages
validateSeuratInput <- function(seurat) {
  # Check required slots and metadata
  # Validate spatial coordinates
  # Check for missing or problematic data
}
```

**Tasks:**
- [ ] Add comprehensive input validation
- [ ] Implement informative error messages
- [ ] Add data quality assessments
- [ ] Create data format conversion utilities

### 3.3 Documentation and Vignettes

**New files:**
- `vignettes/tls-detection-workflow.Rmd`
- `vignettes/parameter-optimization.Rmd`
- `vignettes/biological-interpretation.Rmd`

**Tasks:**
- [ ] Create comprehensive workflow vignette
- [ ] Add biological interpretation guidelines
- [ ] Document parameter selection strategies
- [ ] Include troubleshooting guide

## Phase 4: Benchmarking and Validation (Weeks 9-10)

### 4.1 Comparative Analysis Framework

**New file:** `R/benchmark.R`

```r
#' Benchmark against existing TLS detection methods
#' 
#' @param seurat Seurat object
#' @param methods Vector of methods to compare
#' @param ground_truth Optional ground truth annotations
#' @return Comparative performance metrics
benchmarkTLSMethods <- function(seurat, methods = c("quantile", "statistical", "signature"), ground_truth = NULL) {
  # Compare multiple detection approaches
  # Calculate comprehensive performance metrics
  # Generate comparison visualizations
}
```

**Tasks:**
- [ ] Implement alternative TLS detection methods for comparison
- [ ] Create standardized evaluation metrics
- [ ] Add cross-dataset validation framework
- [ ] Generate performance comparison reports

### 4.2 External Validation Pipeline

**New file:** `R/externalValidation.R`

**Tasks:**
- [ ] Create framework for pathologist annotation comparison
- [ ] Implement inter-rater reliability analysis
- [ ] Add visualization tools for validation results
- [ ] Create validation report generation

## Phase 5: Advanced Features and Integration (Weeks 11-12)

### 5.1 Machine Learning Integration

**New file:** `R/mlIntegration.R`

```r
#' Prepare features for machine learning models
#' 
#' @param seurat Seurat object with TLS annotations
#' @return Feature matrix suitable for ML training
prepareTLSFeatures <- function(seurat) {
  # Extract comprehensive feature set
  # Include spatial, expression, and morphological features
  # Prepare data for external ML frameworks
}
```

### 5.2 Visualization Enhancements

**New file:** `R/visualization.R`

**Tasks:**
- [ ] Create TLS-specific plotting functions
- [ ] Add interactive visualization options
- [ ] Implement confidence visualization
- [ ] Create summary report generation

## Implementation Guidelines for Claude Code

### Setup and Environment
```bash
# Initialize project structure
/init
# This creates a comprehensive R package development environment

# Add to CLAUDE.md:
# This is an R package for TLS detection in spatial transcriptomics data
# Focus on biological accuracy, statistical rigor, and reproducibility
# Priority: biological validation over algorithmic complexity
```

### Development Workflow
1. **Start with Phase 1, Task 1.1** - Statistical rigor is the foundation
2. **Use test-driven development** - Write tests before implementing features
3. **Maintain backward compatibility** - Don't break existing user workflows
4. **Document biological rationale** - Every parameter needs biological justification
5. **Regular validation** - Test against synthetic data at each step

### Code Quality Standards
- **Function length**: Maximum 50 lines per function
- **Documentation**: Every function needs examples with biological interpretation
- **Error handling**: Informative messages that guide users to solutions
- **Performance**: Profile memory usage for large datasets
- **Reproducibility**: Set and document random seeds

### Key Dependencies to Add
```r
# Add to DESCRIPTION
Imports:
    RCTD,
    SPOTlight,
    spatstat,
    FNN,
    cluster,
    testthat,
    knitr,
    rmarkdown
```

### Priority Order for Implementation
1. **Statistical threshold calculation** (highest impact on accuracy)
2. **Confidence scoring** (addresses computer science team's label noise concerns)
3. **Cell type deconvolution** (biological foundation)
4. **Comprehensive testing** (prevents regressions)
5. **Parameter optimization** (reduces arbitrary choices)
6. **Multi-scale analysis** (captures TLS complexity)

### Success Metrics
- [ ] Reproducible results across parameter ranges
- [ ] Correlation with pathologist annotations >0.8
- [ ] Cross-dataset F1 scores >0.75
- [ ] Package passes R CMD check with zero warnings
- [ ] Comprehensive test coverage >90%

This plan addresses the fundamental issues identified in the review while providing a clear roadmap for enhancing the package's biological validity and statistical rigor.