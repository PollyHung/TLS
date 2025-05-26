# annotateTLS: An R Package for Identifying Tumor-associated Lymphoid Structures in Spatial Transcriptomics Data

## Overview

The **annotateTLS** package provides a streamlined workflow for processing spatial transcriptomics data and identifying Tumor-associated Lymphoid Structures (TLS) using adaptive spatial neighborhood analysis. Designed for 10X Genomics Visium data, this package integrates preprocessing, quality control, batch correction, and TLS detection into a reproducible pipeline.

## Features

- **Complete Preprocessing**: Merge multiple samples, calculate QC metrics, perform normalization, and batch correction with Harmony.
- **Adaptive TLS Detection**: Identify TLS regions using K-nearest neighbors with sample-specific distance thresholds.
- **Visualization Tools**: Generate publication-ready plots for QC, feature expression, and spatial TLS distribution.

## Installation

```r
# Install from GitHub
devtools::install_github("PollyHung/annotateTLS")

# Install dependencies
install.packages(c("Seurat", "harmony", "igraph", "RANN", "scCustomize"))
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

## Basic Workflow

### 1. Preprocessing
```r
library(annotateTLS)
library(Seurat)
library(scCustomize)

# Set data path
path <- "~/Desktop/TLS/data/raw-spatial/E-MTAB-13084"
samples <- list.dirs(path, recursive = FALSE, full.names = FALSE)

# Run preprocessing pipeline
seurat <- preprocess(samples = samples, path = path, folder = "spatial")
```

### 2. Quality Control & Visualization
```r
# Generate QC plots
VlnPlot(seurat, features = c("nFeature_SCT", "nCount_SCT"), 
        group.by = "orig.ident", pt.size = 0)
ggsave("InterSampleDifference.png", width = 8, height = 4)

# Visualization
DimPlot_scCustom(seurat, group.by = "orig.ident", colors_use = "Set2")
FeaturePlot_scCustom(seurat, features = c("CD19", "CD3D", "CXCL13"))
```

### 3. TLS Identification
```r
# Detect TLS with adaptive KNN
seurat <- runKNN(seurat = seurat, exp_threshold = 0.96, min_spots = 5)

# Visualize results
SpatialDimPlot(seurat, group.by = "TLS_identity", 
               pt.size.factor = 2, ncol = 4, crop = FALSE)
```

## Expected Outputs

Here's the revised **Expected Outputs** section incorporating your request:

---

## Expected Outputs

### Core Output  
The pipeline returns **a fully integrated Seurat object** containing:  

| Component | Description |  
|-----------|-------------|  
| `TLS_identity` | Binary classification ("TLS"/"not TLS") in metadata from KNN analysis |  
| `TLS_score` | Continuous TLS probability score (optional, if implemented) |  
| Spatial coordinates | Original spot/cell locations in `images` slot |  
| Batch-corrected embeddings | Harmony/PCA/UMAP embeddings for visualization |  
| QC metrics | `nFeature_SCT`, `percent.mt`, `Phase`, etc. |  

### Visualization Outputs  
Automatically generated plots showing:  

| Plot | Description | Visualization Source |  
|------|-------------|----------------------|  
| <img src="examples/SpatialFeaturePlot-TLS-identity.png" width="200"> | Spatial TLS classification | `TLS_identity` metadata |  
| <img src="examples/FeaturePlot.png" width="200"> | Marker expression gradients | `TLS_score` + gene expression |  
| <img src="examples/InterSampleDifference.png" width="200"> | QC metrics across samples | Preprocessing metrics |  


## Advanced Configuration

### Key Parameters for `runKNN()`
- `exp_threshold`: Quantile cutoff for high-expression spots (default: 0.98)
- `min_spots`: Minimum cluster size for TLS validation (default: 3)
- `distance_multiplier`: Scaling factor for adaptive distance threshold (default: 2)

## Citation

If using this package in your research, please cite:
```
Hung et al. annotateTLS: Spatial TLS detection in tumor microenvironments. 2024.
```

## License
!!remember to add license 
