#' Preprocess Spatial Transcriptomics Data
#'
#' Reads, merges, and processes spatial transcriptomics data through a complete preprocessing pipeline
#' including quality control, normalization, integration, and dimension reduction.
#'
#' @param samples A character vector of sample names/directories containing spatial data
#' @param path Base directory path where sample directories are located
#' @param folder Optional subdirectory within each sample directory containing data files (e.g., "filtered_feature_bc_matrix")
#'
#' @return A processed \code{Seurat} object containing:
#' \itemize{
#'   \item Merged data from all samples with unique cell IDs
#'   \item Quality control metrics (mitochondrial/ribosomal percentages)
#'   \item Cell cycle phase scores (S.Score, G2M.Score)
#'   \item Normalized and scaled expression data
#'   \item PCA, Harmony integration, UMAP, and t-SNE embeddings
#'   \item Nearest neighbor graph for downstream analysis
#' }
#'
#' @details The pipeline performs these steps in sequence:
#' \enumerate{
#'   \item Read and merge 10X Genomics data from multiple samples
#'   \item Calculate QC metrics (mitochondrial/ribosomal genes)
#'   \item Normalize data and identify variable features
#'   \item Score cell cycle phases using package-provided S/G2M gene sets
#'   \item Scale data with regression of technical variables
#'   \item Dimensionality reduction (PCA)
#'   \item Batch correction using Harmony
#'   \item Non-linear dimension reduction (UMAP/t-SNE)
#' }
#'
#' @section Required Data:
#' Uses internal package data:
#' \itemize{
#'   \item \code{s.genes}: S-phase marker genes
#'   \item \code{g2m.genes}: G2/M-phase marker genes
#' }
#' See \code{?s.genes} and \code{?g2m.genes} for details.
#'
#' @examples
#' \donttest{
#' samples <- c("sample1", "sample2")
#' seurat_obj <- preprocess(
#'   samples = samples,
#'   path = "~/Desktop/TLS/data/raw-spatial/E-MTAB-13084",
#'   folder = "spatial"
#' )
#' }
#'
#' @importFrom Seurat Read10X CreateSeuratObject RenameCells PercentageFeatureSet
#' @importFrom Seurat CellCycleScoring NormalizeData FindVariableFeatures ScaleData
#' @importFrom Seurat RunPCA FindNeighbors RunUMAP RunTSNE
#' @importFrom SeuratObject JoinLayers
#' @importFrom harmony RunHarmony
#' @importFrom magrittr %>%
#' @export

preprocess <- function(samples,
                       path,
                       folder = NULL){

  message("Read in single-cell spots for all samples")
  myList <- list()
  for (i in samples){
    dir <- file.path(path, i, folder) # get directory
    data <- Seurat::Read10X(data.dir = dir) # load file
    obj <- Seurat::CreateSeuratObject(counts = data, project = i) # create seurat
    obj <- Seurat::RenameCells(obj, new.names = paste0(i, "_",Cells(obj))) # add name
    myList[[dir]] <- obj
  }

  message("Merge seurat object")
  seurat <- Reduce(function(x, y) merge(x, y), myList)
  seurat <- SeuratObject::JoinLayers(seurat)

  message("Add Quality Control Metrics")
  seurat[["mito_percent"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^MT-", assay = "RNA")
  seurat[["ribo_percent"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^RP[SL]", assay = "RNA")

  message("Normalisation and find variable features")
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 5000)

  message("Add Cell Cycle Scoring")
  seurat <- Seurat::CellCycleScoring(seurat,
                                     s.features = annotateTLS::s.genes,
                                     g2m.features = annotateTLS::g2m.genes,
                                     set.ident = FALSE)

  message("Add TLS Score")
  t_cell_markers <- c("CD3D", "CD3E", "CD8A", "CD8B", "CD4", "GZMA", "GZMB",
                      "")
  seurat <- AddModuleScore(object = seurat,
                           features = list(annotateTLS::tls_50_genes,
                                           ))


  return(seurat)
}







