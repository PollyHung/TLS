#' Preprocess Spatial Transcriptomics Data
#'
#' Reads, merges, and processes spatial transcriptomics data through a complete preprocessing pipeline
#' including quality control, normalization, integration, and dimension reduction.
#'
#' @param samples A character vector of sample names/directories containing spatial data.
#' @param path Base directory path where sample directories are located.
#' @param folder Optional subdirectory within each sample directory containing data files (e.g., "filtered_feature_bc_matrix").
#'
#' @return A processed \code{Seurat} object containing:
#' \itemize{
#'   \item Merged data from all samples with unique cell IDs.
#'   \item Quality control metrics (mitochondrial/ribosomal percentages).
#'   \item Cell cycle phase scores (S.Score, G2M.Score).
#'   \item Normalized and scaled expression data.
#'   \item PCA, Harmony integration, UMAP, and t-SNE embeddings.
#'   \item Nearest neighbor graph for downstream analysis.
#' }
#'
#' @details The pipeline performs these steps in sequence:
#' \enumerate{
#'   \item Read and merge 10X Genomics data from multiple samples.
#'   \item Calculate QC metrics (mitochondrial/ribosomal genes).
#'   \item Normalize data and identify variable features.
#'   \item Score cell cycle phases using package-provided S/G2M gene sets.
#'   \item Scale data with regression of technical variables.
#'   \item Dimensionality reduction (PCA).
#'   \item Batch correction using Harmony.
#'   \item Non-linear dimension reduction (UMAP/t-SNE).
#' }
#'
#' @section Required Data:
#' Uses internal package data:
#' \itemize{
#'   \item \code{s.genes}: S-phase marker genes.
#'   \item \code{g2m.genes}: G2/M-phase marker genes.
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
                       path){

  message("________________Start of Preprocessing________________✧*｡٩(ˊᗜˋ*)و✧*｡________________")
  message("Read in all samples")
  myList <- list()
  variableFeatures <- list()
  for (i in samples){
    message(paste0("Processing sample ", i))
    dir <- file.path(path, i)                                                   # Get directory
    obj <- Seurat::Load10X_Spatial(dir, slice = i)                              # Load Object

    Idents(obj) <- i                                                            # add identity
    obj$orig.ident <- i                                                         # add identity
    obj <- Seurat::RenameCells(obj, new.names = paste0(i, "_",Cells(obj)))      # add name
    obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE,
                       return.only.var.genes = FALSE, ,
                       variable.features.n = nrow(obj)) # perform SCT

    DefaultAssay(obj) <- "SCT"
    variableFeatures[[i]] <- VariableFeatures(obj)
    myList[[i]] <- obj
  }

  message("Merge seurat object")
  seurat <- Reduce(function(x, y) merge(x, y), myList)

  message("________________((꜆꜄꜆ ˙꒳˙)꜆꜄꜆ｵﾗｵﾗｵﾗｵﾗｲ________________")
  message("Add Quality Control Metrics")
  seurat[["mito_percent"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat[["ribo_percent"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^RP[SL]")
  seurat <- Seurat::CellCycleScoring(seurat,
                                     s.features = annotateTLS::s.genes,
                                     g2m.features = annotateTLS::g2m.genes,
                                     set.ident = FALSE)

  message("Run PCA, Find Neighbors, Find Clusters, Run UMAP, and Run TSNE")
  VariableFeatures(seurat) <- unique(unname(unlist(variableFeatures)))          # Manually Add the Variable Features Aggregated from Previous Attemps
  seurat <- RunPCA(seurat, verbose = FALSE)
  seurat <- RunHarmony(seurat, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony", lambda = 1)
  seurat <- FindNeighbors(seurat, dims = 1:30)
  seurat <- FindClusters(seurat, verbose = FALSE)
  seurat <- RunUMAP(seurat, dims = 1:30)
  seurat <- RunTSNE(seurat, dims = 1:30)

  message("Add Module Scores for key markers")
  seurat <- AddModuleScore(object = seurat,
                           features = list(annotateTLS::tls_50_genes,
                                           annotateTLS::b_cell_markers,
                                           annotateTLS::t_cell_markers,
                                           annotateTLS::follicular_dc),
                           name = c("TLS", "B.cell", "T.cell", "fDC"),
                           slot = "data",
                           group.by = "orig.ident")
  seurat@meta.data <- seurat@meta.data %>%
    dplyr::rename(TLS = TLS1, B.cell = B.cell2, T.cell = T.cell3,  fDC = fDC4)

  ## This is a separate Method
  # seurat@meta.data$TLS <- rowMeans(FetchData(seurat, annotateTLS::tls_50_genes, layer = "scale.data"))
  # seurat@meta.data$B.cell <- rowMeans(FetchData(seurat, annotateTLS::b_cell_markers, layer = "scale.data"))
  # seurat@meta.data$T.cell <- rowMeans(FetchData(seurat, annotateTLS::t_cell_markers, layer = "scale.data"))
  # seurat@meta.data$fDC <- rowMeans(FetchData(seurat, annotateTLS::follicular_dc, layer = "scale.data"))

  message("________________End of Preprocessing________________*ଘ(੭*ˊᵕˋ)੭* ੈ✩‧₊˚________________")
  return(seurat)
}



