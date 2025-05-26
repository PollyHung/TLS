#' Preprocess Spatial Transcriptomics Data
#'
#' Reads, merges, and adds quality control metrics to spatial transcriptomics data from multiple samples.
#'
#' @param samples A character vector of sample names/directories containing spatial data.
#' @param path Base directory path where sample directories are located.
#' @param folder Optional subdirectory within each sample directory containing the data (e.g., "spatial").
#'               Defaults to NULL.
#'
#' @return A merged \code{Seurat} object with the following additions:
#' \itemize{
#'   \item Cell-level metadata with mitochondrial/ribosomal percentages
#'   \item Cell cycle phase scores (S and G2M phases)
#'   \item Unique cell identifiers incorporating sample names
#' }
#'
#' @details This function:
#' \enumerate{
#'   \item Reads 10X Genomics data for each sample using \code{Seurat::Read10X}
#'   \item Creates and merges Seurat objects
#'   \item Adds quality control metrics:
#'   \itemize{
#'     \item Mitochondrial gene percentage (\code{mito_percent})
#'     \item Ribosomal gene percentage (\code{ribo_percent})
#'     \item Cell cycle scoring using pre-defined S and G2M phase genes
#'   }
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
#' \dontrun{
#' samples <- c("sample1", "sample2")
#' seurat_obj <- preprocess(
#'   samples = samples,
#'   path = "~/Desktop/TLS/data/raw-spatial/E-MTAB-13084",
#'   folder = "spatial"
#' )
#' }
#'
#' @importFrom Seurat Read10X CreateSeuratObject RenameCells
#' @importFrom Seurat PercentageFeatureSet CellCycleScoring
#' @importFrom SeuratObject JoinLayers
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
  seurat <- Seurat::CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)



}

