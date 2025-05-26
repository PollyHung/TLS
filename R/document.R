#' Cell cycle marker genes (S phase)
#'
#' A character vector of genes associated with the S phase of the cell cycle.
#'
#' @format A character vector.
#' @source Processed from [original source, e.g., "Tirosh et al. 2016, Science"]
"s.genes"

#' Cell cycle marker genes (G2/M phase)
#'
#' A character vector of genes associated with the G2/M phase of the cell cycle.
#'
#' @format A character vector.
#' @source Processed from [original source, e.g., "Tirosh et al. 2016, Science"]
"g2m.genes"

#' TLS Signature Genes
#'
#' A character vector of 50 genes associated with tertiary lymphoid structures (TLS).
#'
#' @format A character vector.
#' @source Processed from [your source, e.g., "public dataset XYZ"]
"tls_50_genes"

#' T-cell marker genes
#'
#' A character vector of genes marking T-cell subsets and functional states.
#'
#' @format A character vector with 45 entries:
#' \describe{
#'   \item{CD3D/CD3E/CD3G}{Pan-T-cell markers}
#'   \item{CD4}{Helper T-cells}
#'   \item{CD8A/CD8B}{Cytotoxic T-cells}
#'   \item{FOXP3}{Regulatory T-cells (Tregs)}
#'   \item{PDCD1}{Exhaustion marker}
#' }
#' @source https://doi.org/10.1186/s13059-016-1070-5
"t_cell_markers"

#' B-cell marker genes
#'
#' A character vector of genes marking B-cell subsets and functional states.
#'
#' @format A character vector with 50 entries:
#' \describe{
#'   \item{CD19}{Pan-B-cell marker}
#'   \item{CD20}{B-cell activation marker}
#'   \item{CD22}{B-cell receptor signaling}
#' }
#' @source https://doi.org/10.1186/s13059-016-1070-5
"b_cell_markers"

#' Fibroblast markers
#'
#' A character vector of genes specific to fibroblast cells.
#'
#' @format A character vector with 10 entries:
#' \describe{
#'   \item{FAP}{Fibroblast activation protein}
#'   \item{COL1A1}{Collagen type I alpha 1 chain}
#' }
#' @source https://doi.org/10.1186/s13059-016-1070-5
"fibroblast"

#' NK cell marker genes
#'
#' A character vector of genes marking natural killer (NK) cell subsets.
#'
#' @format A character vector with 12 entries:
#' \describe{
#'   \item{NKG7}{NK cell activation marker}
#'   \item{CD56}{NK cell subset marker}
#' }
#' @source https://doi.org/10.1186/s13059-016-1070-5
"nk_cell_markers"

#' Myeloid markers
#'
#' A character vector of genes marking myeloid cell types.
#'
#' @format A character vector with 15 entries:
#' \describe{
#'   \item{CD14}{Monocyte marker}
#'   \item{CD68}{Macrophage marker}
#' }
#' @source https://doi.org/10.1186/s13059-016-1070-5
"myeloid_markers"

#' Endothelial markers
#'
#' A character vector of genes specific to endothelial cells.
#'
#' @format A character vector with 5 entries:
#' \describe{
#'   \item{CD31}{Endothelial cell adhesion molecule}
#'   \item{VEGFR}{Vascular endothelial growth factor receptor}
#' }
#' @source https://doi.org/10.1186/s13059-016-1070-5
"endothelial"

#' Follicular Dendritic Cell (FDC) Markers
#'
#' A character vector of genes specific to follicular dendritic cells.
#'
#' @format A character vector with 8 entries:
#' \describe{
#'   \item{CD21 (CR2)}{Complement receptor}
#'   \item{CD23 (FCER2)}{IgE receptor}
#'   \item{CXCL13}{Chemokine for B-cell follicle organization}
#' }
#' @source From DeepSeek
"follicular_dc"

#' Myeloid Dendritic Cell Markers
#'
#' A character vector of genes marking myeloid dendritic cell subsets.
#'
#' @format A character vector with 6 entries:
#' \describe{
#'   \item{CD1A/CD1B}{Antigen presentation}
#'   \item{CLEC10A}{cDC2 subset marker}
#' }
#' @source https://doi.org/10.1186/s13059-016-1070-5
"myeloid_dc"





