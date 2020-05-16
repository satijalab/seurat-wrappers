#' @include internal.R
#'
NULL

#' @inheritParams Seurat::CreateSeuratObject
#' Load Alevin Counts data from a EDS file
#'
#' This is a wrapper around ReadCounts from AlevinRtools
#'
#' @param file Path to `quants_mat.gz` file
#'
#' @importFrom utils capture.output
#'
#' @export
#'
#'
ReadAlevin <- function(file) {
    CheckPackage(package = 'k3yavi/alevin-Rtools', repository = 'github')
    counts <- AlevinRtools::ReadCounts(file)

    return(CreateSeuratObject(counts))
}
