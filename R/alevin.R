#' @include internal.R
#'
NULL

#' @inheritParams Seurat::CreateSeuratObject
#' Load alevin Counts data from an EDS file
#'
#' This is a wrapper around ReadCounts from AlevinRtools
#'
#' @param file Path to `quants_mat.gz` file
#'
#' @return Returns a Seurat object with alevin counts
#' @seealso \code{\link[alevin]{alevin}}
#' @author Avi Srivastava
#' @references Srivastava, Avi, et al. "Alevin efficiently 
#' estimates accurate gene abundances from dscRNA-seq data." 
#' Genome biology 20.1 (2019): 65.
#' @export
ReadAlevin <- function(file) {
    CheckPackage(package = 'k3yavi/alevin-Rtools', repository = 'github')
    counts <- AlevinRtools::ReadCounts(file)

    return(CreateSeuratObject(counts))
}
