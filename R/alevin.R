#' @include internal.R
#'
NULL

#' Load alevin quantification data
#'
#' A wrapper around tximport to create a \code{SeuratObject}
#' from alevin quantification data.
#'
#' @param file path to \code{quants_mat.gz} file within
#' alevin directory
#' @param ... extra arguments passed to \code{tximport},
#' for example,
#' \code{alevinArgs=list(filterBarcodes=TRUE, tierImport=TRUE)}.
#'
#' @return returns a Seurat object with alevin counts
#' @seealso \code{\link[alevin]{alevin}}
#' @author Avi Srivastava
#' @references Srivastava, Avi, et al. "Alevin efficiently 
#' estimates accurate gene abundances from dscRNA-seq data." 
#' Genome biology 20.1 (2019): 65.
#' @export
ReadAlevin <- function(file, ...) {
  CheckPackage(package = 'tximport', repository = 'bioconductor')
  txi <- tximport(file, type="alevin", ...)
  # TODO 1 - here we can also bring along other assays optionally
  # TODO 2 - here we can see if tximeta is available, and if so
  #          we could bring along the range information (once
  #          we've sorted out how to stash that properly)
  return(CreateSeuratObject(txi$counts))
}
