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
#' @param getMeta logical, option to use \code{tximeta} to
#' programmatically obtain gene range information, default
#' is FALSE
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
ReadAlevin <- function(file, getMeta=FALSE, ...) {
  CheckPackage(package = 'tximport', repository = 'bioconductor')
  CheckPackage(package = 'fishpond', repository = 'bioconductor')
  hasTximeta <- requireNamespace("tximeta", quietly=TRUE)
  if (getMeta) {
    if (!hasTximeta) stop("tximeta is not installed, use BiocManager::install()")
    # TODO ... use tximeta here and stash ranges
  } else {
    txi <- tximport(file, type="alevin", ...)
  }
  obj <- CreateSeuratObject(txi$counts)
  if ("mean" %in% names(txi)) {
    obj <- obj
  }
  if ("variance" %in% names(txi)) {
    obj <- obj
  }
  if ("tier" %in% names(txi)) {
    obj <- obj
  }
  return(obj)
}
