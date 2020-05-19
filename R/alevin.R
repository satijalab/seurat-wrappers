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
#' is FALSE. Ranges are stored in \code{chr}, \code{start},
#' and \code{end} in the \code{meta.features} slot.
#' @param meanAndVariance logical, should mean and variance
#' of counts be returned in \code{counts} and \code{data}
#' slots, respectively
#' @param ... extra arguments passed to \code{tximport},
#' for example,
#' \code{alevinArgs=list(filterBarcodes=TRUE)}.
#'
#' @return returns a Seurat object with alevin counts
#' @seealso \code{\link[alevin]{alevin}}
#' @author Avi Srivastava
#' @references Srivastava, Avi, et al. "Alevin efficiently 
#' estimates accurate gene abundances from dscRNA-seq data." 
#' Genome biology 20.1 (2019): 65.
#' @export
ReadAlevin <- function(file, getMeta=FALSE, meanAndVariance=FALSE, ...) {
  CheckPackage(package = 'tximport', repository = 'bioconductor')
  CheckPackage(package = 'fishpond', repository = 'bioconductor')
  hasTximeta <- requireNamespace("tximeta", quietly=TRUE)
  metaSuccess <- FALSE
  if (getMeta) {
    if (!hasTximeta)
      stop("tximeta is not installed, use BiocManager::install()")
    se <- tximeta::tximeta(file, type="alevin", ...)
    metaSuccess <- !is.null(SummarizedExperiment::rowRanges(se))
    if (meanAndVariance &
        all(c("mean","variance") %in% SummarizedExperiment::assayNames(se))) {
      txi <- list(mean=SummarizedExperiment::assays(se)[["mean"]],
                  variance=SummarizedExperiment::assays(se)[["variance"]])
    } else {
      txi <- list(counts=SummarizedExperiment::assays(se)[["counts"]])
    }
  } else {
    txi <- tximport::tximport(file, type="alevin", ...)
  }
  if (meanAndVariance) {
    if (!all(c("mean","variance") %in% names(txi)))
      stop("mean and variance not present in alevin directory")
    obj <- CreateSeuratObject(counts=txi$mean)
    obj <- Seurat::SetAssayData(obj, "data", txi$variance)
  } else {
    obj <- CreateSeuratObject(counts=txi$counts)
  }
  if (metaSuccess) {
    r <- SummarizedExperiment::rowRanges(se)
    obj[["RNA"]][["chr"]] <- as.character(GenomicRanges::seqnames(r))
    obj[["RNA"]][["start"]] <- GenomicRanges::start(r)
    obj[["RNA"]][["end"]] <- GenomicRanges::end(r)
  }
  return(obj)
}
