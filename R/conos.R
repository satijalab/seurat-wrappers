#' Extra conversions to Seurat objects
#'
#' @inheritParams Seurat::as.Seurat
#'
#' @rdname as.Seurat.extras
#' @name as.Seurat
#'
#' @seealso \code{\link[Seurat]{as.Seurat}}
#'
#' @aliases as.Seurat
#'
NULL

#' @param method Name of matching method graph was built using
#' @param reduction Name of graph embedding, if calculated
#' @param idents Name of clutering method to set as identity class
#'
#' @details
#' The \code{Conos} method for \code{\link[Seurat]{as.Seurat}} only works if all
#' samples are \code{Seurat} objects. The object is initially constructed by merging
#' all samples together using \code{\link[Seurat]{merge}}, any sample-level dimensional
#' reductions and graphs will be lost during the merge. Extra information is added
#' to the resulting Seurat object as follows:
#' \itemize{
#'   \item Pairwise alignments will be stored in miscellaneous data, as will any
#'   other miscellaneous information
#'   \item If a graph is present in the \code{graph} field, it will be stored as
#'   a \code{Graph} object, reordered to match cell order in the new \code{Seurat}
#'   object. It will be named "\code{DefaultAssay(SeuratObject)}_\code{method}"
#'   \item If an embedding is present in the \code{embedding} field as a
#'   \code{\link{matrix}}, it will be stored as a \code{DimReduc} object with the
#'   name \code{reduction} and a key value of "\code{toupper(reduction)}_"
#'   \item If the length of the \code{clusters} field is greater than zero,
#'   clustering information (\code{groups} field) will be added to object metadata.
#'   Extra information (\code{result} field) will be added to miscellaneous data
#'   with the name "conos.\code{clustering}.result"
#'   \item If present, the first clustering entry in the \code{clusters} field
#'   will be set as object identity classes
#' }
#'
#' @importFrom igraph get.adjacency
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Seurat as.Seurat Misc<- DefaultAssay as.Graph
#' CreateDimReducObject Idents<-
#'
#' @rdname as.Seurat.extras
#' @export
#' @method as.Seurat Conos
#'
as.Seurat.Conos <- function(
  x,
  method = 'mnn',
  reduction = 'largeVis',
  idents = names(x = x$clusters)[1],
  verbose = TRUE,
  ...
) {
  if (!all(sapply(X = x$samples, FUN = inherits, what = 'Seurat'))) {
    stop(
      "Converting a Conos object to a Seurat object requires that all samples are Seurat v3 objects",
      call. = FALSE
    )
  }
  if (verbose) {
    message("Merging ", length(x = x$samples), " samples")
  }
  object <- merge(x = x$samples[[1]], x$samples[2:length(x = x$samples)])
  # Add pairs
  if (length(x = x$pairs) > 0) {
    if (verbose) {
      message("Adding pairwise alignments to 'conos.pairs' in miscellaneous data")
    }
    Misc(object = object, slot = 'conos.pairs') <- x$pairs
  }
  # Add graph
  if (!is.null(x = x$graph)) {
    graph <- paste(DefaultAssay(object = object), method, sep = '_')
    message("Adding graph as '", graph, "'")
    object[[graph]] <- as.Graph(x = get.adjacency(graph = x$graph)[colnames(x = object), colnames(x = object)])
  }
  # Add graph embedding
  if (is.matrix(x = x$embedding)) {
    if (verbose) {
      message("Adding graph embedding as ", reduction)
    }
    object[[reduction]] <- suppressWarnings(expr = CreateDimReducObject(
      embeddings = x$embedding,
      assay = DefaultAssay(object = object),
      key = paste0(toupper(x = reduction), '_')
    ))
  }
  # Add clustering information
  if (length(x = x$clusters) > 0) {
    if (verbose) {
      message("Adding clustering information")
      pb <- txtProgressBar(min = 0, max = length(x = x$clusters), style = 3, file = stderr())
    }
    for (clustering in names(x = x$clusters)) {
      object[[clustering]] <- x$clusters[[clustering]]$groups
      clustering.misc <- paste('conos', clustering, 'result', sep = '.')
      Misc(object = object, slot = clustering.misc) <- x$clusters[[clustering]]$result
      if (clustering == idents) {
        Idents(object = object) <- clustering
      }
      if (verbose) {
        setTxtProgressBar(pb = pb, value = 1 + pb$getVal())
      }
    }
  }
  if (verbose) {
    close(con = pb)
  }
  # Add miscellaneous information
  if (length(x = x$misc) > 0) {
    if (verbose) {
      message("Adding extra information to 'conos.misc' in miscellaneous data")
    }
    Misc(object = object, slot = 'conos.misc') <- x$misc
  }
  return(object)
}
