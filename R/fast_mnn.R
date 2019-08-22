#' @include internal.R
#'
NULL

#' Run fastMNN
#'
#' @param object.list A list of Seurat objects
#' @param assay Assay to use, defaults to the default assay of the first object
#' @param features Either a list of features to use when calculating batch correction, or a number (2000 by default) of variable features to select.
#' @param reduction.name Name to store resulting DimReduc object as
#' @param reduction.key Key for resulting DimReduc
#' @param verbose Print messages from \code{\link[Seurat]{SelectIntegrationFeatures}}
#' @param ... Extra parameters passed to \code{\link[batchelor]{fastMNN}}
#'
#' @return A Seurat object merged from the objects in \code{object.list} and a new
#' DimReduc of name \code{reduction.name} (key set to \code{reduction.key}) with
#' corrected embeddings matrix; all other return values from
#' \code{\link[batchelor]{fastMNN}} are stored in the \code{tool} slot, accessible with
#' \code{\link[Seurat]{Tool}}
#'
#' @importFrom Seurat DefaultAssay DefaultAssay<- SelectIntegrationFeatures
#' as.SingleCellExperiment CreateDimReducObject Tool<- LogSeuratCommand
#'
#' @export
#'
#' @seealso \code{\link[batchelor]{fastMNN}} \code{\link[Seurat]{Tool}}
#'
RunFastMNN <- function(
  object.list,
  assay = NULL,
  features = 2000,
  reduction.name = 'mnn',
  reduction.key = 'mnn_',
  verbose = TRUE,
  ...
) {
  CheckPackage(package = 'batchelor', repository = 'bioconductor')
  if (!all(sapply(X = object.list, FUN = inherits, what = 'Seurat'))) {
    stop("'object.list' must be a list of Seurat objects", call. = FALSE)
  }
  if (length(x = object.list) < 2) {
    stop("'object.list' must contain multiple Seurat objects for integration", call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object.list[[1]])
  for (i in 1:length(x = object.list)) {
    DefaultAssay(object = object.list[[i]]) <- assay
  }
  if (is.numeric(x = features)) {
    if (verbose) {
      message(paste("Computing", features, "integration features"))
    }
    features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = features,
      assay = rep(assay,length(object.list))
    )
  }
  objects.sce <- lapply(
    X = object.list,
    FUN = function(x, f) {
      return(as.SingleCellExperiment(x = subset(x = x, features = f)))
    },
    f = features
  )
  integrated <- merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  )
  out <- do.call(
    what = batchelor::fastMNN,
    args = c(
      objects.sce,
      list(...)
    )
  )
  rownames(x = SingleCellExperiment::reducedDim(x = out)) <- colnames(x = integrated)
  colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key, 1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
  integrated[[reduction.name]] <- CreateDimReducObject(
    embeddings = SingleCellExperiment::reducedDim(x = out),
    assay = DefaultAssay(object = integrated),
    key = reduction.key
  )
  Tool(object = integrated) <- out
  integrated <- LogSeuratCommand(object = integrated)
  return(integrated)
}
