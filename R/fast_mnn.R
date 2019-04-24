#' @include internal.R
#'
NULL

#' Run fastMNN
#'
#' @param object.list A list of Seurat objects
#' @param assay Assay to use, defaults to the default assay of the first object
#' @param nfeatures Number of features to return
#' @param reduction.name Name to store resulting DimReduc object as
#' @param reduction.key Key for resulting DimReduc
#' @param verbose Print messages from \code{\link[Seurat]{SelectIntegrationFeatures}}
#' @param ... Extra parameters passed to \code{\link[scran]{fastMNN}}
#'
#' @return A Seurat object merged from the objects in \code{object.list} and a new
#' DimReduc of name \code{reduction.name} (key set to \code{reduction.key}) with
#' corrected embeddings matrix; all other return values from
#' \code{\link[scran]{fastMNN}} are stored in the \code{tool} slot, accessible with
#' \code{\link[Seurat]{Tool}}
#'
#' @importFrom Seurat DefaultAssay SelectIntegrationFeatures
#' as.SingleCellExperiment CreateDimReducObject Tool<-
#'
#' @export
#'
#' @seealso \code{\link[scran]{fastMNN}} \code{\link[Seurat]{Tool}}
#'
RunFastMNN <- function(
  object.list,
  assay = NULL,
  nfeatures = 2000,
  reduction.name = 'mnn',
  reduction.key = 'DIM_',
  verbose = TRUE,
  ...
) {
  if (!requireNamespace(package = 'scran', quietly = TRUE)) {
    stop("Please install scran for fastMNN")
  }
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
  features <- SelectIntegrationFeatures(
    object.list = object.list,
    nfeatures = nfeatures,
    verbose = verbose
  )
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
    what = scran::fastMNN,
    args = c(
      objects.sce,
      list(...)
    )
  )
  rownames(x = out$corrected) <- colnames(x = integrated)
  integrated[[reduction.name]] <- CreateDimReducObject(
    embeddings = out$corrected,
    assay = DefaultAssay(object = integrated),
    key = reduction.key
  )
  out$corrected <- NULL
  Tool(object = integrated) <- out
  return(integrated)
}
