#' @include internal.R
#'
NULL

#' Run optimizeALS on a Seurat object
#'
#' @inheritParams liger::optimizeALS
#' @inheritParams RunFastMNN
#' @param object A merged Seurat object
#' @param split.by Attribute for splitting, defaults to "orig.ident"
#' @param ... Arguments passed to other methods
#'
#' @return A Seurat object with embeddings and loadings from \code{\link[liger]{optimizeALS}}
#' stored as a DimReduc object with name \code{reduction.name} (key set to \code{reduction.key});
#' per-dataset feature loadings matrices stored in the \code{tool} slot, accessible with
#' \code{\link[Seurat]{Tool}}
#'
#' @importFrom liger optimizeALS
#' @importFrom Seurat DefaultAssay SplitObject GetAssayData VariableFeatures
#' CreateDimReducObject Tool<- LogSeuratCommand
#'
#' @aliases optimizeALS
#' @seealso \code{\link[liger]{optimizeALS}} \code{\link[Seurat]{Tool}}
#'
#' @export
#' @method optimizeALS Seurat
#'
optimizeALS.Seurat <- function(
  object,
  k,
  assay = NULL,
  split.by = 'orig.ident',
  lambda = 5,
  thresh = 1e-4,
  max.iters = 100,
  reduction.name = '',
  reduction.key = '',
  nrep = 1,
  H.init = NULL,
  W.init = NULL,
  V.init = NULL,
  rand.seed = 1,
  print.obj = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = 'scale.data'))) {
    stop("Data is unscaled, splease scale before running", call. = FALSE)
  }
  scale.data <- sapply(
    X = SplitObject(object = object, split.by = split.by),
    FUN = GetAssayData,
    slot = 'scale.data',
    assay = assay,
    simplify = FALSE
  )
  scale.data <- sapply(X = scale.data, FUN = t, simplify = FALSE)
  out <- optimizeALS(
    object = scale.data,
    k = k,
    lambda = lambda,
    thresh = thresh,
    max.iters = max.iters,
    nrep = nrep,
    H.init = H.init,
    W.init = W.init,
    V.init = V.init,
    rand.seed = rand.seed,
    print.obj = print.obj
  )
  colnames(x = out$W) <- VariableFeatures(object = object)
  out$V <- sapply(
    X = out$V,
    FUN = function(x) {
      colnames(x = x) <- VariableFeatures(object = object)
      return(x)
    },
    simplify = FALSE
  )
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = do.call(what = 'rbind', args = out$H),
    loadings = out$W,
    assay = assay,
    key = reduction.key
  )
  Tool(object = object) <- out$V
  object <- LogSeuratCommand(object = object)
  return(object)
}
