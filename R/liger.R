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
  reduction.name = 'iNMF_raw',
  reduction.key = 'riNMF_',
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
  # scale.data <- sapply(
  #   X = SplitObject(object = object, split.by = split.by),
  #   FUN = GetAssayData,
  #   slot = 'scale.data',
  #   assay = assay,
  #   simplify = FALSE
  # )
  if (is.character(x = split.by) && length(x = split.by) == 1) {
    split.by <- object[[split.by]]
  }
  split.cells <- split(x = colnames(x = object), f = split.by)
  scale.data <- lapply(
    X = split.cells,
    FUN = function(x) {
      return(t(x = GetAssayData(
        object = object,
        slot = 'scale.data',
        assay = assay
      )[, x]))
    }
  )
  # scale.data <- sapply(X = scale.data, FUN = t, simplify = FALSE)
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
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = do.call(what = 'rbind', args = out$H),
    loadings = t(x = out$W),
    assay = assay,
    key = reduction.key
  )
  Tool(object = object) <- sapply(
    X = out$V,
    FUN = function(x) {
      colnames(x = x) <- VariableFeatures(object = object)
      rownames(x = x) <- colnames(x = object[[reduction.name]])
      return(t(x = x))
    },
    simplify = FALSE
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Generate shared factor neighborhood graph
#'
#' @inheritParams liger::SNF
#' @inheritParams optimizeALS.Seurat
#' @param reduction Name of reduction to use
#'
#' @return A Seurat object with the SNF list stored in the \code{tool} slot,
#' accessible with \code{\link[Seurat]{Tool}}
#'
#' @importFrom liger SNF
#' @importFrom Seurat SplitObject Embeddings Tool<- LogSeuratCommand
#'
#' @aliases SNF
#' @seealso \code{\link[liger]{SNF}} \code{\link[Seurat]{Tool}}
#'
#' @export
#' @method SNF Seurat
#'
SNF.Seurat <- function(
  object,
  split.by = 'orig.ident',
  reduction = 'iNMF_raw',
  dims.use = NULL,
  dist.use = 'CR',
  center = FALSE,
  knn_k = 20,
  k2 = 500,
  small.clust.thresh = knn_k,
  ...
) {
  cells <- sapply(
    X = SplitObject(object = object, split.by = split.by),
    FUN = colnames,
    simplify = FALSE
  )
  dims.use <- dims.use %||% 1:length(x = object[[reduction]])
  embeddings <- sapply(
    X = cells,
    FUN = function(x) {
      return(Embeddings(object = object[[reduction]])[x, ])
    }
  )
  snf <- SNF(
    object = embeddings,
    dims.use = dims.use,
    dist.use = dist.use,
    center = center,
    knn_k = knn_k,
    k2 = k2,
    small.clust.thresh = small.clust.thresh,
    ...
  )
  Tool(object = object) <- snf
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Run quantileAlignSNF on a Seurat object
#'
#' @inheritParams SNF.Seurat
#' @inheritParams optimizeALS.Seurat
#' @inheritParams liger::quantileAlignSNF
#' @param recalc.snf Recalculate \code{\link{SNF}}
#' @param ... Arguments passed to other methods, and to
#' \code{\link[seurat.wrappers]{SNF}} if \code{recalc.snf = TRUE} or
#' \code{\link[seurat.wrappers]{SNF}} hasn't been run
#'
#' @return A Seurat object with embeddings from \code{\link[liger]{quantileAlignSNF}}
#' stored as a DimReduc object with name \code{reduction.name} (key set to \code{reduction.key})
#'
#' @importFrom liger quantileAlignSNF
#' @importFrom Seurat Tool SplitObject Embeddings CreateDimReducObject
#' DefaultAssay Tool<- Idents<- LogSeuratCommand
#'
#' @aliases quantileAlignSNF
#' @seealso \code{\link[liger]{quantileAlignSNF}}
#'
#' @export
#' @method quantileAlignSNF Seurat
#'
quantileAlignSNF.Seurat <- function(
  object,
  split.by = 'orig.ident',
  reduction = 'iNMF_raw',
  reduction.name = 'iNMF',
  reduction.key = 'iNMF_',
  recalc.snf = FALSE,
  ref_dataset = NULL,
  prune.thresh = 0.2,
  min_cells = 2,
  quantiles = 50,
  nstart = 10,
  resolution = 1,
  center = FALSE,
  id.number = NULL,
  print.mod = FALSE,
  print.align.summary = FALSE,
  ...
) {
  if (recalc.snf || is.null(x = Tool(object = object, slot = 'SNF'))) {
    object <- SNF(object = object, ...)
  }
  embeddings <- sapply(
    X = SplitObject(object = object, split.by = split.by),
    FUN = function(x) {
      return(Embeddings(object = x[[reduction]]))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (is.null(x = ref_dataset)) {
    num.samples <- vapply(
      X = embeddings,
      FUN = nrow,
      FUN.VALUE = integer(length = 1L)
    )
    ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  } else if (is.numeric(x = ref_dataset)) {
    ref_dataset <- names(x = embeddings)[ref_dataset]
  }
  if (is.character(x = ref_dataset) && !ref_dataset %in% names(x = embeddings)) {
    stop("Cannot find reference dataset '", ref_dataset, "' in the split", call. = FALSE)
  }
  out <- quantileAlignSNF(
    object = embeddings,
    snf = Tool(object = object, slot = 'SNF'),
    cell.names = colnames(x = object),
    ref_dataset = ref_dataset,
    prune.thresh = prune.thresh,
    min_cells = min_cells,
    quantiles = quantiles,
    nstart = nstart,
    resolution = resolution,
    center = center,
    id.number = id.number,
    print.mod = print.mod,
    print.align.summary = print.align.summary,
    ...
  )
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = out$H.norm,
    assay = DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  out <- as.data.frame(x = out[names(x = out) != 'H.norm'])
  object[[colnames(x = out)]] <- out
  Idents(object = object) <- 'clusters'
  object <- LogSeuratCommand(object = object)
  return(object)
}
