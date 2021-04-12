#' @include internal.R
#'
NULL

#' Run optimizeALS on a Seurat object
#'
#' @inheritParams rliger::optimizeALS
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
# @importFrom rliger optimizeALS
#' @importFrom Seurat DefaultAssay SplitObject GetAssayData VariableFeatures
#' CreateDimReducObject Tool<- LogSeuratCommand
#'
#' @aliases optimizeALS
#' @seealso \code{\link[rliger]{optimizeALS}} \code{\link[Seurat]{Tool}}
#'
#' @export
# @method optimizeALS Seurat
#'
RunOptimizeALS <- function(
  object,
  k,
  assay = NULL,
  split.by = 'orig.ident',
  lambda = 5,
  thresh = 1e-6,
  max.iters = 30,
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
  CheckPackage(package = 'rliger', repository = 'cran')
  assay <- assay %||% DefaultAssay(object = object)
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = 'scale.data'))) {
    stop("Data is unscaled, splease scale before running", call. = FALSE)
  }
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
  out <- rliger::optimizeALS(
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
#' This is a deprecated function. Call 'RunQuantileNorm' instead.
#'
# @inheritParams rliger::SNF
#' @inheritParams RunOptimizeALS
#' @param reduction Name of reduction to use
#'
#' @return A Seurat object with the SNF list stored in the \code{tool} slot,
#' accessible with \code{\link[Seurat]{Tool}}
#'
#' @importFrom Seurat SplitObject Embeddings Tool<- LogSeuratCommand
#'
#' @aliases SNF
#' @seealso \code{\link[rliger]{RunQuantileNorm}} \code{\link[Seurat]{Tool}}
#'
#' @export
# @method SNF Seurat
#'
RunSNF <- function(
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
  CheckPackage(package = 'rliger', repository = 'cran')
  # cells <- sapply(
  #   X = SplitObject(object = object, split.by = split.by),
  #   FUN = colnames,
  #   simplify = FALSE
  # )
  # dims.use <- dims.use %||% 1:length(x = object[[reduction]])
  # embeddings <- sapply(
  #   X = cells,
  #   FUN = function(x) {
  #     return(Embeddings(object = object[[reduction]])[x, ])
  #   },
  #   simplify = FALSE,
  #   USE.NAMES = TRUE
  # )
  # snf <- liger::SNF(
  #   object = embeddings,
  #   dims.use = dims.use,
  #   dist.use = dist.use,
  #   center = center,
  #   knn_k = knn_k,
  #   k2 = k2,
  #   small.clust.thresh = small.clust.thresh,
  #   ...
  # )
  # Tool(object = object) <- snf
  # object <- LogSeuratCommand(object = object)
  # return(object)
  .Deprecated(
    new = 'RunQuantileNorm',
    msg = paste(
      "This is a deprecated function. Call 'RunQuantileNorm' instead."
    )
  )
}

#' Run quantileAlignSNF on a Seurat object
#'
#' This is a deprecated function. Call 'RunQuantileNorm' instead.
#' 
#' @inheritParams RunSNF
#' @inheritParams RunOptimizeALS
#' @inheritParams rliger::quantileAlignSNF
#' @param recalc.snf Recalculate \code{\link{SNF}}
#' @param ... Arguments passed to other methods, and to
#' \code{\link[seurat.wrappers]{SNF}} if \code{recalc.snf = TRUE} or
#' \code{\link[seurat.wrappers]{SNF}} hasn't been run
#'
#' @return A Seurat object with embeddings from \code{\link[liger]{quantileAlignSNF}}
#' stored as a DimReduc object with name \code{reduction.name} (key set to \code{reduction.key})
#'
# @importFrom rliger quantileAlignSNF
#' @importFrom Seurat Tool SplitObject Embeddings CreateDimReducObject
#' DefaultAssay Tool<- Idents<- LogSeuratCommand
#'
#' @aliases quantileAlignSNF
#' @seealso \code{\link[rliger]{RunQuantileNorm}}
#'
#' @export
# @method quantileAlignSNF Seurat
#'
RunQuantileAlignSNF <- function(
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
  # CheckPackage(package = 'rliger', repository = 'cran')
  # if (recalc.snf || is.null(x = Tool(object = object, slot = 'RunSNF'))) {
  #   object <- RunSNF(
  #     object = object,
  #     split.by = split.by,
  #     reduction = reduction,
  #     center = center,
  #     ...
  #   )
  # }
  # embeddings <- sapply(
  #   X = SplitObject(object = object, split.by = split.by),
  #   FUN = function(x) {
  #     return(Embeddings(object = x[[reduction]]))
  #   },
  #   simplify = FALSE,
  #   USE.NAMES = TRUE
  # )
  # if (is.null(x = ref_dataset)) {
  #   num.samples <- vapply(
  #     X = embeddings,
  #     FUN = nrow,
  #     FUN.VALUE = integer(length = 1L)
  #   )
  #   ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  # } else if (is.numeric(x = ref_dataset)) {
  #   ref_dataset <- names(x = embeddings)[ref_dataset]
  # }
  # if (is.character(x = ref_dataset) && !ref_dataset %in% names(x = embeddings)) {
  #   stop("Cannot find reference dataset '", ref_dataset, "' in the split", call. = FALSE)
  # }
  # out <- rliger::quantileAlignSNF(
  #   object = embeddings,
  #   snf = Tool(object = object, slot = 'RunSNF'),
  #   cell.names = colnames(x = object),
  #   ref_dataset = ref_dataset,
  #   prune.thresh = prune.thresh,
  #   min_cells = min_cells,
  #   quantiles = quantiles,
  #   nstart = nstart,
  #   resolution = resolution,
  #   center = center,
  #   id.number = id.number,
  #   print.mod = print.mod,
  #   print.align.summary = print.align.summary,
  #   ...
  # )
  # object[[reduction.name]] <- CreateDimReducObject(
  #   embeddings = out$H.norm,
  #   assay = DefaultAssay(object = object[[reduction]]),
  #   key = reduction.key
  # )
  # out <- as.data.frame(x = out[names(x = out) != 'H.norm'])
  # object[[colnames(x = out)]] <- out
  # Idents(object = object) <- 'clusters'
  # object <- LogSeuratCommand(object = object)
  # return(object)
  message(paste(
      "This is a deprecated function. Calling 'RunQuantileNorm' instead.",
      "Note that not all parameters can be passed to 'RunQuantileNorm'.",
      "It's suggested to run 'louvainCluster' subsequently as well."
    ))

  .Deprecated(
    new = 'RunQuantileNorm',
    msg = paste(
      "This is a deprecated function. Calling 'quantile_norm' instead.",
      "Note that not all parameters can be passed to 'quantile_norm'.",
      "It's suggested to run 'louvainCluster' subsequently as well."
    )
  )

  return(RunQuantileNorm(object,
    split.by = split.by,
    reduction = reduction,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    quantiles = quantiles,
    ref_dataset = NULL,
    min_cells = 20,
    knn_k = 20,
    dims.use = NULL,
    do.center = FALSE,
    max_sample = 1000,
    eps = 0.9,
    refine.knn = TRUE,
    ...
  ))
}

#' Run quantile_norm on a Seurat object
#'
#' @inheritParams RunOptimizeALS
#' @inheritParams rliger::quantile_norm
#' @param ... Arguments passed to other methods
#'
#' @return A Seurat object with embeddings from \code{\link[liger]{quantile_norm}}
#' stored as a DimReduc object with name \code{reduction.name} (key set to \code{reduction.key})
#'
# @importFrom rliger quantile_norm
#' @importFrom Seurat Tool SplitObject Embeddings CreateDimReducObject
#' DefaultAssay Tool<- Idents<- LogSeuratCommand
#'
#' @aliases quantile_norm
#' @seealso \code{\link[rliger]{quantile_norm}}
#'
#' @export
# @method quantile_norm Seurat
#'
RunQuantileNorm <- function(
  object,
  split.by = 'orig.ident',
  reduction = 'iNMF_raw',
  reduction.name = 'iNMF',
  reduction.key = 'iNMF_',
  quantiles = 50,
  ref_dataset = NULL,
  min_cells = 20,
  knn_k = 20,
  dims.use = NULL,
  do.center = FALSE,
  max_sample = 1000,
  eps = 0.9,
  refine.knn = TRUE,
  ...
) {
  CheckPackage(package = 'rliger', repository = 'cran')
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
  out <- rliger::quantile_norm(
    object = embeddings,
    quantiles = quantiles,
    ref_dataset = ref_dataset,
    min_cells = min_cells,
    knn_k = knn_k,
    dims.use = dims.use,
    do.center = do.center,
    max_sample = max_sample,
    eps = eps,
    refine.knn = refine.knn,
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
