#' Run LocalMAP (Local Manifold Approximation Projection)
#'
#' Runs LocalMAP (Wang, Rudin & Shaposhnik), a variant of PaCMAP that adds a
#' local graph-adjustment stage in phase 3: further-pair partners are
#' resampled to points already close in the low-D embedding (within
#' \code{low_dist_thres}) and the NN-term gradient is scaled by
#' \code{low_dist_thres / (2 sqrt(d_ij))}. Compared with base PaCMAP this
#' tightens local structure.
#'
#' This wrapper uses the native R + Rcpp package \pkg{pacmapr}
#' (\url{https://github.com/williamsyy/pacmap-for-R}). No Python / conda
#' installation is required.
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param ... Additional arguments forwarded to the default method / \code{pacmapr::localmap}.
#'
#' @author Yiyang Sun
#' @rdname RunLocalMAP
#' @export
RunLocalMAP <- function(object, ...) {
  if (inherits(object, "Seurat")) {
    RunLocalMAP.Seurat(object, ...)
  } else {
    RunLocalMAP.default(object, ...)
  }
}


#' @rdname RunLocalMAP
#' @method RunLocalMAP Seurat
#' @param reduction A character string specifying the reduction to be used as input. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param layer A character string specifying the layer name to be used. Default is "data".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "localmap".
#' @param reduction.key A character string specifying the prefix for the column names of the LocalMAP embeddings. Default is "LocalMAP_".
#'
#' @importFrom Seurat LogSeuratCommand DefaultAssay GetAssayData Embeddings
#' @export
RunLocalMAP.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                               assay = NULL, layer = "data",
                               n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                               distance_method = "euclidean",
                               lr = 1, num_iters = 250L, apply_pca = TRUE, init = "random",
                               low_dist_thres = 10,
                               reduction.name = "localmap", reduction.key = "LocalMAP_",
                               n_threads = NULL,
                               verbose = TRUE, seed.use = 11L, ...) {
  if (is.null(dims) && is.null(features)) {
    stop("Please specify one of `dims` or `features`.")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)

    data.use <- t(as.matrix(x = GetAssayData(object = object, layer = layer, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n_components) {
      stop(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " LocalMAP components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    if (!is.null(x = assay) && assay != DefaultAssay(object = object[[reduction]])) {
      warning("If both `assay` and `dims` are specified, the value of `assay` will get ignored.")
    }
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      stop(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " LocalMAP components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims or features")
  }
  object[[reduction.name]] <- RunLocalMAP(
    object = data.use, assay = assay,
    n_components = n_components, n.neighbors = n.neighbors,
    MN_ratio = MN_ratio, FP_ratio = FP_ratio,
    distance_method = distance_method,
    lr = lr, num_iters = num_iters, apply_pca = apply_pca, init = init,
    low_dist_thres = low_dist_thres,
    reduction.key = reduction.key, n_threads = n_threads,
    verbose = verbose, seed.use = seed.use, ...
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' @rdname RunLocalMAP
#' @method RunLocalMAP default
#' @importFrom Seurat CreateDimReducObject
#' @param n_components An integer specifying the number of LocalMAP components. Default is 2.
#' @param n.neighbors An integer specifying the number of neighbors considered in the k-Nearest Neighbor graph. Defaults to 10 for datasets with n <= 10000. For larger datasets the default is \code{round(10 + 15 * (log10(n) - 4))}.
#' @param MN_ratio A numeric value specifying the ratio of mid-near pairs to neighbor pairs. Default is 0.5.
#' @param FP_ratio A numeric value specifying the ratio of further pairs to neighbor pairs. Default is 2.
#' @param distance_method A character string specifying the distance metric to be used. One of "euclidean", "manhattan", "angular", "hamming". Default is "euclidean".
#' @param lr A numeric value specifying the Adam learning rate. Default is 1.
#' @param num_iters An integer or a length-3 integer vector giving the iterations for each of the three phases. A scalar is expanded to \code{c(100, 100, num_iters)}. Default is 250L (100 + 100 + 250 = 450 total iterations).
#' @param apply_pca A logical value indicating whether to apply PCA-to-100 preprocessing when the input has more than 100 features. Default is TRUE.
#' @param init A character string ("pca" or "random") or a numeric matrix used to initialize the low-dimensional embedding. Default is "random".
#' @param low_dist_thres Low-D distance threshold used for phase-3 FP resampling and as the LocalMAP NN-gradient coefficient. Default is 10 (Wang et al.).
#' @param reduction.key A character string specifying the prefix for the column names of the LocalMAP embeddings. Default is "LocalMAP_".
#' @param n_threads Number of threads for the ANN and gradient steps. Default (\code{NULL}) uses \code{parallel::detectCores() - 1L}.
#' @param verbose A logical value indicating whether to print progress. Default is TRUE.
#' @param seed.use An integer specifying the random seed (passed to \code{random_state}). Default is 11.
#' @param ... Additional arguments forwarded to \code{pacmapr::localmap}.
#' @export
RunLocalMAP.default <- function(object, assay = NULL,
                                n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                                distance_method = "euclidean",
                                lr = 1, num_iters = 250L, apply_pca = TRUE, init = "random",
                                low_dist_thres = 10,
                                reduction.key = "LocalMAP_",
                                n_threads = NULL,
                                verbose = TRUE, seed.use = 11L, ...) {
  if (!requireNamespace("pacmapr", quietly = TRUE)) {
    stop("Package 'pacmapr' is required. Install from GitHub with:\n",
         "  remotes::install_github(\"williamsyy/pacmap-for-R\", subdir = \"pacmapr\")",
         call. = FALSE)
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  X <- as.matrix(object)
  storage.mode(X) <- "double"

  res <- pacmapr::localmap(
    X              = X,
    n_components   = as.integer(n_components),
    n_neighbors    = if (is.null(n.neighbors)) NULL else as.integer(n.neighbors),
    MN_ratio       = MN_ratio,
    FP_ratio       = FP_ratio,
    distance       = distance_method,
    lr             = lr,
    num_iters      = num_iters,
    init           = init,
    apply_pca      = apply_pca,
    low_dist_thres = low_dist_thres,
    n_threads      = n_threads,
    random_state   = if (is.null(seed.use)) NULL else as.integer(seed.use),
    verbose        = verbose,
    ...
  )

  embedding <- res$embedding
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(X)

  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(reduction)
}
