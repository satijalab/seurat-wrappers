#' Run Adaptively-thresholded Low Rank Approximation (ALRA)
#'
#' Runs ALRA, a method for imputation of dropped out values in scRNA-seq data.
#' Computes the k-rank approximation to A_norm and adjusts it according to the
#' error distribution learned from the negative values. Described in
#' Linderman, G. C., Zhao, J., Kluger, Y. (2018). "Zero-preserving imputation
#' of scRNA-seq data using low rank approximation." (bioRxiv:138677)
#'
#' @param object An object
#' @param k  The rank of the rank-k approximation. Set to NULL for automated choice of k.
#' @param q  The number of additional power iterations in randomized SVD when
#' computing rank k approximation. By default, q=10.
#' @param quantile.prob The quantile probability to use when calculating threshold.
#' By default, quantile.prob = 0.001.
#' @param use.mkl Use the Intel MKL based implementation of SVD. Needs to be
#' installed from https://github.com/KlugerLab/rpca-mkl.
#' @param mkl.seed Only relevant if use.mkl=T. Set the seed for the random
#' generator for the Intel MKL implementation of SVD. Any number <0 will
#' use the current timestamp. If use.mkl=F, set the seed using
#' set.seed() function as usual.
#' @param assay Assay to use
#' @param slot slot to use
#' @param setDefaultAssay If TRUE, will set imputed results as default Assay
#' @param genes.use genes to impute
#' @param K Number of singular values to compute when choosing k. Must be less
#' than the smallest dimension of the matrix. Default 100 or smallest dimension.
#' @param p.val.th  The threshold for ''significance'' when choosing k. Default 1e-10.
#' @param noise.start  Index for which all smaller singular values are considered noise.
#' Default K - 20.
#' @param q.k  Number of additional power iterations when choosing k. Default 2.
#' @param k.only If TRUE, only computes optimal k WITHOUT performing ALRA
#'
#' @param ... Arguments passed to other methods
#'
#' @importFrom rsvd rsvd
#' @importFrom Matrix Matrix
#' @importFrom stats pnorm sd setNames quantile
#' @importFrom Seurat DefaultAssay Tool GetAssayData Tool<- CreateAssayObject
#' DefaultAssay<-
#'
#'
#' @rdname RunALRA
#' @export RunALRA
#'
#' @author Jun Zhao, George Linderman
#' @references Linderman, G. C., Zhao, J., Kluger, Y. (2018). "Zero-preserving imputation
#' of scRNA-seq data using low rank approximation." (bioRxiv:138677)
#' @seealso \code{\link{ALRAChooseKPlot}}
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' # Example 1: Simple usage, with automatic choice of k.
#' pbmc_small_alra <- RunALRA(object = pbmc_small)
#' # Example 2: Visualize choice of k, then run ALRA
#' # First, choose K
#' pbmc_small_alra <- RunALRA(pbmc_small, k.only=TRUE)
#' # Plot the spectrum, spacings, and p-values which are used to choose k
#' ggouts <- ALRAChooseKPlot(pbmc_small_alra)
#' do.call(gridExtra::grid.arrange, c(ggouts, nrow=1))
#' # Run ALRA with the chosen k
#' pbmc_small_alra <- RunALRA(pbmc_small_alra)
#' }
#'
RunALRA <- function(object, ...) {
  UseMethod(generic = 'RunALRA', object = object)
}

#' @rdname RunALRA
#' @export
#'
RunALRA.default <- function(
  object,
  k = NULL,
  q = 10,
  quantile.prob = 0.001,
  use.mkl = FALSE,
  mkl.seed = -1,
  ...
) {
  A.norm <- t(x = as.matrix(x = object))
  message("Identifying non-zero values")
  originally.nonzero <- A.norm > 0
  message("Computing Randomized SVD")
  if (use.mkl) {
    CheckPackage(package = 'KlugerLab/rpca-mkl/fastRPCA', repository = 'github')
    fastDecomp.noc <- setNames(
      object = vector(mode = "list", length = 3),
      nm = c("u", "d", "v")
    )
    fastPCAOut <- fastRPCA::fastPCA(
      inputMatrix = A.norm,
      k = k,
      its = q,
      l = (k + 10),
      seed = mkl.seed
    )
    fastDecomp.noc$u <- fastPCAOut$U
    fastDecomp.noc$v <- fastPCAOut$V
    fastDecomp.noc$d <- diag(x = fastPCAOut$S)
  } else {
    fastDecomp.noc <- rsvd(A = A.norm, k = k, q = q)
  }
  A.norm.rank.k <- fastDecomp.noc$u[, 1:k] %*%
    diag(x = fastDecomp.noc$d[1:k]) %*%
    t(x = fastDecomp.noc$v[,1:k])
  message(sprintf("Find the %f quantile of each gene", quantile.prob))
  # A.norm.rank.k.mins <- abs(x = apply(X = A.norm.rank.k, MARGIN = 2, FUN = min))
  A.norm.rank.k.mins <- abs(x = apply(
    X = A.norm.rank.k,
    MARGIN = 2,
    FUN = function(x) {
      return(quantile(x = x, probs = quantile.prob))
    }
  ))
  message("Thresholding by the most negative value of each gene")
  A.norm.rank.k.cor <- replace(
    x = A.norm.rank.k,
    list = A.norm.rank.k <= A.norm.rank.k.mins[col(A.norm.rank.k)],
    values = 0
  )
  sd.nonzero <- function(x) {
    return(sd(x = x[!x == 0]))
  }
  sigma.1 <- apply(X = A.norm.rank.k.cor, MARGIN = 2, FUN = sd.nonzero)
  sigma.2 <- apply(X = A.norm, MARGIN = 2, FUN = sd.nonzero)
  mu.1 <- colSums(x = A.norm.rank.k.cor) / colSums(x = !!A.norm.rank.k.cor)
  mu.2 <- colSums(x = A.norm) / colSums(x = !!A.norm)
  toscale <- !is.na(sigma.1) & !is.na(sigma.2) & !(sigma.1 == 0 & sigma.2 == 0) & !(sigma.1 == 0)
  message(sprintf(fmt = "Scaling all except for %d columns", sum(!toscale)))
  sigma.1.2 <- sigma.2 / sigma.1
  toadd <- -1 * mu.1 * sigma.2 / sigma.1 + mu.2
  A.norm.rank.k.temp <- A.norm.rank.k.cor[, toscale]
  A.norm.rank.k.temp <- sweep(
    x = A.norm.rank.k.temp,
    MARGIN = 2,
    STATS = sigma.1.2[toscale],
    FUN = "*"
  )
  A.norm.rank.k.temp <- sweep(
    x = A.norm.rank.k.temp,
    MARGIN = 2,
    STATS = toadd[toscale],
    FUN = "+"
  )
  A.norm.rank.k.cor.sc <- A.norm.rank.k.cor
  A.norm.rank.k.cor.sc[, toscale] <- A.norm.rank.k.temp
  A.norm.rank.k.cor.sc[A.norm.rank.k.cor == 0] <- 0
  lt0 <- A.norm.rank.k.cor.sc < 0
  A.norm.rank.k.cor.sc[lt0] <- 0
  message(sprintf(
    fmt = "%.2f%% of the values became negative in the scaling process and were set to zero",
    100 * sum(lt0) / prod(dim(x = A.norm))
  ))
  A.norm.rank.k.cor.sc[originally.nonzero & A.norm.rank.k.cor.sc == 0] <-
    A.norm[originally.nonzero & A.norm.rank.k.cor.sc == 0]
  colnames(x = A.norm.rank.k) <- colnames(x = A.norm.rank.k.cor.sc) <-
    colnames(x = A.norm.rank.k.cor) <- colnames(x = A.norm)
  original.nz <- sum(A.norm > 0) / prod(dim(x = A.norm))
  completed.nz <- sum(A.norm.rank.k.cor.sc > 0) / prod(dim(x = A.norm))
  message(sprintf(
    fmt = "The matrix went from %.2f%% nonzero to %.2f%% nonzero",
    100 * original.nz,
    100 * completed.nz
  ))
  return(A.norm.rank.k.cor.sc)
}

#' @rdname RunALRA
#' @export
#' @method RunALRA Seurat
#'
RunALRA.Seurat <- function(
  object,
  k = NULL,
  q = 10,
  quantile.prob = 0.001,
  use.mkl = FALSE,
  mkl.seed=-1,
  assay = NULL,
  slot = "data",
  setDefaultAssay = TRUE,
  genes.use = NULL,
  K = NULL,
  thresh=6,
  noise.start = NULL,
  q.k = 2,
  k.only = FALSE,
  ...
) {
  if (!is.null(x = k) && k.only) {
    warning("Stop: k is already given, set k.only = FALSE or k = NULL")
  }
  genes.use <- genes.use %||% rownames(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  alra.previous <- Tool(object = object, slot = 'RunALRA')
  alra.info <- list()
  # Check if k is already stored
  if (is.null(x = k) & !is.null(alra.previous[["k"]])) {
    k <- alra.previous[["k"]]
    message("Using previously computed value of k")
  }
  data.used <- GetAssayData(object = object, assay = assay, slot = slot)[genes.use,]
  # Choose k with heuristics if k is not given
  if (is.null(x = k)) {
    # set K based on data dimension
    if (is.null(x = K)) {
      K <- 100
      if (K > min(dim(x = data.used))) {
        K <- min(dim(x = data.used))
        warning("For best performance, we recommend using ALRA on expression matrices larger than 100 by 100")
      }
    }
    if (K > min(dim(x = data.used))) {
      stop("For an m by n data, K must be smaller than the min(m,n)")
    }
    # set noise.start based on K
    if (is.null(x = noise.start)) {
      noise.start <- K - 20
      if (noise.start <= 0) {
        noise.start <- max(K - 5, 1)
      }
    }
    if (noise.start > K - 5) {
      stop("There need to be at least 5 singular values considered noise")
    }
    noise.svals <- noise.start:K
    if (use.mkl) {
      CheckPackage(package = 'KlugerLab/rpca-mkl/fastRPCA', repository = 'github')
      L <- min(K + 10, min(dim(x = data.used)))
      rsvd.out <- setNames(
        object = vector(mode = "list", length = 3),
        nm = c("u", "d", "v")
      )
      fastPCAOut <- fastRPCA::fastPCA(
        inputMatrix = as.matrix(x = data.used),
        k = K,
        its = q.k,
        l = L,
        seed = mkl.seed
      )
      rsvd.out$u <- fastPCAOut$U
      rsvd.out$v <- fastPCAOut$V
      rsvd.out$d <- diag(x = fastPCAOut$S)
    } else {
      rsvd.out <- rsvd(A = t(x = as.matrix(x = data.used)), k = K, q = q.k)
    }
    diffs <- rsvd.out$d[1:(length(x = rsvd.out$d) - 1)] - rsvd.out$d[2:length(x = rsvd.out$d)]
    mu <- mean(x = diffs[noise.svals - 1])
    sigma <- sd(x = diffs[noise.svals - 1])
    num_of_sds <- (diffs - mu) / sigma
    k <- max(which(x = num_of_sds > thresh))
    alra.info[["d"]] <- rsvd.out$d
    alra.info[["k"]] <- k
    alra.info[["diffs"]] <- diffs
    Tool(object = object) <- alra.info
  }
  if (k.only) {
    message("Chose rank k = ", k, ", WITHOUT performing ALRA")
    return(object)
  }
  message("Rank k = ", k)
  # Perform ALRA on data.used
  output.alra <- RunALRA(
    object = data.used,
    k = k,
    q = q,
    quantile.prob = quantile.prob,
    use.mkl = use.mkl,
    mkl.seed = mkl.seed
  )
  # Save ALRA data in object@assay
  data.alra <- Matrix(data = t(x = output.alra), sparse = TRUE)
  rownames(x = data.alra) <- genes.use
  colnames(x = data.alra) <- colnames(x = object)
  assay.alra <- CreateAssayObject(data = data.alra)
  object[["alra"]] <- assay.alra
  if (setDefaultAssay) {
    message("Setting default assay as alra")
    DefaultAssay(object = object) <- "alra"
  }
  return(object)
}


#' ALRA Approximate Rank Selection Plot
#'
#' Plots the results of the approximate rank selection process for ALRA.
#'
#'
#' @param object Seurat object
#' @param start Index to start plotting singular value spacings from.
#' The transition from "signal" to "noise" in the is hard to see because the
#' first singular value spacings are so large. Nicer visualizations result from
#' skipping the first few. If set to 0 (default) starts from k/2.
#' @param combine Combine plots into a single gg object; note that if TRUE,
#' themeing will not work when plotting multiple features
#'
#' @return A list of 3 ggplot objects splotting the singular values, the
#' spacings of the singular values, and the p-values of the singular values.
#'
#' @author Jun Zhao, George Linderman
#' @seealso \code{\link{RunALRA}}
#'
#' @importFrom Seurat CombinePlots
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line
#' geom_vline scale_x_continuous labs
#' @export
#'
ALRAChooseKPlot <- function(object, start = 0, combine = TRUE) {
  alra.data <- Tool(object = object, slot = 'RunALRA')
  if (is.null(x = alra.data)) {
    stop('RunALRA should be run prior to using this function.')
  }
  d <- alra.data[["d"]]
  diffs <- alra.data[["diffs"]]
  pvals <- alra.data[["pvals"]]
  k <- alra.data[["k"]]
  if (start == 0) {
    start <- floor(x = k / 2)
  }
  if (start > k) {
    stop("Plots should include k (i.e. starting.from should be less than k)")
  }
  breaks <- seq(from = 10, to = length(x = d), by = 10)
  ggdata <- data.frame(x = 1:length(x = d), y = d)
  gg1 <- ggplot(data = ggdata, mapping = aes_string(x = 'x', y = 'y')) +
    geom_point(size = 1) +
    geom_line(size = 0.5) +
    geom_vline(xintercept = k) +
    theme_cowplot() +
    scale_x_continuous(breaks = breaks) +
    labs(x = NULL, y = 's_i', title = 'Singular values')
  ggdata <- data.frame(x = 2:length(x = d), y = diffs)[-(1:(start - 1)), ]
  gg2 <- ggplot(data = ggdata, mapping = aes_string(x = 'x', y = 'y')) +
    geom_point(size = 1) +
    geom_line(size = 0.5) +
    geom_vline(xintercept = k + 1) +
    theme_cowplot() +
    scale_x_continuous(breaks = breaks) +
    labs(x = NULL, y = 's_{i} - s_{i-1}', title = 'Singular value spacings')
  ggdata <- data.frame(x = 2:length(x = d), y = pvals)
  gg3 <- ggplot(data = ggdata, mapping = aes_string(x = 'x', y = 'y')) +
    geom_point(size = 1) +
    geom_vline(xintercept = k + 1) +
    theme_cowplot() +
    scale_x_continuous(breaks = breaks) +
    labs(x = NULL, y = 'p.val', title = 'Singular value spacing p-values')
  plots <- list(spectrum = gg1, spacings = gg2, pvals = gg3)
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}
