
#' @include internal.R
#'
NULL

#' Run Banksy on a Seurat Object
#'
#' @param object A Seurat object
#' @param lambda (numeric) Spatial weight parameter
#' @param assay (character) Assay in Seurat object to use
#' @param slot (character) Slot in Seurat assay to use
#' @param use_agf (boolean) Whether to use the AGF
#' @param dimx (character) Column name of spatial x dimension (must be in metadata)
#' @param dimy (character) Column name of spatial y dimension (must be in metadata)
#' @param dimz (character) Column name of spatial z dimension (must be in metadata)
#' @param ndim (integer) Number of spatial dimensions to extract
#' @param features (character) Features to compute. Can be 'all', 'variable' or
#'   a vector of feature names
#' @param group (character) Column name of a grouping variable (must be in metadata)
#' @param split.scale (boolean) Whether to separate scaling by group
#' @param k_geom (numeric) kNN parameter - number of neighbors to use
#' @param n (numeric) kNN_rn parameter - exponent of radius
#' @param sigma (numeric) rNN parameter - standard deviation of Gaussian kernel
#' @param alpha (numeric) rNN parameter - determines radius used
#' @param k_spatial (numeric) rNN parameter - number of neighbors to use
#' @param spatial_mode (character) Kernel for neighborhood computation
#' \itemize{
#'  \item{kNN_median: k-nearest neighbors with median-scaled Gaussian kernel}
#'  \item{kNN_r: k-nearest neighbors with $1/r$ kernel}
#'  \item{kNN_rn: k-nearest neighbors with $1/r^n$ kernel}
#'  \item{kNN_rank: k-nearest neighbors with rank Gaussian kernel}
#'  \item{kNN_unif: k-nearest neighbors wth uniform kernel}
#'  \item{rNN_gauss: radial nearest neighbors with Gaussian kernel}
#' }
#' @param assay_name (character) Name for output. When lazy=FALSE, names the
#'   BANKSY assay. When lazy=TRUE, names the dimensionality reduction.
#' @param M (numeric) Advanced usage. Highest azimuthal harmonic
#' @param chunk_size A integer scalar specifying the number of rows / genes of
#'   the neighborhood cell matrix to compute. Must be less than floor of
#'   2e31-1 / number of cells. This is automatically computed but can be
#'   specified.
#' @param parallel A logical scalar specifying whether to compute chunks in
#'   parallel using bplapply. Not implemented for Windows.
#' @param num_cores A integer scalar specifying the number of cores to use
#'   if parallel is TRUE.
#' @param lazy (boolean) If TRUE, compute PCA directly without materializing
#'   the full BANKSY matrix. Enables analysis of very large datasets (millions
#'   of cells) that would otherwise exceed available memory. Currently only
#'   supported for M=0 (no AGF). Supports within-group scaling with
#'   \code{split.scale=TRUE}. Requires the irlba package.
#' @param npcs (integer) Number of PCs to compute when lazy=TRUE (default 50)
#' @param scale_max (numeric) Maximum absolute z-score for clipping when
#'   lazy=TRUE. Matches Seurat's FastRowScale default of 10.
#' @param verbose (boolean) Print messages
#'
#' @return A Seurat object. If lazy=FALSE, contains a new assay holding the
#'   BANKSY matrix. If lazy=TRUE, contains a dimensionality reduction named
#'   by assay_name with PCA embeddings computed on the BANKSY matrix.
#'
#' @seealso \code{\link[Banksy]{ComputeBanksy}}
#'
#' @author Joseph Lee, Vipul Singhal
#'
#' @references Vipul Singhal, Nigel Chou et. al. BANKSY: A Spatial Omics
#' Algorithm that Unifies Cell Type Clustering and Tissue Domain Segmentation
#'
#' @export
RunBanksy <- function(object, lambda, assay='RNA', slot='data', use_agf=FALSE,
                      dimx=NULL, dimy=NULL, dimz=NULL, ndim=2,
                      features='variable',
                      group=NULL, split.scale=TRUE,
                      k_geom=15, n=2, sigma=1.5,
                      alpha=0.05, k_spatial=10, spatial_mode='kNN_median',
                      assay_name='BANKSY', M=NULL, chunk_size=NULL,
                      parallel=FALSE, num_cores=NULL,
                      lazy=FALSE, npcs=50L, scale_max=10,
                      verbose=TRUE) {
    # Check packages
    SeuratWrappers:::CheckPackage(package = 'data.table', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Matrix', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Banksy', repository = 'github')

    # Check lambda param
    if (lambda < 0 || lambda > 1) stop('Lambda must be between 0 and 1')

    # Get data
    data_own <- get_data(object, assay, slot, features, verbose)

    # Get locs
    locs <- get_locs(object, dimx, dimy, dimz, ndim, data_own, group, verbose)
    if (!is.null(group)) {
        object <- AddMetaData(
            object, metadata = locs,
            col.name = paste0('staggered_', colnames(locs)))
    }

    # Compute neighbor matrix
    knn_list <- lapply(k_geom, function(kg) {
      Banksy:::computeNeighbors(locs,
                                spatial_mode = spatial_mode, k_geom = kg, n = n,
                                sigma=sigma, alpha=alpha, k_spatial=k_spatial,
                                verbose=verbose)
    })

    # Resolve harmonics
    M <- seq(0, max(Banksy:::getM(use_agf, M)))

    if (!lazy && nrow(data_own) * ncol(data_own) > 1e8) {
        message('Note: dataset is large (', nrow(data_own), ' features x ',
                ncol(data_own), ' cells). Consider using lazy=TRUE for ',
                'memory-efficient PCA without materializing the full BANKSY matrix.')
    }

    if (lazy) {
        object <- .banksy_lazy_pca(
            object, data_own, knn_list, M, lambda, group, split.scale,
            assay, assay_name, npcs, scale_max, verbose)
    } else {
        object <- .banksy_standard(
            object, data_own, knn_list, M, lambda, group, split.scale,
            assay, slot, assay_name, verbose, chunk_size, parallel, num_cores)
    }

    # Log commands
    object <- Seurat::LogSeuratCommand(object = object)
    return(object)
}

# Lazy PCA path
.banksy_lazy_pca <- function(object, data_own, knn_list, M, lambda, group,
                             split.scale, assay, assay_name, npcs, scale_max,
                             verbose) {
    # Guards
    if (max(M) > 0) stop('lazy=TRUE currently only supports M=0 (no AGF)')
    if (!is.null(group) && !split.scale) {
        if (verbose) warning('lazy=TRUE: group is used for spatial ',
                             'staggering only; scaling is performed globally')
    }
    SeuratWrappers:::CheckPackage(package = 'irlba', repository = 'CRAN')

    n_genes <- nrow(data_own)
    n_cells <- ncol(data_own)
    lambdas <- Banksy:::getLambdas(lambda, n_harmonics = 1)

    # Determine split-scale grouping
    split_scale <- !is.null(group) && split.scale
    if (split_scale) {
        groups <- unlist(object[[group]])
        ugroups <- unique(groups)
        group_idx <- lapply(ugroups, function(g) which(g == groups))
    } else {
        group_idx <- NULL
    }

    # Validate npcs
    max_npcs <- min(2L * n_genes, n_cells) - 1L
    if (npcs >= max_npcs) {
        stop('npcs (', npcs, ') must be < min(2*n_genes, n_cells) = ',
             max_npcs + 1L)
    }
    if (min(2L * n_genes, n_cells) < 6L) {
        stop('lazy=TRUE requires at least 3 genes and 6 cells')
    }

    # Build sparse weight matrix
    if (verbose) message('Building sparse weight matrix')
    W <- Banksy:::.buildWeightMatrix(knn_list[[1]], n_cells)

    # Compute scaling parameters and clipping excess
    own_result <- .lazy_own_scaling(data_own, split_scale, group_idx,
                                    groups = if (split_scale) groups else NULL,
                                    ugroups = if (split_scale) ugroups else NULL,
                                    n_genes, n_cells, scale_max, verbose)
    h0_result <- .lazy_h0_scaling(data_own, W, split_scale, group_idx,
                                  n_genes, n_cells, scale_max, verbose)

    if (verbose) {
        n_own_clip <- if (!is.null(own_result$excess)) length(own_result$excess@x) else 0
        n_h0_clip <- if (!is.null(h0_result$excess)) length(h0_result$excess@x) else 0
        message('Clipping corrections: own=', n_own_clip, ' H0=', n_h0_clip, ' entries')
    }

    # Build operator struct
    banksy_op <- structure(list(
        gcm = data_own, W = W,
        mu = list(own_result$mu, h0_result$mu),
        sd = list(own_result$sd, h0_result$sd),
        lam = lambdas,
        excess = list(own_result$excess, h0_result$excess),
        valid = list(own_result$valid, h0_result$valid),
        split_scale = split_scale,
        group_idx = group_idx,
        n_genes = n_genes, n_cells = n_cells
    ), class = 'BanksyLazy')

    # Run irlba and store result
    if (verbose) message('Computing BANKSY PCA (', npcs, ' PCs) via lazy operator')
    pca <- irlba::irlba(banksy_op, nv = npcs, mult = .banksy_lazy_mult)

    # Cell embeddings: V * D
    embeddings <- sweep(pca$v, 2, pca$d, `*`)
    rownames(embeddings) <- colnames(data_own)
    colnames(embeddings) <- paste0(assay_name, '_', seq_len(npcs))

    # Feature loadings
    loadings <- pca$u
    feat_names <- c(rownames(data_own), paste0(rownames(data_own), '.m0'))
    rownames(loadings) <- feat_names
    colnames(loadings) <- paste0(assay_name, '_', seq_len(npcs))

    # Percent variance
    total_var <- sum(pca$d^2)
    stdev <- pca$d / sqrt(max(1, n_cells - 1))

    # Store as DimReduc
    dimreduc <- Seurat::CreateDimReducObject(
        embeddings = embeddings,
        loadings = loadings,
        stdev = stdev,
        key = paste0(assay_name, '_'),
        assay = assay,
        misc = list(total.variance = total_var)
    )
    object[[assay_name]] <- dimreduc

    if (verbose) message('Done. Access reduction with Reductions(obj, "',
                         assay_name, '")')
    object
}

# Own expression scaling params + clipping excess
.lazy_own_scaling <- function(data_own, split_scale, group_idx, groups,
                              ugroups, n_genes, n_cells, scale_max, verbose) {
    if (verbose) message('Computing scaling parameters for own expression')

    if (split_scale) {
        mu <- matrix(0, nrow = n_genes, ncol = length(group_idx))
        sd <- matrix(0, nrow = n_genes, ncol = length(group_idx))
        for (gr in seq_along(group_idx)) {
            cid <- group_idx[[gr]]
            n_c <- as.double(length(cid))
            if (inherits(data_own, 'sparseMatrix')) {
                mu[, gr] <- Matrix::rowMeans(data_own[, cid, drop = FALSE])
                sd[, gr] <- sqrt(pmax(
                    n_c / (n_c - 1) * (
                        Matrix::rowMeans(data_own[, cid, drop = FALSE]^2) -
                            mu[, gr]^2
                    ), 0
                ))
            } else {
                mu[, gr] <- rowMeans(data_own[, cid, drop = FALSE])
                sd[, gr] <- sqrt(pmax(
                    n_c / (n_c - 1) * (
                        rowMeans(data_own[, cid, drop = FALSE]^2) -
                            mu[, gr]^2
                    ), 0
                ))
            }
        }
        valid <- rowSums(sd == 0) == 0
        sd[sd == 0] <- 1
    } else {
        n_c <- as.double(n_cells)
        if (inherits(data_own, 'sparseMatrix')) {
            mu <- Matrix::rowMeans(data_own)
            sd <- sqrt(pmax(
                n_c / (n_c - 1) * (Matrix::rowMeans(data_own^2) - mu^2), 0
            ))
        } else {
            mu <- rowMeans(data_own)
            sd <- sqrt(pmax(
                n_c / (n_c - 1) * (rowMeans(data_own^2) - mu^2), 0
            ))
        }
        sd[sd == 0] <- 1
        valid <- NULL
    }

    # Compute clipping excess
    if (verbose) message('Computing clipping excess for own expression')
    excess <- .lazy_own_excess(data_own, mu, sd, valid, split_scale, group_idx,
                               groups, ugroups, n_genes, n_cells, scale_max)

    list(mu = mu, sd = sd, valid = valid, excess = excess)
}

.lazy_own_excess <- function(data_own, mu, sd, valid, split_scale, group_idx,
                             groups, ugroups, n_genes, n_cells, scale_max) {
    if (inherits(data_own, 'sparseMatrix')) {
        own_i <- data_own@i + 1L
        own_j <- rep(seq_len(n_cells), diff(data_own@p))
        if (split_scale) {
            own_gr <- match(groups[own_j], ugroups)
            own_mu <- mu[cbind(own_i, own_gr)]
            own_sd <- sd[cbind(own_i, own_gr)]
        } else {
            own_mu <- mu[own_i]
            own_sd <- sd[own_i]
        }
        exceed <- data_own@x > own_mu + scale_max * own_sd
        if (split_scale) exceed <- exceed & valid[own_i]
        if (any(exceed)) {
            ei <- own_i[exceed]
            z_vals <- (data_own@x[exceed] - own_mu[exceed]) / own_sd[exceed]
            Matrix::sparseMatrix(
                i = ei, j = own_j[exceed],
                x = z_vals - scale_max,
                dims = c(n_genes, n_cells))
        } else {
            NULL
        }
    } else {
        if (split_scale) {
            excess <- Matrix::sparseMatrix(
                i = integer(0), j = integer(0),
                dims = c(n_genes, n_cells))
            for (gr in seq_along(group_idx)) {
                cid <- group_idx[[gr]]
                z_own <- (data_own[, cid, drop = FALSE] - mu[, gr]) / sd[, gr]
                z_own[!valid, ] <- 0
                exceed_mask <- z_own > scale_max
                if (any(exceed_mask)) {
                    excess[, cid] <- Matrix::Matrix(
                        (z_own - scale_max) * exceed_mask, sparse = TRUE)
                }
            }
            if (length(excess@x) == 0) NULL else excess
        } else {
            z_own <- (data_own - mu) / sd
            exceed_mask <- z_own > scale_max
            if (any(exceed_mask)) {
                Matrix::Matrix((z_own - scale_max) * exceed_mask, sparse = TRUE)
            } else {
                NULL
            }
        }
    }
}

# H0 scaling params + two-pass clipping excess
.lazy_h0_scaling <- function(data_own, W, split_scale, group_idx,
                             n_genes, n_cells, scale_max, verbose) {
    if (verbose) message('Computing scaling params and clipping for H0')
    chunk_sz <- 100L

    # Pass 1: compute mu, ss, and per-gene max of H0 = data_own %*% W
    if (split_scale) {
        n_groups <- length(group_idx)
        mu <- matrix(0, nrow = n_genes, ncol = n_groups)
        ss <- matrix(0, nrow = n_genes, ncol = n_groups)
        max_h0 <- matrix(-Inf, nrow = n_genes, ncol = n_groups)
        for (ch_start in seq(1L, n_genes, by = chunk_sz)) {
            ch_end <- min(ch_start + chunk_sz - 1L, n_genes)
            ri <- ch_start:ch_end
            chunk <- as.matrix(data_own[ri, , drop = FALSE] %*% W)
            for (gr in seq_along(group_idx)) {
                cid <- group_idx[[gr]]
                chunk_group <- chunk[, cid, drop = FALSE]
                mu[ri, gr] <- rowMeans(chunk_group)
                ss[ri, gr] <- rowSums(chunk_group * chunk_group)
                max_h0[ri, gr] <- apply(chunk_group, 1, max)
            }
        }
        sd <- matrix(0, nrow = n_genes, ncol = n_groups)
        for (gr in seq_along(group_idx)) {
            n_c <- as.double(length(group_idx[[gr]]))
            sd[, gr] <- sqrt(pmax(
                n_c / (n_c - 1) * (ss[, gr] / n_c - mu[, gr]^2), 0
            ))
        }
        valid <- rowSums(sd == 0) == 0
        sd[sd == 0] <- 1
        thresh <- mu + scale_max * sd
    } else {
        mu <- numeric(n_genes)
        ss <- numeric(n_genes)
        max_h0 <- rep(-Inf, n_genes)
        for (ch_start in seq(1L, n_genes, by = chunk_sz)) {
            ch_end <- min(ch_start + chunk_sz - 1L, n_genes)
            ri <- ch_start:ch_end
            chunk <- as.matrix(data_own[ri, , drop = FALSE] %*% W)
            mu[ri] <- rowMeans(chunk)
            ss[ri] <- rowSums(chunk * chunk)
            max_h0[ri] <- apply(chunk, 1, max)
        }
        n_c <- as.double(n_cells)
        sd <- sqrt(pmax(n_c / (n_c - 1) * (ss / n_c - mu * mu), 0))
        sd[sd == 0] <- 1
        thresh <- mu + scale_max * sd
        valid <- NULL
    }

    # Identify genes that need clipping
    if (split_scale) {
        clip_genes <- which(valid & rowSums(max_h0 > thresh) > 0)
    } else {
        clip_genes <- which(max_h0 > thresh)
    }
    if (verbose) message('H0 genes requiring clipping: ', length(clip_genes),
                         ' / ', n_genes)

    # Pass 2: compute excess for clipped genes only
    excess <- .lazy_h0_excess(data_own, W, mu, sd, split_scale, group_idx,
                              clip_genes, chunk_sz, n_genes, n_cells, scale_max)

    list(mu = mu, sd = sd, valid = valid, excess = excess)
}

.lazy_h0_excess <- function(data_own, W, mu, sd, split_scale, group_idx,
                            clip_genes, chunk_sz, n_genes, n_cells, scale_max) {
    if (length(clip_genes) == 0) return(NULL)

    exc_cap <- max(1024L, length(clip_genes) * 10L)
    exc_i <- integer(exc_cap)
    exc_j <- integer(exc_cap)
    exc_x <- numeric(exc_cap)
    exc_n <- 0L

    clip_chunks <- unique((clip_genes - 1L) %/% chunk_sz)
    for (ch_idx in clip_chunks) {
        ch_start <- ch_idx * chunk_sz + 1L
        ch_end <- min(ch_start + chunk_sz - 1L, n_genes)
        ri <- ch_start:ch_end
        ri_clip <- ri[ri %in% clip_genes]
        chunk <- as.matrix(data_own[ri_clip, , drop = FALSE] %*% W)
        scan_idx <- if (split_scale) group_idx else list(seq_len(n_cells))
        for (gr in seq_along(scan_idx)) {
            cid <- scan_idx[[gr]]
            if (split_scale) {
                z_chunk <- (chunk[, cid, drop = FALSE] -
                    mu[ri_clip, gr]) / sd[ri_clip, gr]
            } else {
                z_chunk <- (chunk - mu[ri_clip]) / sd[ri_clip]
            }
            wh <- which(z_chunk > scale_max, arr.ind = TRUE)
            if (nrow(wh) > 0) {
                new_n <- nrow(wh)
                while (exc_n + new_n > length(exc_i)) {
                    exc_cap <- exc_cap * 2L
                    length(exc_i) <- exc_cap
                    length(exc_j) <- exc_cap
                    length(exc_x) <- exc_cap
                }
                idx <- seq(exc_n + 1L, exc_n + new_n)
                exc_i[idx] <- ri_clip[wh[, 1]]
                exc_j[idx] <- cid[wh[, 2]]
                exc_x[idx] <- z_chunk[wh] - scale_max
                exc_n <- exc_n + new_n
            }
        }
    }

    Matrix::sparseMatrix(
        i = exc_i[1:exc_n], j = exc_j[1:exc_n],
        x = exc_x[1:exc_n], dims = c(n_genes, n_cells))
}

# Lazy operator multiply (forward + adjoint)
#' @export
dim.BanksyLazy <- function(x) c(x$n_genes * 2L, x$n_cells)

.as_base <- function(x) {
    if (inherits(x, 'Matrix')) x <- as.matrix(x)
    if (!is.matrix(x)) x <- as.matrix(x)
    storage.mode(x) <- 'double'
    x
}

.banksy_lazy_mult <- function(A, x, transpose = FALSE) {
    if (inherits(x, 'BanksyLazy')) {
        tmp <- A; A <- x; x <- tmp; transpose <- TRUE
    }
    x <- .as_base(x)
    k <- ncol(x)

    if (!transpose) {
        .banksy_forward(A, x, k)
    } else {
        .banksy_adjoint(A, x, k)
    }
}

.banksy_forward <- function(A, x, k) {
    # A %*% x: x is n_cells x k
    # clip(Z) %*% x = Z %*% x - excess %*% x
    if (A$split_scale) {
        own <- matrix(0, nrow = A$n_genes, ncol = k)
        h0 <- matrix(0, nrow = A$n_genes, ncol = k)
        for (gr in seq_along(A$group_idx)) {
            cid <- A$group_idx[[gr]]
            x_group <- x[cid, , drop = FALSE]
            cs <- colSums(x_group)
            own_group <- .as_base(A$gcm[, cid, drop = FALSE] %*% x_group)
            own_group <- (own_group -
                outer(A$mu[[1]][, gr], cs)) / A$sd[[1]][, gr]
            if (!is.null(A$excess[[1]])) {
                own_group <- own_group - .as_base(
                    A$excess[[1]][, cid, drop = FALSE] %*% x_group)
            }
            own <- own + own_group

            Wx <- .as_base(A$W[, cid, drop = FALSE] %*% x_group)
            h0_group <- .as_base(A$gcm %*% Wx)
            h0_group <- (h0_group -
                outer(A$mu[[2]][, gr], cs)) / A$sd[[2]][, gr]
            if (!is.null(A$excess[[2]])) {
                h0_group <- h0_group - .as_base(
                    A$excess[[2]][, cid, drop = FALSE] %*% x_group)
            }
            h0 <- h0 + h0_group
        }
        own[!A$valid[[1]], ] <- 0
        h0[!A$valid[[2]], ] <- 0
        own <- A$lam[1] * own
        h0 <- A$lam[2] * h0
    } else {
        cs <- colSums(x)
        own <- .as_base(A$gcm %*% x)
        own <- A$lam[1] * (own - outer(A$mu[[1]], cs)) / A$sd[[1]]
        if (!is.null(A$excess[[1]])) {
            own <- own - A$lam[1] * .as_base(A$excess[[1]] %*% x)
        }
        Wx <- .as_base(A$W %*% x)
        h0 <- .as_base(A$gcm %*% Wx)
        h0 <- A$lam[2] * (h0 - outer(A$mu[[2]], cs)) / A$sd[[2]]
        if (!is.null(A$excess[[2]])) {
            h0 <- h0 - A$lam[2] * .as_base(A$excess[[2]] %*% x)
        }
    }
    rbind(own, h0)
}

.banksy_adjoint <- function(A, x, k) {
    # t(A) %*% x: x is (2*n_genes) x k
    # t(clip(Z)) %*% x = t(Z) %*% x - t(excess) %*% x
    ng <- A$n_genes
    xo <- x[1:ng, , drop = FALSE]
    xh <- x[(ng+1):(2*ng), , drop = FALSE]
    if (A$split_scale) {
        xo[!A$valid[[1]], ] <- 0
        xh[!A$valid[[2]], ] <- 0
        r <- matrix(0, nrow = A$n_cells, ncol = k)
        for (gr in seq_along(A$group_idx)) {
            cid <- A$group_idx[[gr]]
            xo_s <- xo / A$sd[[1]][, gr]
            xh_s <- xh / A$sd[[2]][, gr]
            adj_o <- colSums(A$mu[[1]][, gr] * xo_s)
            own <- .as_base(crossprod(
                A$gcm[, cid, drop = FALSE], xo_s)) -
                matrix(adj_o, length(cid), k, byrow = TRUE)
            if (!is.null(A$excess[[1]])) {
                own <- own - .as_base(crossprod(
                    A$excess[[1]][, cid, drop = FALSE], xo))
            }
            adj_h <- colSums(A$mu[[2]][, gr] * xh_s)
            ht <- .as_base(crossprod(A$gcm, xh_s))
            ht <- .as_base(crossprod(
                A$W[, cid, drop = FALSE], ht))
            ht <- ht - matrix(adj_h, length(cid), k, byrow = TRUE)
            if (!is.null(A$excess[[2]])) {
                ht <- ht - .as_base(crossprod(
                    A$excess[[2]][, cid, drop = FALSE], xh))
            }
            r[cid, ] <- A$lam[1] * own + A$lam[2] * ht
        }
    } else {
        xo_s <- xo / A$sd[[1]]
        xh_s <- xh / A$sd[[2]]
        adj_o <- colSums(A$mu[[1]] * xo_s)
        r <- A$lam[1] * (.as_base(crossprod(A$gcm, xo_s)) -
             matrix(adj_o, A$n_cells, k, byrow = TRUE))
        if (!is.null(A$excess[[1]])) {
            r <- r - A$lam[1] * .as_base(crossprod(A$excess[[1]], xo))
        }
        adj_h <- colSums(A$mu[[2]] * xh_s)
        ht <- .as_base(crossprod(A$gcm, xh_s))
        ht <- .as_base(crossprod(A$W, ht))
        ht <- ht - matrix(adj_h, A$n_cells, k, byrow = TRUE)
        r <- r + A$lam[2] * ht
        if (!is.null(A$excess[[2]])) {
            r <- r - A$lam[2] * .as_base(crossprod(A$excess[[2]], xh))
        }
    }
    r
}

# Standard path: materialize full BANKSY matrix as a Seurat assay
.banksy_standard <- function(object, data_own, knn_list, M, lambda, group,
                             split.scale, assay, slot, assay_name, verbose,
                             chunk_size, parallel, num_cores) {
    # Compute harmonics
    center <- rep(TRUE, length(M))
    center[1] <- FALSE
    har <- Map(function(knn_df, M, center) {
      x <- Banksy:::computeHarmonics(gcm=data_own,
                                     knn_df=knn_df,
                                     M=M, center=center,
                                     verbose=verbose,
                                     chunk_size=chunk_size,
                                     parallel=parallel,
                                     num_cores=num_cores)
      rownames(x) <- paste0(rownames(x), '.m', M)
      x
    }, knn_list, M, center)

    # Scale by lambdas
    lambdas <- Banksy:::getLambdas(lambda, n_harmonics = length(har))

    # Merge with own expression (coerce to dense for Seurat operations)
    if (verbose) message('Creating Banksy matrix')
    data_own <- as.matrix(data_own)
    data_banksy <- c(list(data_own), har)
    if (verbose) message('Scaling BANKSY matrix. Do not call ScaleData on assay ', assay_name)
    data_scaled <- lapply(data_banksy, fast_scaler,
                          object, group, split.scale, verbose)

    # Multiple by lambdas
    data_banksy <- Map(function(lam, mat) lam * mat, lambdas, data_banksy)
    data_scaled <- Map(function(lam, mat) lam * mat, lambdas, data_scaled)

    # Rbind
    data_banksy <- do.call(rbind, data_banksy)
    data_scaled <- do.call(rbind, data_scaled)

    # Create an assay object
    if (grepl(pattern = 'counts', x = slot)) {
        banksy_assay <- Seurat::CreateAssayObject(counts = data_banksy)
    } else {
        banksy_assay <- Seurat::CreateAssayObject(data = data_banksy)
    }

    # Add assay to Seurat object and set as default
    if (verbose) message('Setting default assay to ', assay_name)
    object[[assay_name]] <- banksy_assay
    DefaultAssay(object) <- assay_name
    object <- SetAssayData(object, layer = 'scale.data', new.data = data_scaled,
                           assay = assay_name)
    object
}

# Helpers
# Get own expression matrix from Seurat object
get_data <- function(object, assay, slot, features, verbose) {
    # Fetch data from Seurat
    if (verbose) message('Fetching data from slot ', slot,' from assay ', assay)
    data_own <- Seurat::GetAssayData(object = object, assay = assay, layer = slot)
    # Feature subset
    if (features[1] != 'all') {
        if (verbose) message('Subsetting by features')
        if (features[1] == 'variable') {
            feat <- Seurat::VariableFeatures(object)
            if (length(feat) == 0) {
                warning('No variable features found. Running Seurat::FindVariableFeatures')
                object <- Seurat::FindVariableFeatures(object)
                feat <- Seurat::VariableFeatures(object)
            }
        } else {
            feat <- features[which(rownames(object) %in% features)]
            if (length(feat) == 0) stop('None of the specified features found. Check if features in Seurat object')
        }
        data_own <- data_own[feat,,drop=FALSE]
    }
    return(data_own)
}

# Get locations from Seurat object
get_locs <- function(object, dimx, dimy, dimz, ndim, data_own, group, verbose) {

    if (!is.null(dimx) & !is.null(dimy)) {
        # Extract locations from metadata
        locs <- data.frame(
            sdimx = unlist(object[[dimx]]),
            sdimy = unlist(object[[dimy]])
        )
        rownames(locs) <- colnames(object)

        # Add z-dim if present
        if (!is.null(dimz)) locs$sdimz = object[[dimz]]

        # Check locations
        obj_samples <- colnames(data_own)
        locs_samples <- rownames(locs)
        if (any(is.na(match(obj_samples, locs_samples)))) {
            na_id <- which(is.na(match(obj_samples, locs_samples)))
            warning('No centroids found for samples: ',
                    paste(obj_samples[na_id], collapse = ', '), '. Dropping samples.')
            data_own <- data_own[, -na_id, drop = FALSE]
        }
        locs <- locs[match(obj_samples, locs_samples),,drop=FALSE]

    } else {
        # Extract locations with Seurat accessor
        locs <- Seurat::GetTissueCoordinates(object)[,seq_len(ndim)]
    }

    dim_names <- paste0('sdim', c('x','y','z'))
    colnames(locs) <- dim_names[seq_len(ncol(locs))]

    if (!is.null(group)) {
        # Stagger locations by group
        if (verbose) message('Staggering locations by ', group)
        locs[,1] = locs[,1] + abs(min(locs[,1]))
        max_x = max(locs[,1]) * 2
        n_groups = length(unique(unlist(object[[group]])))
        shift = seq(from = 0, length.out = n_groups, by = max_x)
        shift_order = match(unique(unlist(object[[group]])), names(table(object[[group]])))
        locs[, 1] = locs[, 1] + rep(shift, table(object[[group]])[shift_order])
    }

    return(locs)
}

# Scaling
fast_scaler = function(data, object, group, split.scale, verbose) {
    # Split scaling by group
    if (!is.null(group) & split.scale) {
        groups = unlist(object[[group]])
        ugroups = unique(groups)
        for (curr_group in ugroups) {
            if (verbose) message('Scaling group: ', curr_group)
            curr_group_id <- which(curr_group == groups)
            data[, curr_group_id] <- Seurat:::FastRowScale(
              data[, curr_group_id])
        }
    } else {
        data <- Seurat::FastRowScale(data)
    }
    data
}
