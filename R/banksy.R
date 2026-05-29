
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
#' @param pca_backend (character) Backend for PCA computation when lazy=TRUE.
#'   'cpp' (default) uses C++ irlba for lower memory and faster runtime.
#'   'r' uses R's irlba package.
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
                      pca_backend=c("cpp", "r"),
                      verbose=TRUE) {
    # Check packages
    SeuratWrappers:::CheckPackage(package = 'data.table', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Matrix', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Banksy', repository = 'github')

    # Check lambda param
    if (lambda < 0 || lambda > 1) stop('Lambda must be between 0 and 1')

    # Get data
    data_own <- get_data(object, assay, slot, features, verbose)

    # Resolve harmonics
    M <- seq(0, max(Banksy:::getM(use_agf, M)))

    # Per-group lazy path: kNN per group on raw coordinates
    if (lazy && !is.null(group)) {
        locs <- get_locs(object, dimx, dimy, dimz, ndim, data_own,
                         group = NULL, verbose)
        stag <- get_locs(object, dimx, dimy, dimz, ndim, data_own,
                         group, verbose)
        object <- AddMetaData(object, metadata = stag,
                              col.name = paste0('staggered_', colnames(stag)))

        groups_vec <- unlist(object[[group]])
        ugroups <- unique(groups_vec)
        group_idx <- lapply(ugroups, function(g) which(g == groups_vec))

        if (verbose) {
            message('Computing per-group neighbors')
            for (gr in seq_along(group_idx))
                message('  ', ugroups[gr], ': ', length(group_idx[[gr]]), ' cells')
        }

        group_locs <- lapply(group_idx, function(cid) locs[cid, , drop = FALSE])

        if (!is.null(num_cores) && num_cores > 1L) {
            if (verbose) message('  Parallel kNN: mclapply (', num_cores, ' cores)')
            group_knn <- parallel::mclapply(group_locs,
                function(loc_slice) {
                    lapply(k_geom, function(kg)
                        Banksy:::computeNeighbors(
                            loc_slice, spatial_mode = spatial_mode,
                            k_geom = kg, n = n, sigma = sigma,
                            alpha = alpha, k_spatial = k_spatial,
                            verbose = FALSE))
                }, mc.cores = num_cores)
        } else {
            group_knn <- lapply(group_locs,
                function(loc_slice) {
                    lapply(k_geom, function(kg)
                        Banksy:::computeNeighbors(
                            loc_slice, spatial_mode = spatial_mode, k_geom = kg,
                            n = n, sigma = sigma, alpha = alpha,
                            k_spatial = k_spatial, verbose = FALSE))
                })
        }
        rm(group_locs)
        rm(locs, stag); gc(verbose = FALSE)

        object <- .banksy_lazy_pca(
            object, data_own, NULL, M, lambda, group, split.scale,
            assay, assay_name, npcs, scale_max, verbose,
            group_knn = group_knn, group_idx = group_idx,
            pca_backend = pca_backend)
        rm(data_own, group_knn)
    } else {
        locs <- get_locs(object, dimx, dimy, dimz, ndim, data_own, group, verbose)
        if (!is.null(group)) {
            object <- AddMetaData(
                object, metadata = locs,
                col.name = paste0('staggered_', colnames(locs)))
        }

        knn_list <- lapply(k_geom, function(kg) {
          Banksy:::computeNeighbors(locs,
                                    spatial_mode = spatial_mode, k_geom = kg,
                                    n = n, sigma = sigma, alpha = alpha,
                                    k_spatial = k_spatial, verbose = verbose)
        })

        if (!lazy && nrow(data_own) * ncol(data_own) > 1e8) {
            message('Note: dataset is large (', nrow(data_own), ' features x ',
                    ncol(data_own), ' cells). Consider using lazy=TRUE for ',
                    'memory-efficient PCA without materializing the full BANKSY matrix.')
        }

        if (lazy) {
            object <- .banksy_lazy_pca(
                object, data_own, knn_list, M, lambda, group, split.scale,
                assay, assay_name, npcs, scale_max, verbose,
                pca_backend = pca_backend)
            rm(data_own, knn_list)
        } else {
            object <- .banksy_standard(
                object, data_own, knn_list, M, lambda, group, split.scale,
                assay, slot, assay_name, verbose, chunk_size, parallel, num_cores)
        }
    }

    # Log commands
    object <- Seurat::LogSeuratCommand(object = object)
    return(object)
}

# Lazy PCA path
.banksy_lazy_pca <- function(object, data_own, knn_list, M, lambda, group,
                             split.scale, assay, assay_name, npcs, scale_max,
                             verbose,
                             group_knn = NULL, group_idx = NULL,
                             pca_backend = "cpp") {
    # Guards
    if (max(M) > 0) stop('lazy=TRUE currently only supports M=0 (no AGF)')
    SeuratWrappers:::CheckPackage(package = 'irlba', repository = 'CRAN')

    n_genes <- nrow(data_own)
    n_cells <- ncol(data_own)
    gene_names <- rownames(data_own)
    cell_names <- colnames(data_own)
    lambdas <- Banksy:::getLambdas(lambda, n_harmonics = 1)
    has_groups <- !is.null(group_knn)
    split_scale <- !is.null(group) && split.scale

    # Validate npcs
    max_npcs <- min(2L * n_genes, n_cells) - 1L
    if (npcs >= max_npcs) {
        stop('npcs (', npcs, ') must be < min(2*n_genes, n_cells) = ',
             max_npcs + 1L)
    }
    if (min(2L * n_genes, n_cells) < 6L) {
        stop('lazy=TRUE requires at least 3 genes and 6 cells')
    }

    if (has_groups) {
        # Per-group decomposition
        if (verbose) message('Building per-group weight matrices')
        W_list <- lapply(seq_along(group_idx), function(gr)
            Banksy:::.buildWeightMatrix(group_knn[[gr]][[1]],
                                        length(group_idx[[gr]])))
        rm(group_knn); gc(verbose = FALSE)

        if (verbose) message('Building per-group expression matrices')

        # Materialize per-group dgCMatrix
        gcm_list <- lapply(group_idx, function(cid)
            as(data_own[, cid], 'dgCMatrix'))

        groups <- if (split_scale) unlist(object[[group]]) else NULL
        ugroups <- if (split_scale) unique(groups) else NULL
        own_result <- .lazy_own_scaling(NULL, split_scale, group_idx,
                                        groups, ugroups,
                                        n_genes, n_cells, scale_max, verbose,
                                        gcm_list = gcm_list)

        h0_result <- .lazy_h0_scaling_grouped(
            NULL, W_list, split_scale, group_idx,
            n_genes, n_cells, scale_max, verbose,
            gcm_list = gcm_list)

        if (verbose) {
            n_own_clip <- if (!is.null(own_result$excess)) length(own_result$excess@x) else 0
            n_h0_clip <- if (!is.null(h0_result$excess)) length(h0_result$excess@x) else 0
            message('Clipping corrections: own=', n_own_clip,
                    ' H0=', n_h0_clip, ' entries')
        }

        # Extract W CSC slots
        w_csc <- lapply(W_list, function(w)
            if (is(w, 'dgCMatrix')) w else as(w, 'dgCMatrix'))
        w_slots <- list(
            i = lapply(w_csc, slot, 'i'),
            p = lapply(w_csc, slot, 'p'),
            x = lapply(w_csc, slot, 'x'),
            ncol = vapply(w_csc, ncol, integer(1)))
        rm(w_csc)

        gcm_slots <- list(
            i = lapply(gcm_list, slot, 'i'),
            p = lapply(gcm_list, slot, 'p'),
            x = lapply(gcm_list, slot, 'x'))
        rm(gcm_list, data_own); gc(verbose = FALSE)

        banksy_op <- structure(list(
            gcm_slots = gcm_slots,
            W_list = W_list, w_slots = w_slots,
            mu = list(own_result$mu, h0_result$mu),
            sd = list(own_result$sd, h0_result$sd),
            lam = lambdas,
            excess = list(own_result$excess, h0_result$excess),
            valid = list(own_result$valid, h0_result$valid),
            has_groups = TRUE,
            split_scale = split_scale,
            group_idx = group_idx,
            n_genes = n_genes, n_cells = n_cells
        ), class = 'BanksyLazy')
        rm(own_result, h0_result)
    } else {
        # Single-matrix path
        if (!is.null(group) && !split.scale) {
            if (verbose) warning('lazy=TRUE: group is used for spatial ',
                                 'staggering only; scaling is performed globally')
        }

        if (split_scale) {
            groups <- unlist(object[[group]])
            ugroups <- unique(groups)
            group_idx <- lapply(ugroups, function(g) which(g == groups))
        }

        if (verbose) message('Building sparse weight matrix')
        W <- Banksy:::.buildWeightMatrix(knn_list[[1]], n_cells)
        rm(knn_list); gc(verbose = FALSE)  # kNN data consumed

        own_result <- .lazy_own_scaling(data_own, split_scale, group_idx,
                                        groups = if (split_scale) groups else NULL,
                                        ugroups = if (split_scale) ugroups else NULL,
                                        n_genes, n_cells, scale_max, verbose)
        h0_result <- .lazy_h0_scaling(data_own, W, split_scale, group_idx,
                                      n_genes, n_cells, scale_max, verbose)

        if (verbose) {
            n_own_clip <- if (!is.null(own_result$excess)) length(own_result$excess@x) else 0
            n_h0_clip <- if (!is.null(h0_result$excess)) length(h0_result$excess@x) else 0
            message('Clipping corrections: own=', n_own_clip,
                    ' H0=', n_h0_clip, ' entries')
        }

        banksy_op <- structure(list(
            gcm = data_own, W = W,
            mu = list(own_result$mu, h0_result$mu),
            sd = list(own_result$sd, h0_result$sd),
            lam = lambdas,
            excess = list(own_result$excess, h0_result$excess),
            valid = list(own_result$valid, h0_result$valid),
            has_groups = FALSE,
            split_scale = split_scale,
            group_idx = group_idx,
            n_genes = n_genes, n_cells = n_cells
        ), class = 'BanksyLazy')
        rm(own_result, h0_result)
    }

    # Reclaim memory before iterative solve
    gc(verbose = FALSE)

    # SVD solve
    pca_backend <- match.arg(pca_backend, c("cpp", "r"))
    if (pca_backend == "cpp") {
        # irlba algorithm with C++ Q buffer management.
        # Ported from irlba::irlba R source. The implicit restart rotates
        # V via lanczos_extract (C++), avoiding the memory spike from
        # R's copy-on-modify on V %*% Bsvd$v.
        # Peak memory = Q (n_cells * work) + extraction buffer (n_cells * k).
        work <- as.integer(min(npcs + 7L, n_cells - 1L))
        maxit <- 1000L
        tol <- 1e-5
        svtol <- tol
        if (verbose) message('Computing BANKSY PCA (', npcs, ' PCs) via ',
                             'C++ irlba (work=', work, ')')

        m <- 2L * n_genes
        n <- n_cells
        k <- npcs
        eps <- .Machine$double.eps
        eps23 <- eps^(2/3)
        sqrteps <- sqrt(eps)

        # C++ managed V buffer (n x work) — never copied by R
        V_buf <- numeric(as.double(n) * work)
        W <- matrix(0.0, m, work)
        F_vec <- numeric(n)
        B <- NULL

        set.seed(42)
        q1 <- rnorm(n); q1 <- q1 / lanczos_norm(q1)
        lanczos_set_col(V_buf, 0L, q1, n)
        rm(q1)

        t0 <- proc.time()["elapsed"]
        mprod <- 0L
        iter <- 0L
        Smax <- 1.0
        Smin <- NULL
        lastsv <- numeric(0)
        converged <- FALSE
        restart <- FALSE

        while (iter < maxit) {
            iter <- iter + 1L
            j <- if (iter == 1L && !restart) 1L else k + 1L

            # ── Lanczos bidiagonalization: steps j..work ──
            # First step of each cycle
            VJ <- lanczos_get_col(V_buf, j - 1L, n)
            avj <- as.numeric(.banksy_lazy_mult(
                banksy_op, matrix(VJ, ncol = 1)))
            rm(VJ)
            W[, j] <- avj
            mprod <- mprod + 1L

            if (iter > 1L && j > 1L)
                W[, j] <- W[, j] - as.numeric(
                    W[, 1:(j-1), drop=FALSE] %*%
                    crossprod(W[, 1:(j-1), drop=FALSE], W[, j]))

            S <- sqrt(sum(W[, j]^2))
            if (S < eps23) {
                W[, j] <- rnorm(m)
                if (j > 1L) W[, j] <- W[, j] - as.numeric(
                    W[, 1:(j-1), drop=FALSE] %*%
                    crossprod(W[, 1:(j-1), drop=FALSE], W[, j]))
                W[, j] <- W[, j] / sqrt(sum(W[, j]^2))
                S <- 0
            } else {
                W[, j] <- W[, j] / S
            }

            # Inner Lanczos loop
            while (j <= work) {
                # Adjoint
                F_vec <- as.numeric(.banksy_lazy_mult(
                    banksy_op, matrix(W[, j], ncol = 1), transpose = TRUE))
                mprod <- mprod + 1L
                F_vec <- F_vec - S * lanczos_get_col(V_buf, j - 1L, n)
                F_vec <- lanczos_reorth(V_buf, j, F_vec, n)

                if (j + 1L <= work) {
                    R <- lanczos_norm(F_vec)
                    if (R < eps23) {
                        F_vec <- rnorm(n)
                        F_vec <- lanczos_reorth(V_buf, j, F_vec, n)
                        lanczos_set_col(V_buf, j,
                            F_vec / lanczos_norm(F_vec), n)
                        R <- 0
                    } else {
                        lanczos_set_col(V_buf, j, F_vec / R, n)
                    }

                    # Build B
                    if (is.null(B)) {
                        B <- matrix(c(S, R), nrow = 1)
                    } else {
                        B <- rbind(cbind(B, 0),
                                   c(rep(0, ncol(B) - 1), S, R))
                    }

                    # Forward for next step
                    VJP1 <- lanczos_get_col(V_buf, j, n)
                    W[, j + 1L] <- as.numeric(.banksy_lazy_mult(
                        banksy_op, matrix(VJP1, ncol = 1)))
                    rm(VJP1)
                    mprod <- mprod + 1L
                    W[, j + 1L] <- W[, j + 1L] - R * W[, j]
                    if (j > 1L)
                        W[, j + 1L] <- W[, j + 1L] - as.numeric(
                            W[, 1:j, drop=FALSE] %*%
                            crossprod(W[, 1:j, drop=FALSE], W[, j + 1L]))
                    S <- sqrt(sum(W[, j + 1L]^2))
                    if (S < eps23) {
                        W[, j + 1L] <- rnorm(m)
                        if (j > 0L) W[, j + 1L] <- W[, j + 1L] - as.numeric(
                            W[, 1:j, drop=FALSE] %*%
                            crossprod(W[, 1:j, drop=FALSE], W[, j + 1L]))
                        W[, j + 1L] <- W[, j + 1L] / sqrt(sum(W[, j + 1L]^2))
                        S <- 0
                    } else {
                        W[, j + 1L] <- W[, j + 1L] / S
                    }
                } else {
                    B <- rbind(B, c(rep(0, j - 1), S))
                }
                j <- j + 1L
            }

            # ── Convergence test ──
            Bsz <- nrow(B)
            R_F <- lanczos_norm(F_vec)
            F_vec <- F_vec / R_F
            Bsvd <- svd(B)
            Smax <- max(Smax, Bsvd$d[1])

            R_vals <- R_F * Bsvd$u[Bsz, , drop = FALSE]
            ct_conv <- all(abs(R_vals[1:k]) < tol * Smax)
            sv_conv <- length(lastsv) >= k &&
                all(abs(Bsvd$d[1:k] - lastsv[1:k]) < svtol * Bsvd$d[1:k])
            lastsv <- Bsvd$d

            if (verbose && (iter <= 2L || iter %% 5 == 0 || ct_conv || sv_conv)) {
                rss <- tryCatch({
                    l <- readLines("/proc/self/status", warn = FALSE)
                    v <- grep("^VmRSS:", l, value = TRUE)
                    as.numeric(gsub("[^0-9]", "", v)) / 1024^2
                }, error = function(e) NA)
                elapsed <- proc.time()["elapsed"] - t0
                message(sprintf(
                    "  iter=%d  mprod=%d  sv[%d]=%.4e  RSS=%.1fGB  t=%.0fs",
                    iter, mprod, k, Bsvd$d[k], rss, elapsed))
            }

            if (ct_conv || sv_conv) {
                converged <- TRUE
                if (verbose) message(sprintf(
                    "  Converged: iter=%d, mprod=%d", iter, mprod))
                break
            }

            # ── Implicit restart (from irlba source lines 434-438) ──
            # V[,1:(k+1)] = cbind(V %*% Bsvd$v[,1:k], F_vec)
            # In R's irlba this is V[,...] <- cbind(V %*% Bsvd$v, F)
            # which creates a full copy of V. We use C++ extraction instead.
            V_rot <- lanczos_extract(V_buf, Bsvd$v[, 1:k, drop = FALSE],
                                     n, Bsz, k)
            for (i in seq_len(k))
                lanczos_set_col(V_buf, i - 1L, V_rot[, i], n)
            rm(V_rot)
            lanczos_set_col(V_buf, k, F_vec, n)

            # B = [diag(d[1:k]) | R[1:k]]
            B <- cbind(diag(Bsvd$d[1:k], nrow = k),
                       R_F * Bsvd$u[Bsz, 1:k])

            # W[,1:k] = W %*% Bsvd$u[,1:k] (tiny, m=10202)
            W[, 1:k] <- W[, 1:Bsz, drop = FALSE] %*%
                Bsvd$u[, 1:k, drop = FALSE]

            restart <- TRUE
        }

        if (!converged)
            warning("C++ irlba did not converge; try increasing maxit or work")

        # Final extraction (C++ managed, no spike)
        Bsvd <- svd(B)
        V <- lanczos_extract(V_buf, Bsvd$v[, 1:npcs, drop = FALSE],
                             n, nrow(B), npcs)
        rm(V_buf); gc(verbose = FALSE)

        U <- W[, 1:nrow(B), drop = FALSE] %*%
             Bsvd$u[, 1:npcs, drop = FALSE]
        rm(W)

        d <- Bsvd$d[1:npcs]
        rm(Bsvd, B, F_vec)
        pca <- list(u = U, v = V, d = d)

    } else {
        if (verbose) message('Computing BANKSY PCA (', npcs, ' PCs) via irlba')
        SeuratWrappers:::CheckPackage(package = 'irlba', repository = 'CRAN')

        .irlba_iter <- 0L
        .irlba_t0 <- proc.time()["elapsed"]
        .instrumented_mult <- function(A, x, transpose = FALSE) {
            .irlba_iter <<- .irlba_iter + 1L
            if (.irlba_iter %% 20 == 0) gc(verbose = FALSE)
            if (.irlba_iter %% 20 == 1) {
                rss <- tryCatch({
                    l <- readLines("/proc/self/status", warn = FALSE)
                    v <- grep("^VmRSS:", l, value = TRUE)
                    as.numeric(gsub("[^0-9]", "", v)) / 1024^2
                }, error = function(e) NA)
                elapsed <- proc.time()["elapsed"] - .irlba_t0
                message(sprintf("  irlba iter=%d  t=%.0fs  RSS=%.1fGB",
                                .irlba_iter, elapsed, rss))
            }
            .banksy_lazy_mult(A, x, transpose)
        }
        pca <- irlba::irlba(banksy_op, nv = npcs, mult = .instrumented_mult)
    }

    # Cell embeddings: V * D
    embeddings <- sweep(pca$v, 2, pca$d, `*`)
    rownames(embeddings) <- cell_names
    colnames(embeddings) <- paste0(assay_name, '_', seq_len(npcs))

    # Feature loadings
    loadings <- pca$u
    feat_names <- c(gene_names, paste0(gene_names, '.m0'))
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
                              ugroups, n_genes, n_cells, scale_max, verbose,
                              gcm_list = NULL) {
    if (verbose) message('Computing scaling parameters for own expression')

    if (split_scale) {
        mu <- matrix(0, nrow = n_genes, ncol = length(group_idx))
        sd <- matrix(0, nrow = n_genes, ncol = length(group_idx))
        for (gr in seq_along(group_idx)) {
            grp <- if (!is.null(gcm_list)) gcm_list[[gr]]
                   else data_own[, group_idx[[gr]], drop = FALSE]
            n_c <- as.double(ncol(grp))
            mu[, gr] <- Matrix::rowMeans(grp)
            sd[, gr] <- sqrt(pmax(
                n_c / (n_c - 1) * (Matrix::rowMeans(grp^2) - mu[, gr]^2), 0
            ))
        }
        valid <- rowSums(sd == 0) == 0
        sd[sd == 0] <- 1
    } else {
        if (!is.null(gcm_list)) {
            # Aggregate global stats from per-group data
            row_sums <- numeric(n_genes)
            row_sq_sums <- numeric(n_genes)
            for (gr in seq_along(gcm_list)) {
                row_sums <- row_sums + Matrix::rowSums(gcm_list[[gr]])
                row_sq_sums <- row_sq_sums + Matrix::rowSums(gcm_list[[gr]]^2)
            }
            n_c <- as.double(n_cells)
            mu <- row_sums / n_c
            sd <- sqrt(pmax(n_c / (n_c - 1) * (row_sq_sums / n_c - mu^2), 0))
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
        }
        sd[sd == 0] <- 1
        valid <- NULL
    }

    # Compute clipping excess via C++ (no R temporaries)
    if (verbose) message('Computing clipping excess for own expression')
    if (!is.null(gcm_list)) {
        mu_mat <- if (is.matrix(mu)) mu else matrix(0, 0, 0)
        sd_mat <- if (is.matrix(sd)) sd else matrix(0, 0, 0)
        mu_vec <- if (!is.matrix(mu)) mu else numeric(0)
        sd_vec <- if (!is.matrix(sd)) sd else numeric(0)
        exc <- own_excess_cpp(
            lapply(gcm_list, slot, 'i'), lapply(gcm_list, slot, 'p'),
            lapply(gcm_list, slot, 'x'), group_idx,
            split_scale, mu_mat, sd_mat, mu_vec, sd_vec,
            valid, scale_max, n_genes, n_cells)
        excess <- if (length(exc$i) > 0L)
            Matrix::sparseMatrix(i = exc$i, j = exc$j, x = exc$x,
                                 dims = c(n_genes, n_cells))
        else NULL
    } else {
        excess <- .lazy_own_excess(data_own, mu, sd, valid, split_scale,
                                   group_idx, groups, ugroups,
                                   n_genes, n_cells, scale_max)
    }

    list(mu = mu, sd = sd, valid = valid, excess = excess)
}

.lazy_own_excess <- function(data_own, mu, sd, valid, split_scale, group_idx,
                             groups, ugroups, n_genes, n_cells, scale_max,
                             gcm_list = NULL) {
    if (!is.null(gcm_list)) {
        # Per-group excess from materialized dgCMatrix list
        exc_i <- integer(0); exc_j <- integer(0); exc_x <- numeric(0)
        for (gr in seq_along(gcm_list)) {
            g <- gcm_list[[gr]]
            cid <- group_idx[[gr]]
            gi <- g@i + 1L
            gj <- rep(cid, diff(g@p))
            if (split_scale) {
                g_mu <- mu[cbind(gi, gr)]
                g_sd <- sd[cbind(gi, gr)]
            } else {
                g_mu <- mu[gi]
                g_sd <- sd[gi]
            }
            exceed <- g@x > g_mu + scale_max * g_sd
            if (split_scale) exceed <- exceed & valid[gi]
            if (any(exceed)) {
                z_vals <- (g@x[exceed] - g_mu[exceed]) / g_sd[exceed]
                exc_i <- c(exc_i, gi[exceed])
                exc_j <- c(exc_j, gj[exceed])
                exc_x <- c(exc_x, z_vals - scale_max)
            }
        }
        if (length(exc_i) > 0)
            return(Matrix::sparseMatrix(i = exc_i, j = exc_j, x = exc_x,
                                        dims = c(n_genes, n_cells)))
        else
            return(NULL)
    }
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

    if (exc_n == 0L) return(NULL)
    Matrix::sparseMatrix(
        i = exc_i[1:exc_n], j = exc_j[1:exc_n],
        x = exc_x[1:exc_n], dims = c(n_genes, n_cells))
}

# Per-group H0 scaling via C++ column-wise sweep (no intermediate allocation)
.lazy_h0_scaling_grouped <- function(data_own, W_list, split_scale, group_idx,
                                     n_genes, n_cells, scale_max, verbose,
                                     gcm_list = NULL) {
    if (verbose) message('Computing scaling params and clipping for H0 (per-group)')
    n_groups <- length(group_idx)

    # Extract gcm CSC slots (always per-group)
    gcm_i <- lapply(gcm_list, slot, 'i')
    gcm_p <- lapply(gcm_list, slot, 'p')
    gcm_x <- lapply(gcm_list, slot, 'x')

    # Extract W CSC slots
    w_csc <- lapply(W_list, function(w)
        if (is(w, 'dgCMatrix')) w else as(w, 'dgCMatrix'))
    w_i <- lapply(w_csc, slot, 'i')
    w_p <- lapply(w_csc, slot, 'p')
    w_x <- lapply(w_csc, slot, 'x')
    rm(w_csc)

    if (verbose) {
        for (gr in seq_along(group_idx))
            message('  Group ', gr, ': ', length(group_idx[[gr]]), ' cells')
    }

    # Pass 1: per-group stats via C++
    stats <- h0_group_stats_cpp(
        gcm_i, gcm_p, gcm_x, nrow(gcm_list[[1]]),
        w_i, w_p, w_x, group_idx)

    mu_g <- stats$mu
    ss_g <- stats$ss
    max_g <- stats$max_val
    rm(stats)

    # Derive mu, sd, valid
    if (split_scale) {
        mu <- mu_g
        sd <- matrix(0, nrow = n_genes, ncol = n_groups)
        for (gr in seq_along(group_idx)) {
            n_c <- as.double(length(group_idx[[gr]]))
            sd[, gr] <- sqrt(pmax(
                n_c / (n_c - 1) * (ss_g[, gr] / n_c - mu[, gr]^2), 0
            ))
        }
        valid <- rowSums(sd == 0) == 0
        sd[sd == 0] <- 1
    } else {
        n_g_vec <- vapply(group_idx, length, integer(1))
        mu <- rowSums(sweep(mu_g, 2, n_g_vec, `*`)) / n_cells
        ss <- rowSums(ss_g)
        n_c <- as.double(n_cells)
        sd <- sqrt(pmax(n_c / (n_c - 1) * (ss / n_c - mu * mu), 0))
        sd[sd == 0] <- 1
        valid <- NULL
    }

    # Identify clip genes
    if (split_scale) {
        thresh <- mu + scale_max * sd
        clip_genes <- which(valid & rowSums(max_g > thresh) > 0)
    } else {
        thresh <- mu + scale_max * sd
        max_h0 <- apply(max_g, 1, max)
        clip_genes <- which(max_h0 > thresh)
    }
    rm(mu_g, ss_g, max_g)
    if (verbose) message('H0 genes requiring clipping: ', length(clip_genes),
                         ' / ', n_genes)

    # Pass 2: excess via C++ (only for clip genes, no dense intermediates)
    if (length(clip_genes) == 0) {
        excess <- NULL
    } else {
        mu_mat <- if (is.matrix(mu)) mu else matrix(mu, ncol = 1)
        sd_mat <- if (is.matrix(sd)) sd else matrix(sd, ncol = 1)
        exc <- h0_group_excess_cpp(
            gcm_i, gcm_p, gcm_x, nrow(gcm_list[[1]]),
            w_i, w_p, w_x, group_idx,
            as.integer(clip_genes), split_scale, mu_mat, sd_mat, scale_max)
        excess <- if (length(exc$i) > 0L)
            Matrix::sparseMatrix(i = exc$i, j = exc$j, x = exc$x,
                                 dims = c(n_genes, n_cells))
        else NULL
    }

    list(mu = mu, sd = sd, valid = valid, excess = excess)
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

# ─── Lazy operator dispatch ──────────────────────────────────────────────────

.banksy_lazy_mult <- function(A, x, transpose = FALSE) {
    if (inherits(x, 'BanksyLazy')) {
        tmp <- A; A <- x; x <- tmp; transpose <- TRUE
    }
    x <- .as_base(x)
    k <- ncol(x)
    if (!transpose) .banksy_forward(A, x, k)
    else .banksy_adjoint(A, x, k)
}

.banksy_forward <- function(A, x, k) {
    if (isTRUE(A$has_groups)) {
        ws <- A$w_slots
        gs_i <- A$gcm_slots$i
        gs_p <- A$gcm_slots$p
        gs_x <- A$gcm_slots$x

        if (A$split_scale) {
            result <- banksy_forward_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx, x,
                TRUE, A$mu[[1]], A$sd[[1]], A$mu[[2]], A$sd[[2]],
                A$lam, A$valid[[1]], A$valid[[2]])

            if (!is.null(A$excess[[1]]) || !is.null(A$excess[[2]])) {
                ng <- A$n_genes
                own <- result[1:ng, , drop = FALSE]
                h0  <- result[(ng+1):(2*ng), , drop = FALSE]
                if (!is.null(A$excess[[1]]))
                    own <- own - A$lam[1] * .as_base(A$excess[[1]] %*% x)
                if (!is.null(A$excess[[2]]))
                    h0 <- h0 - A$lam[2] * .as_base(A$excess[[2]] %*% x)
                result <- rbind(own, h0)
            }
        } else {
            result <- banksy_forward_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx, x,
                FALSE,
                matrix(0, 0, 0), matrix(0, 0, 0),
                matrix(0, 0, 0), matrix(0, 0, 0),
                A$lam, NULL, NULL)

            ng <- A$n_genes
            own <- result[1:ng, , drop = FALSE]
            h0  <- result[(ng+1):(2*ng), , drop = FALSE]
            cs <- colSums(x)
            own <- A$lam[1] * (own - outer(A$mu[[1]], cs)) / A$sd[[1]]
            h0  <- A$lam[2] * (h0  - outer(A$mu[[2]], cs)) / A$sd[[2]]
            if (!is.null(A$excess[[1]]))
                own <- own - A$lam[1] * .as_base(A$excess[[1]] %*% x)
            if (!is.null(A$excess[[2]]))
                h0 <- h0 - A$lam[2] * .as_base(A$excess[[2]] %*% x)
            result <- rbind(own, h0)
        }
        result
    } else {
        cs <- colSums(x)
        own <- .as_base(A$gcm %*% x)
        own <- A$lam[1] * (own - outer(A$mu[[1]], cs)) / A$sd[[1]]
        if (!is.null(A$excess[[1]]))
            own <- own - A$lam[1] * .as_base(A$excess[[1]] %*% x)
        Wx <- .as_base(A$W %*% x)
        h0 <- .as_base(A$gcm %*% Wx)
        h0 <- A$lam[2] * (h0 - outer(A$mu[[2]], cs)) / A$sd[[2]]
        if (!is.null(A$excess[[2]]))
            h0 <- h0 - A$lam[2] * .as_base(A$excess[[2]] %*% x)
        rbind(own, h0)
    }
}

.banksy_adjoint <- function(A, x, k) {
    ng <- A$n_genes
    xo <- x[1:ng, , drop = FALSE]
    xh <- x[(ng+1):(2*ng), , drop = FALSE]

    if (isTRUE(A$has_groups)) {
        ws <- A$w_slots
        gs_i <- A$gcm_slots$i
        gs_p <- A$gcm_slots$p
        gs_x <- A$gcm_slots$x

        if (A$split_scale) {
            if (!is.null(A$valid[[1]])) xo[!A$valid[[1]], ] <- 0
            if (!is.null(A$valid[[2]])) xh[!A$valid[[2]], ] <- 0

            r <- banksy_adjoint_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx,
                xo, xh,
                numeric(k), numeric(k),
                TRUE, A$mu[[1]], A$sd[[1]], A$mu[[2]], A$sd[[2]],
                A$lam)

            if (!is.null(A$excess[[1]]))
                r <- r - A$lam[1] * .as_base(crossprod(A$excess[[1]], xo))
            if (!is.null(A$excess[[2]]))
                r <- r - A$lam[2] * .as_base(crossprod(A$excess[[2]], xh))
        } else {
            xo_s <- xo / A$sd[[1]]
            xh_s <- xh / A$sd[[2]]
            adj_o <- colSums(A$mu[[1]] * xo_s)
            adj_h <- colSums(A$mu[[2]] * xh_s)

            r <- banksy_adjoint_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx,
                xo_s, xh_s, adj_o, adj_h,
                FALSE,
                matrix(0, 0, 0), matrix(0, 0, 0),
                matrix(0, 0, 0), matrix(0, 0, 0),
                A$lam)

            if (!is.null(A$excess[[1]]))
                r <- r - A$lam[1] * .as_base(crossprod(A$excess[[1]], xo))
            if (!is.null(A$excess[[2]]))
                r <- r - A$lam[2] * .as_base(crossprod(A$excess[[2]], xh))
        }
        r
    } else {
        xo_s <- xo / A$sd[[1]]
        xh_s <- xh / A$sd[[2]]
        adj_o <- colSums(A$mu[[1]] * xo_s)
        r <- A$lam[1] * (.as_base(crossprod(A$gcm, xo_s)) -
             matrix(adj_o, A$n_cells, k, byrow = TRUE))
        if (!is.null(A$excess[[1]]))
            r <- r - A$lam[1] * .as_base(crossprod(A$excess[[1]], xo))
        adj_h <- colSums(A$mu[[2]] * xh_s)
        ht <- .as_base(crossprod(A$gcm, xh_s))
        ht <- .as_base(crossprod(A$W, ht))
        ht <- ht - matrix(adj_h, A$n_cells, k, byrow = TRUE)
        r <- r + A$lam[2] * ht
        if (!is.null(A$excess[[2]]))
            r <- r - A$lam[2] * .as_base(crossprod(A$excess[[2]], xh))
        r
    }
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
        # Skip row subset if all features selected in original order (avoids copy)
        if (length(feat) < nrow(data_own) ||
            !identical(feat, rownames(data_own))) {
            data_own <- data_own[feat, , drop = FALSE]
        }
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
        groups_vec <- unlist(object[[group]])
        ugroups <- unique(groups_vec)
        shift <- setNames(
            seq(from = 0, length.out = length(ugroups), by = max_x),
            ugroups)
        locs[, 1] = locs[, 1] + shift[groups_vec]
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
