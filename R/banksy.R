
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
#' @param assay_name (character) Name for Banksy assay in Seurat object
#' @param M (numeric) Advanced usage. Highest azimuthal harmonic
#' @param chunk_size A integer scalar specifying the number of rows / genes of
#'   the neighborhood cell matrix to compute. Must be less than floor of
#'   2e31-1 / number of cells. This is automatically computed but can be
#'   specified.
#' @param parallel A logical scalar specifying whether to compute chunks in
#'   parallel using bplapply. Not implemented for Windows.
#' @param num_cores A integer scalar specifying the number of cores to use
#'   if parallel is TRUE.
#' @param compute_pca (boolean) If TRUE, compute PCA directly without materializing
#'   the full BANKSY matrix. Enables analysis of very large datasets (millions
#'   of cells) that would otherwise exceed available memory. Currently only
#'   supported for M=0 (no AGF). Requires the irlba package.
#' @param npcs (integer) Number of PCs to compute when compute_pca=TRUE (default 50)
#' @param verbose (boolean) Print messages
#'
#' @return A Seurat object. If compute_pca=FALSE, contains a new assay holding the
#'   BANKSY matrix. If compute_pca=TRUE, contains a dimensionality reduction named
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
                      compute_pca=FALSE, npcs=50L,
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

    if (compute_pca) {
        # -----------------------------------------------------------------
        # Lazy PCA path: compute PCA directly without forming the full
        # BANKSY matrix. Peak memory ~ sparse input + sparse W + vectors.
        # Currently M=0 only (no AGF).
        # -----------------------------------------------------------------
        if (max(M) > 0) stop('compute_pca=TRUE currently only supports M=0 (no AGF)')
        if (!is.null(group)) {
            if (split.scale) {
                stop('compute_pca=TRUE does not yet support split.scale with group')
            }
            if (verbose) warning('compute_pca=TRUE: group is used for spatial ',
                                 'staggering only; scaling is performed globally')
        }

        SeuratWrappers:::CheckPackage(package = 'irlba', repository = 'CRAN')

        n_genes <- nrow(data_own)
        n_cells <- ncol(data_own)
        lambdas <- Banksy:::getLambdas(lambda, n_harmonics = 1)

        # Validate npcs
        max_npcs <- min(2L * n_genes, n_cells) - 1L
        if (npcs >= max_npcs) {
            stop('npcs (', npcs, ') must be < min(2*n_genes, n_cells) = ',
                 max_npcs + 1L)
        }
        if (min(2L * n_genes, n_cells) < 6L) {
            stop('compute_pca=TRUE requires at least 3 genes and 6 cells')
        }

        # Build sparse weight matrix
        if (verbose) message('Building sparse weight matrix')
        W <- Banksy:::.buildWeightMatrix(knn_list[[1]], n_cells)

        # Compute scaling parameters for data_own (sample sd, n-1 denom,
        # to match Seurat::FastRowScale)
        if (verbose) message('Computing scaling parameters for own expression')
        n_c <- as.double(n_cells)
        if (inherits(data_own, 'sparseMatrix')) {
            mu_own <- Matrix::rowMeans(data_own)
            sd_own <- sqrt(pmax(
                n_c / (n_c - 1) * (Matrix::rowMeans(data_own^2) - mu_own^2), 0
            ))
        } else {
            mu_own <- rowMeans(data_own)
            sd_own <- sqrt(pmax(
                n_c / (n_c - 1) * (rowMeans(data_own^2) - mu_own^2), 0
            ))
        }
        sd_own[sd_own == 0] <- 1

        # Compute scaling parameters for H0 = data_own %*% W (chunk-by-chunk)
        if (verbose) message('Computing scaling parameters for H0')
        h0_params <- Banksy:::.computeH0ScalingParams(data_own, W)

        # Define lazy linear operator for irlba.
        # A = rbind(lam0 * S(data_own), lam1 * S(H0))
        # where S(X) = row-wise z-score.
        # Forward:  A %*% v  (features x 1)
        # Adjoint:  t(A) %*% u  (cells x 1)
        banksy_op <- structure(list(
            gcm = data_own, W = W,
            mu = list(mu_own, h0_params$mu),
            sd = list(sd_own, h0_params$sd),
            lam = lambdas,
            n_genes = n_genes, n_cells = n_cells
        ), class = 'BanksyLazy')

        # Ensure base numeric matrix (no Matrix S4 classes for irlba compat)
        .as_base <- function(x) {
            if (inherits(x, 'Matrix')) x <- as.matrix(x)
            if (!is.matrix(x)) x <- as.matrix(x)
            storage.mode(x) <- 'double'
            x
        }

        banksy_mult <- function(A, x, transpose = FALSE) {
            # irlba may swap A and x for transpose multiply
            if (inherits(x, 'BanksyLazy')) {
                tmp <- A; A <- x; x <- tmp; transpose <- TRUE
            }
            x <- .as_base(x)
            k <- ncol(x)

            if (!transpose) {
                # A %*% x: x is n_cells x k
                cs <- colSums(x)
                own <- .as_base(A$gcm %*% x)
                own <- A$lam[1] * (own - outer(A$mu[[1]], cs)) / A$sd[[1]]
                Wx <- .as_base(A$W %*% x)
                h0 <- .as_base(A$gcm %*% Wx)
                h0 <- A$lam[2] * (h0 - outer(A$mu[[2]], cs)) / A$sd[[2]]
                rbind(own, h0)
            } else {
                # t(A) %*% x: x is (2*n_genes) x k
                ng <- A$n_genes
                xo <- x[1:ng, , drop = FALSE]
                xh <- x[(ng+1):(2*ng), , drop = FALSE]
                xo_s <- xo / A$sd[[1]]
                xh_s <- xh / A$sd[[2]]
                adj_o <- colSums(A$mu[[1]] * xo_s)
                r <- A$lam[1] * (.as_base(crossprod(A$gcm, xo_s)) -
                     matrix(adj_o, A$n_cells, k, byrow = TRUE))
                adj_h <- colSums(A$mu[[2]] * xh_s)
                ht <- .as_base(crossprod(A$gcm, xh_s))
                ht <- .as_base(crossprod(A$W, ht))
                ht <- ht - matrix(adj_h, A$n_cells, k, byrow = TRUE)
                r + A$lam[2] * ht
            }
        }

        if (verbose) message('Computing BANKSY PCA (', npcs, ' PCs) via lazy operator')
        pca <- irlba::irlba(banksy_op, nv = npcs, mult = banksy_mult)

        # Cell embeddings: V * D
        embeddings <- sweep(pca$v, 2, pca$d, `*`)
        rownames(embeddings) <- colnames(data_own)
        colnames(embeddings) <- paste0(assay_name, '_', seq_len(npcs))

        # Feature loadings
        loadings <- pca$u
        feat_names <- c(rownames(data_own),
                        paste0(rownames(data_own), '.m0'))
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
    } else {
        # -----------------------------------------------------------------
        # Standard path: materialize full BANKSY matrix as a Seurat assay.
        # -----------------------------------------------------------------

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
    }

    # Log commands
    object <- Seurat::LogSeuratCommand(object = object)

  return(object)
}

# S3 method for irlba to get dimensions of lazy BANKSY operator
#' @export
dim.BanksyLazy <- function(x) c(x$n_genes * 2L, x$n_cells)

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
