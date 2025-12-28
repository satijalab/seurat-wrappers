
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
#' @param verbose (boolean) Print messages
#'
#' @return A Seurat object with new assay holding a Banksy matrix
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
                      assay_name='BANKSY', M=NULL, verbose=TRUE) {
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

    # Create Banksy matrix
    M <- seq(0, max(Banksy:::getM(use_agf, M)))
    # Compute harmonics
    center <- rep(TRUE, length(M))
    # Only center higher harmonics
    center[1] <- FALSE
    har <- Map(function(knn_df, M, center) {
      x <- Banksy:::computeHarmonics(data_own, knn_df, M, center, verbose)
      rownames(x) <- paste0(rownames(x), '.m', M)
      x
    }, knn_list, M, center)

    # Scale by lambdas
    lambdas <- Banksy:::getLambdas(lambda, n_harmonics = length(har))

    # Merge with own expression
    if (verbose) message('Creating Banksy matrix')
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

    # Log commands
    object <- Seurat::LogSeuratCommand(object = object)

  return(object)
}

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
    data_own <- as.matrix(x = data_own)
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
        locs[,1] = locs[,1] + rep(shift, table(object[[group]]))
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
