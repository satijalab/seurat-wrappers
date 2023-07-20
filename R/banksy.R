
#' @include internal.R
#'
NULL

#' Run Banksy on a Seurat Object
#'
#' @param object A Seurat object
#' @param lambda (numeric) Spatial weight parameter
#' @param assay (character) Assay in Seurat object to use
#' @param slot (character) Slot in Seurat assay to use
#' @param dimx (character) Column name of spatial x dimension (must be in metadata)
#' @param dimy (character) Column name of spatial y dimension (must be in metadata)
#' @param features (character) Features to compute. Can be 'all', 'variable' or
#'   a vector of feature names
#' @param M (numeric) highest azimuthal harmonic
#' @param k_geom (numeric) kNN parameter - number of neighbors to use
#' @param n (numeric) kNN_rn parameter - exponent of radius
#' @param sigma (numeric) rNN parameter - standard deviation of Gaussian kernel
#' @param alpha (numeric) rNN parameter - determines radius used
#' @param k_spatial (numeric) rNN parameter - number of neighbors to use
#' @param spatial_mode (character) spatial mode to use (kNN_r, kNN_rn, kNN_rank,
#'   kNN_unif, rNN_gauss)
#' @param assay_name (character) Name for Banksy assay in Seurat object
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
RunBanksy <- function(object, lambda, assay='RNA', slot='data',
                      dimx=NULL, dimy=NULL, features='variable',
                      M=0, k_geom=15, n=2, sigma=1.5,
                      alpha=0.05, k_spatial=10, spatial_mode='kNN_r',
                      assay_name='BANKSY', verbose=TRUE) {
    # Check packages
    SeuratWrappers:::CheckPackage(package = 'data.table', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Matrix', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Banksy', repository = 'github')

    # Check lambda param
    if (lambda < 0 || lambda > 1) stop('Lambda must be between 0 and 1')

    # Get data
    data_own <- get_data(object, assay, slot, features, verbose)

    # Get locs
    locs <- get_locs(object, dimx, dimy, data_own, verbose)

    # Compute neighbor matrix
    if (verbose) message('Computing neighbors')
    knn_list <- lapply(k_geom, function(kg) {
      Banksy:::computeNeighbors(locs,
                                spatial_mode = spatial_mode, k_geom = kg, n = n,
                                sigma=sigma, alpha=alpha, k_spatial=k_spatial,
                                verbose=verbose)
    })


    # Create Banksy matrix
    M <- seq(0, M)
    # Compute harmonics
    center <- rep(TRUE, length(M))
    # Only center higher harmonics
    center[1] <- FALSE
    har <- Map(function(knn_df, M, center) {
      Banksy:::computeHarmonics(data_own, knn_df, M, center, verbose)
    }, knn_list, M, center)

    # Scale by lambdas
    lambdas <- Banksy:::getLambdas(lambda, n_harmonics = length(har))

    # Merge with own expression
    if (verbose) message('Creating Banksy matrix')
    data_banksy <- c(list(data_own), har)
    data_scaled <- lapply(data_banksy, Seurat:::FastRowScale)

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
    object <- SetAssayData(object, slot = 'scale.data', new.data = data_scaled,
                           assay = assay_name)

    # Log commands
    object <- Seurat::LogSeuratCommand(object = object)

  return(object)
}

# Get own expression matrix from Seurat object
get_data <- function(object, assay, slot, features, verbose) {
    # Fetch data from Seurat
    if (verbose) message('Fetching data from slot ', slot,' from assay ', assay)
    data_own <- Seurat::GetAssayData(object = object, assay = assay, slot = slot)
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
get_locs <- function(object, dimx, dimy, data_own, verbose) {
    if (!is.null(dimx) & !is.null(dimy)) {
        # Convert locations data to data table
        locations <- data.frame(
            sdimx = unlist(object[[dimx]]),
            sdimy = unlist(object[[dimy]])
        )
        rownames(locations) <- colnames(object)
        locs <- data.table::data.table(locations, keep.rownames = TRUE)
        data.table::setnames(locs, 'rn', 'cell_ID')

        # Check locations
        obj_samples <- colnames(data_own)
        locs_samples <- locs[['cell_ID']]
        if (any(is.na(match(obj_samples, locs_samples)))) {
            na_id <- which(is.na(match(obj_samples, locs_samples)))
            warning('No centroids found for samples: ', paste(obj_samples[na_id], collapse = ', '), '. Dropping samples.')
            data_own <- data_own[, -na_id, drop = FALSE]
        }
        locs <- locs[match(obj_samples, locs_samples),,drop=FALSE]
    } else {
        locations <- Seurat::GetTissueCoordinates(object)[,1:2]
        locs <- data.table::data.table(locations, keep.rownames = TRUE)
    }

    colnames(locs) <- c('cell_ID', 'sdimx', 'sdimy')

    return(locs)
}

