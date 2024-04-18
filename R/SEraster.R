#' @include internal.R
#'
NULL

#' rasterizeGeneExpression
#' 
#' @description Function to rasterize feature x observation matrix in spatially-resolved 
#' omics data represented as SpatialExperiment class.
#'  
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} 
#' object or a \code{list} of \code{SpatialExperiment} objects.

#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom Matrix colSums
#' 
#' @export
#' 
rasterizeGeneExpression <- function(object, assay_name = NULL, resolution = 100, square = TRUE, fun = "mean", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
    ## create bbox
    pos <- GetTissueCoordinates(object)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## rasterize
    if (is.null(assay_name)) {
      out <- SEraster::rasterizeMatrix(object[[assay_name]], pos, bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
    } else {
      stopifnot(is.character(assay_name))
      stopifnot("assay_name does not exist in the input Seurat object"= assay_name %in% Assays(object)
      out <- SEraster::rasterizeMatrix(object[[assay_name]], pos, bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
    }
    data_rast <- out$data_rast
    pos_rast <- out$pos_rast
    meta_rast <- out$meta_rast[,c("num_cell","type","resolution")]
    
    ## construct a new Seurat object as output
    output <- CreateSeuratObject(
      counts = data_rast,
      assay = paste0(assay_name,".rast"),
      meta.data = meta_rast
    )

    coords <- CreateFOV(assay = paste0("rast.",assay_name), coords = as.data.frame(pos_rast), type = c("centroids"), molecules = NULL)
    output[['rasterized']] <- coords
    
    ## return a Seurat object
    return(output)
  }


#' rasterizeCellType
#' 
#' @description Function to rasterize cell type labels in spatially-resolved 
#' omics data represented as SpatialExperiment class.
#'
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} 
#' object or a \code{list} of \code{SpatialExperiment} objects.
#' 
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom Matrix sparse.model.matrix t colSums
#' 
#' @export
#' 
rasterizeCellType <- function(object, col_name, resolution = 100, square = TRUE, fun = "sum", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
    ## create bbox
    pos <- GetTissueCoordinates(object)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## extract cell type labels 
    stopifnot(is.character(col_name))
    stopifnot("col_name does not exist in the input Seurat object"= col_name %in% colnames(object[[]]))
    # to try: object$fake_ct <- rep(c('A', 'B', 'C', 'D'), length.out = ncol(object))
    cellTypes <- as.factor(object@meta.data[,col_name])
    
    ## one-hot encode cell type labels as sparse matrix
    mat_ct <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
    rownames(mat_ct) <- levels(cellTypes)
    colnames(mat_ct) <- rownames(pos)
    
    ## rasterize
    out <- rasterizeMatrix(mat_ct, pos, bbox, resolution = resolution, square = square, fun = fun, n_threads = 1, BPPARAM = BPPARAM, verbose = verbose)
    data_rast <- out$data_rast
    pos_rast <- out$pos_rast
    meta_rast <- out$meta_rast[,c("num_cell","type","resolution")]
    
    ## construct a new Seurat object as output
    output <- CreateSeuratObject(
      counts = data_rast,
      assay = paste0(assay_name,".rast"),
      meta.data = meta_rast
    )

    coords <- CreateFOV(assay = paste0("rast.",assay_name), coords = as.data.frame(pos_rast), type = c("centroids"), molecules = NULL)
    output[['rasterized']] <- coords
    
    ## return a new Seurat object
    return(output)
  }


#' permutateByRotation
#' 
#' @description Function to permutate a given input SpatialExperiment object(s) 
#' by rotating the x,y coordinates around the midrange point.
#' 
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} 
#' object or a \code{list} of \code{SpatialExperiment} objects.
#'
#' @description When the input is a \code{list} of \code{SpatialExperiment} objects, 
#' all \code{SpatialExperiment} objects will be rotated around a common midrange 
#' point computed based on combined x,y coordinates.
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assays colData
#' @importFrom rearrr rotate_2d midrange
#' 
#' @export
#' 
permutateByRotation <- function(object, n_perm = 1, verbose = FALSE) {
    ## compute rotation degrees based on the required number of permutation
    angles <- seq(0, 360, by = 360/n_perm)[1:n_perm]

    if (verbose) {
        message(paste0("Number of permutations = ", n_perm))
        message(paste0("Angles used for rotations are ", paste(angles, collapse = ", "), " degrees"))
    }
  
    stopifnot("input must be a SpatialExperiment object or a list of SpatialExperiment objects"=class(input) == "SpatialExperiment")

    ## get original x,y coordinates
    assay_name <- DefaultAssay(object)
    pos_orig <- GetTissueCoordinates(object)
    colnames(pos_orig) <- c("x","y")
    stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively. Please change column names accordingly."=colnames(pos_orig)[1:2] == c("x", "y"))

    for(angle in angles) {
        ## rotate around the midrange point
        pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin_fn = rearrr::midrange, overwrite = TRUE)
        
        pos_rotated <- as.matrix(pos_rotated[,c("x_rotated", "y_rotated")])
        colnames(pos_rotated) <- c("x", "y")
        rownames(pos_rotated) <- rownames(pos_orig)

        coords <- CreateFOV(assay = assay_name, coords = as.data.frame(pos_rotated), type = c("centroids"), molecules = NULL)
        object[[paste0("rotated_", angle)]] <- coords
    }

    ## return the Seurat object with new fields of view
    return(output)
}
