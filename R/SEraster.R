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
rasterizeGeneExpression <- function(input, assay_name = NULL, resolution = 100, square = TRUE, fun = "mean", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
    if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- GetTissueCoordinates(input[[i]])
      if (!is.null(names(input))) {
        dataset <- names(input)[[i]]
      } else {
        dataset <- i
      }
      return(data.frame(dataset = dataset, xmin = min(pos[,1]), xmax = max(pos[,1]), ymin = min(pos[,2]), ymax = max(pos[,2])))
    }))

    bbox_common <- sf::st_bbox(c(
      xmin = floor(min(bbox_mat$xmin)-resolution/2), 
      xmax = ceiling(max(bbox_mat$xmax)+resolution/2), 
      ymin = floor(min(bbox_mat$ymin)-resolution/2), 
      ymax = ceiling(max(bbox_mat$ymax)+resolution/2)
    ))
    
    ## rasterize iteratively
    output_list <- lapply(seq_along(input), function(i) {
      ## get Seurat object of the given index
      spe <- input[[i]]
      coords <- GetTissueCoordinates(spe)

      if (is.null(assay_name)) {
        assay_name <- DefaultAssay(spe)
        out <- SEraster::rasterizeMatrix(spe[[assay_name]], coords, bbox = bbox_common, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
      } else {
        stopifnot(is.character(assay_name))
        stopifnot("assay_name does not exist in the input SpatialExperiment object"= assay_name %in% Assays(spe))
        out <- SEraster::rasterizeMatrix(spe[[assay_name]], coords, bbox = bbox_common, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
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

      coords <- CreateFOV(coords = as.data.frame(pos_rast), assay = paste0(assay_name, ".rast"), type = c("centroids"), molecules = NULL)
      output[['rasterized']] <- coords
      return(output)
    })
    
    if (!is.null(names(input))) {
      names(output_list) <- names(input)
    }
    
    ## return a list of Seurat objects
    return(output_list)
  } else {
    ## create bbox
    pos <- GetTissueCoordinates(input)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## rasterize
    if (is.null(assay_name)) {
      assay_name <- DefaultAssay(input)
      out <- SEraster::rasterizeMatrix(input[[assay_name]], pos, bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
    } else {
      stopifnot(is.character(assay_name))
      stopifnot("assay_name does not exist in the input Seurat object"= assay_name %in% Assays(object))
      out <- SEraster::rasterizeMatrix(input[[assay_name]], pos, bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
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

    coords <- CreateFOV(coords = as.data.frame(pos_rast), assay = paste0(assay_name, ".rast"), type = c("centroids"), molecules = NULL)
    output[['rasterized']] <- coords
    
    ## return a Seurat object
    return(output)
  }
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
rasterizeCellType <- function(input, col_name, resolution = 100, square = TRUE, fun = "sum", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
    if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- GetTissueCoordinates(input[[i]])
      if (!is.null(names(input))) {
        dataset <- names(input)[[i]]
      } else {
        dataset <- i
      }
      return(data.frame(dataset = dataset, xmin = min(pos[,1]), xmax = max(pos[,1]), ymin = min(pos[,2]), ymax = max(pos[,2])))
    }))
    bbox_common <- sf::st_bbox(c(
      xmin = floor(min(bbox_mat$xmin)-resolution/2), 
      xmax = ceiling(max(bbox_mat$xmax)+resolution/2), 
      ymin = floor(min(bbox_mat$ymin)-resolution/2), 
      ymax = ceiling(max(bbox_mat$ymax)+resolution/2)
    ))
    
    ## rasterize iteratively
    output_list <- lapply(seq_along(input), function(i) {
      ## get Seurat object of the given index
      spe <- input[[i]]
      pos <- GetTissueCoordinates(input[[i]])
      
      ## extract cell type labels from SpatialExperiment
      stopifnot(is.character(col_name))
      stopifnot("col_name does not exist in the input SpatialExperiment object"= col_name %in% colnames(spe[[]]))
      cellTypes <- as.factor(spe@meta.data[,col_name])
      
      ## one-hot encode cell type labels as sparse matrix
      mat_ct <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
      rownames(mat_ct) <- levels(cellTypes)
      colnames(mat_ct) <- rownames(pos)
      
      ## rasterize
      out <- SEraster::rasterizeMatrix(mat_ct, pos, bbox_common, resolution = resolution, square = square, fun = fun, n_threads = 1, BPPARAM = BPPARAM, verbose = verbose)
      data_rast <- out$data_rast
      pos_rast <- out$pos_rast
      meta_rast <- out$meta_rast[,c("num_cell","type","resolution")]
      
      ## construct a new Seurat object[[ as output
      assay_name <- DefaultAssay(input[[i]])
      output <- CreateSeuratObject(
        counts = data_rast,
        assay = paste0(assay_name,".rast"),
        meta.data = meta_rast
      )

      coords <- CreateFOV(assay = paste0(assay_name, ".rast"), coords = as.data.frame(pos_rast), type = c("centroids"), molecules = NULL)
      output[['rasterized']] <- coords
      
      ## return a new Seurat object
      return(output)
    })
    
    if (!is.null(names(input))) {
      names(output_list) <- names(input)
    }
    
    ## return a list of SpatialExperiment object
    return(output_list)
    
  } else {
    ## create bbox
    pos <- GetTissueCoordinates(input)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## extract cell type labels 
    stopifnot(is.character(col_name))
    stopifnot("col_name does not exist in the input Seurat object"= col_name %in% colnames(input[[]]))
    cellTypes <- as.factor(input@meta.data[,col_name])
    
    ## one-hot encode cell type labels as sparse matrix
    mat_ct <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
    rownames(mat_ct) <- levels(cellTypes)
    colnames(mat_ct) <- rownames(pos)
    
    ## rasterize
    out <- SEraster::rasterizeMatrix(mat_ct, pos, bbox, resolution = resolution, square = square, fun = fun, n_threads = 1, BPPARAM = BPPARAM, verbose = verbose)
    data_rast <- out$data_rast
    pos_rast <- out$pos_rast
    meta_rast <- out$meta_rast[,c("num_cell","type","resolution")]
    
    ## construct a new Seurat object as output
    assay_name <- DefaultAssay(input)
    output <- CreateSeuratObject(
      counts = data_rast,
      assay = paste0(assay_name,".rast"),
      meta.data = meta_rast
    )

    coords <- CreateFOV(coords = as.data.frame(pos_rast), assay = paste0(assay_name, ".rast"), type = c("centroids"), molecules = NULL)
    output[['rasterized']] <- coords
    
    ## return a new Seurat object
    return(output)
  }
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
permutateByRotation <- function(input, n_perm = 1, verbose = FALSE) {
  ## compute rotation degrees based on the required number of permutation
  angles <- seq(0, 360, by = 360/n_perm)[1:n_perm]

  if (verbose) {
      message(paste0("Number of permutations = ", n_perm))
      message(paste0("Angles used for rotations are ", paste(angles, collapse = ", "), " degrees"))
  }

  if (is.list(input)) {
    ## combine all x,y coordinates
    pos_comb <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- GetTissueCoordinates(input[[i]])
      if (!is.null(names(input))) {
        dataset <- names(input)[[i]]
      } else {
        dataset <- i
      }
      return(data.frame(dataset = dataset, x = pos[,1], y = pos[,2]))
    }))
    ## find the midrange point across combined x,y coordinates
    midrange_pt <- rearrr::midrange(pos_comb, cols = c("x", "y"))
    
    if (verbose) {
      message(paste0("Rotating all datasets around (x, y) = (", midrange_pt$x, ", ", midrange_pt$y, ")."))
    }
    
    output <- unlist(lapply(input, function(spe) {      
      ## get original x,y coordinates
      pos_orig <- data.frame(GetTissueCoordinates(spe))
      stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively. Please change column names accordingly."=colnames(pos_orig)[1:2] == c("x", "y"))
      
      output2 <- lapply(angles, function(angle) {
        ## rotate around the midrange point
        pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin = as.numeric(midrange_pt), overwrite = TRUE)
        
        pos_rotated <- as.matrix(pos_rotated[,c("x_rotated", "y_rotated")])
        colnames(pos_rotated) <- c("x", "y")
        rownames(pos_rotated) <- rownames(pos_orig)
        
        ## update SpatialExperiment object
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = SummarizedExperiment::assays(spe),
          spatialCoords = pos_rotated,
          colData = SummarizedExperiment::colData(spe)
        )
        return(spe_rotated)
      })
      return(output2)
    }))
    
    ## assign list names
    if (!is.null(names(input))) {
      names(output) <- paste0(rep(names(input), each = length(angles)), "_rotated_", angles)
    }

    ## return a list of SpatialExperiment objects
    return(output)
  
  } else {

    ## get original x,y coordinates
    assay_name <- DefaultAssay(input)
    pos_orig <- GetTissueCoordinates(input)
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
}


## test script
library(Seurat)
library(SeuratData)
brain <- LoadData("stxBrain", type = "anterior1")
brain2 <- LoadData("stxBrain", type = "posterior1")

fake_col1 <- rep(c("A", "B", "C"), length.out = ncol(brain))
fake_col2 <- rep(c("D", "E", "F"), length.out = ncol(brain2))

brain$fake_celltypes <- fake_col1
brain2$fake_celltypes <- fake_col2

brains <- list(brain, brain2)

brains_rast <- rasterizeGeneExpression(brains, assay_name = NULL, resolution = 100, square = TRUE, fun = "mean", n_threads = 1, BPPARAM = NULL, verbose = FALSE) 
brains_cells_rast <- rasterizeCellType(brains, col_name = "fake_celltypes", resolution = 100, square = TRUE, fun = "sum", n_threads = 1, BPPARAM = NULL, verbose = FALSE)