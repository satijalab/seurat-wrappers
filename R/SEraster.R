#' @include internal.R
#'
NULL
#' createRasterizedObject
#' @keyword internal
#' 
createRasterizedObject <- function(input, out, name) {
  image <- ifelse(name %in% Assays(input), Images(input, assay = name)[1], Images(input, assay = DefaultAssay(input))[1])
  input_fov <- input[[image]]

  data_rast <- out$data_rast
  meta_rast <- out$meta_rast[,c("num_cell","type","resolution")]
  resolution <- meta_rast[["resolution"]][1]
  for (i in seq_along(out$meta_rast$cellID_list)) {
    meta_rast$cellID_list[i] <- paste(unlist(out$meta_rast$cellID_list[[i]]), collapse = ", ")
  }

  output_image_name <- paste0("ras.", resolution,".", image)
  output_coordinates <- as.data.frame(out$pos_rast)
  
  output <- CreateSeuratObject(
    counts = data_rast,
    assay = name,
    meta.data = meta_rast
  )

  input_molecules <- tryCatch(input_fov[["molecules"]], error = function(e) NULL)

  # `scale_factors` will be set to `NULL` unless there a matching 
  # implementation of the `ScaleFactors` generic available for `input_fov` 
  scale_factors <- tryCatch(
    ScaleFactors(input_fov),
    error = function(e) {
      return (NULL)
    }
  )

  output_radius <- sqrt(nrow(input[[]]) / nrow(meta_rast)) * ifelse(!is.null(Radius(input_fov, scale = NULL)), Radius(input_fov, scale = NULL), Radius(input_fov[['centroids']], scale = NULL))
  
  output_centroids <- CreateCentroids(
    coords = output_coordinates,
    radius = output_radius
  )

  output_fov <- CreateFOV(
    coords = output_centroids,
    type = 'centroids',
    molecules = input_molecules,
    assay = name,
    key = Key(output_image_name, quiet = TRUE)
  )
  
  if (!is.null(scale_factors)) {
    scale_factors$spot <- output_radius
    output_fov <- new(
      Class = "VisiumV2",
      boundaries = output_fov@boundaries,
      molecules = output_fov@molecules,
      assay = output_fov@assay,
      key = output_fov@key,
      image = input_fov@image,
      scale.factors = scale_factors
    )
  } 

  output[[output_image_name]] <- output_fov
  return(output)
}

#' rasterizeGeneExpression
#' @export
#' 
rasterizeGeneExpression <- function(
  input, 
  assay_name = NULL, 
  image = NULL, 
  slot = "counts", 
  resolution = 100, 
  square = TRUE, 
  fun = "mean", 
  n_threads = 1, 
  BPPARAM = NULL, 
  verbose = FALSE
) {
  if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- GetTissueCoordinates(input[[i]],scale = NULL)
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
      image <- ifelse(is.null(image), Images(spe, assay=assay_name)[1], image)
      coords <- GetTissueCoordinates(spe, image = image, scale = NULL)
      if("cell" %in% colnames(coords)){
        rownames(coords) <- colnames(spe)
      }
      
      if (is.null(assay_name)) {
        assay_name <- DefaultAssay(spe)
        counts <- LayerData(spe[[assay_name]], layer = slot)
        out <- SEraster::rasterizeMatrix(counts, coords, bbox = bbox_common, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
      } else {
        stopifnot(is.character(assay_name))
        stopifnot("assay_name does not exist in the input Seurat object"= assay_name %in% Assays(spe))
        counts <- LayerData(spe[[assay_name]], layer = slot)
        out <- SEraster::rasterizeMatrix(counts, coords, bbox = bbox_common, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
      }
      
      output <- createRasterizedObject(spe, out, assay_name)
      return(output)
    })
    
    if (!is.null(names(input))) {
      names(output_list) <- names(input)
    }
    ## return a list of Seurat objects
    return(output_list)
  } else {
    ## create bbox
    image <- ifelse(is.null(image), Images(input, assay=assay_name)[1], image)
    pos <- GetTissueCoordinates(input, image = image, scale = NULL)
    if("cell" %in% colnames(pos)){
      rownames(pos) <- pos$cell
    }

    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## rasterize
    if (is.null(assay_name)) {
      assay_name <- DefaultAssay(input)
      counts <- LayerData(input[[assay_name]], layer = slot)
      out <- SEraster::rasterizeMatrix(counts, pos, bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
    } else {
      stopifnot(is.character(assay_name))
      stopifnot("assay_name does not exist in the input Seurat object"= assay_name %in% Assays(input))
      counts <- LayerData(input[[assay_name]], layer = slot)
      out <- SEraster::rasterizeMatrix(counts, pos, bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
    }
    
    output <- createRasterizedObject(input=input, out=out, name=assay_name)
    return(output)
  }
}

#' rasterizeCellType
#' 
#' @description Function to rasterize cell type labels in spatially-resolved 
#' omics data represented as Seurat object class.
#'
#' @description This function assumes that the input is provided as a \code{Seurat object} 
#' object or a \code{list} of \code{Seurat} objects.
#' 
#' @importFrom Matrix sparse.model.matrix t colSums
#' 
#' @export
#' 
rasterizeCellType <- function(
  input, 
  col_name, 
  assay_name = NULL, 
  image = NULL, 
  resolution = 100, 
  square = TRUE, 
  fun = "sum", 
  n_threads = 1, 
  BPPARAM = NULL, 
  verbose = FALSE
) {
  if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- GetTissueCoordinates(input[[i]], scale = NULL)
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
      image <- ifelse(is.null(image), Images(spe, assay=assay_name)[1], image)
      coords <- GetTissueCoordinates(spe, image = image, scale = NULL)
      if("cell" %in% colnames(coords)){
        rownames(coords) <- coords$cell
      }
      
      ## extract cell type labels from Seurat object 
      stopifnot(is.character(col_name))
      stopifnot("col_name does not exist in the input Seurat object"= col_name %in% colnames(spe[[]]))
      cellTypes <- as.factor(spe@meta.data[,col_name])
      
      ## one-hot encode cell type labels as sparse matrix
      mat_ct <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
      rownames(mat_ct) <- levels(cellTypes)
      colnames(mat_ct) <- rownames(coords)
      
      ## rasterize
      out <- SEraster::rasterizeMatrix(mat_ct, coords, bbox_common, resolution = resolution, square = square, fun = fun, n_threads = 1, BPPARAM = BPPARAM, verbose = verbose)
      output <- createRasterizedObject(spe, out, col_name)
      return(output)
    })
    
    if (!is.null(names(input))) {
      names(output_list) <- names(input)
    }
    
    ## return a list of Seurat objects
    return(output_list)
    
  } else {
    ## create bbox
    image <- ifelse(is.null(image), Images(input, assay=assay_name)[1], image)
    pos <- GetTissueCoordinates(input, image = image, scale = NULL)
    if("cell" %in% colnames(pos)){
        rownames(pos) <- pos$cell
    }
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
    output <- createRasterizedObject(input=input, out=out, name=col_name)
    return(output)
  }
}

#' permutateByRotation
#' 
#' @description Function to permutate a given input Seurat object(s) 
#' by rotating the x,y coordinates around the midrange point.
#' 
#' @description This function assumes that the input is provided as a \code{Seurat} 
#' object or a \code{list} of \code{Seurat} objects.
#'
#' @description When the input is a \code{list} of \code{Seurat} objects, 
#' all \code{Seurat} objects will be rotated around a common midrange 
#' point computed based on combined x,y coordinates.
#' 
#' @importFrom rearrr rotate_2d midrange
#' 
#' @export
#' 
permutateByRotation <- function(input, n_perm = 1, verbose = FALSE) {
  angles <- seq(0, 360, by = 360/n_perm)[1:n_perm]
  
  if (verbose) {
    message(paste0("Number of permutations = ", n_perm))
    message(paste0("Angles used for rotations are ", paste(angles, collapse = ", "), " degrees"))
  }
  
  if (is.list(input)) {
    pos_comb <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- GetTissueCoordinates(input[[i]], scale = NULL)
      dataset <- ifelse(!is.null(names(input)), names(input)[[i]], i)
      return(data.frame(dataset = dataset, x = pos[, 1], y = pos[, 2]))
    }))
    midrange_pt <- rearrr::midrange(pos_comb, cols = c("x", "y"))
    
    if (verbose) {
      message(paste0("Rotating all datasets around (x, y) = (", midrange_pt$x, ", ", midrange_pt$y, ")."))
    }
    
    output <- unlist(lapply(input, function(spe) {
      assay_name <- DefaultAssay(spe)
      pos_orig <- data.frame(GetTissueCoordinates(spe, scale = NULL))
      colnames(pos_orig) <- c("x", "y")
      stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively." = colnames(pos_orig)[1:2] == c("x", "y"))
      
      lapply(angles, function(angle) {
        pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin = as.numeric(midrange_pt), overwrite = TRUE)
        pos_rotated <- as.data.frame(pos_rotated[, c("x_rotated", "y_rotated")])
        colnames(pos_rotated) <- c("x", "y")
        rownames(pos_rotated) <- rownames(pos_orig)
        
        image_name <- Images(spe, assay = assay_name)[[1]]
        class <- class(spe[[image_name]])
        
        if (class == 'VisiumV1') {
          input <- updateVisiumV1(spe, pos_rotated, assay_name, angle, image_name)
        } else if (class == 'VisiumV2') {
          input <- updateVisiumV2(spe, pos_rotated, assay_name, angle, image_name)
        } else if (class == 'FOV') {
          input <- updateFOV(spe, pos_rotated, assay_name, angle, image_name)
        }
      })
    }))
    
    return(output)
    
  } else {
    assay_name <- DefaultAssay(input)
    pos_orig <- GetTissueCoordinates(input, scale = NULL)
    colnames(pos_orig) <- c("x", "y")
    stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively." = colnames(pos_orig)[1:2] == c("x", "y"))
    
    output_list <- list()
    for (angle in angles) {
      image_name <- Images(input, assay = assay_name)[[1]]
      image <- input[[image_name]]
      midrange_pt <- rearrr::midrange(pos_orig, cols = c("x", "y"))
      pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin = as.numeric(midrange_pt), overwrite = TRUE)
      colnames(pos_rotated) <- c("x", "y")
      rownames(pos_rotated) <- rownames(pos_orig)
      
      class <- class(input[[image_name]])
      
      if (class == 'VisiumV1') {
        output <- updateVisiumV1(input, pos_rotated, assay_name, angle, image_name)
      } else if (class == 'VisiumV2') {
        output <- updateVisiumV2(input, pos_rotated, assay_name, angle, image_name)
      } else if (class == 'FOV') {
        output <- updateFOV(input, pos_rotated, assay_name, angle, image_name)
      }
      output_list <- append(output_list, output)
    }
    return(output_list)
  }
}

#' @keyword internal
#' 
updateVisiumV1 <- function(input, pos_rotated, assay_name, angle, image_name) {
  scale.use <- ScaleFactors(input[[image_name]])[['lowres']]
  pos_old <- slot(input[[image_name]], name = "coordinates")
  pos_new <- as.data.frame(cbind(tissue = pos_old$tissue, 
                                     row = pos_rotated$x, 
                                     col = pos_rotated$y,
                                     imagerow = ceiling(pos_rotated$x / scale.use),
                                     imagecol = ceiling(pos_rotated$y / scale.use)))
  rownames(pos_new) <- rownames(pos_rotated)
  input_fov <- input[[image_name]]

  visium.fov <- new(
    Class = "VisiumV1",
    coordinates = pos_new,
    assay = assay_name,
    key = paste0("rotated", angle, "_"),
    image = rotate(input_fov@image, angle, pos_new),
    scale.factors = input_fov@scale.factors,
    spot.radius = input_fov@spot.radius
  )
  input[[image_name]] <- NULL
  input@images[[paste0("rotated", angle)]] <- visium.fov
  return(input)
}

#' @keyword internal
#' 
updateVisiumV2 <- function(input, pos_rotated, assay_name, angle, image_name) {
  pos_new <- pos_rotated[, c("x", "y")]
  input_fov <- input[[image_name]]

  fov <- CreateFOV(
    pos_new[, c("x", "y")],
    type = "centroids",
    radius = input_fov@scale.factors[["spot"]],
    assay = assay_name,
    theta = angle,
    key = paste0("rotated", angle, "_")
  )
  visium.fov <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = rotate(input_fov@image, angle, pos_new),
    scale.factors = input_fov@scale.factors
  )
  input[[image_name]] <- NULL
  input@images[[paste0("rotated", angle)]] <- visium.fov
  return(input)
}

#' @keyword internal
#' 
updateFOV <- function(input, pos_rotated, assay_name, angle, image_name) {
  pos_new <- pos_rotated[, c("x", "y")]
  pos_new$cell <- rownames(pos_rotated)
  input_fov <- input[[image_name]]
  radius <- slot(input_fov$centroids, name = 'radius')

  output_boundaries <- list(
    "centroids" = CreateCentroids(coords = pos_new[, c("x", "y")], nsides = input_fov[["centroids"]]@nsides, radius = radius, theta = angle),
    "segmentation" = input_fov[["segmentation"]]
  )

  new.fov <- CreateFOV(
    coords = output_boundaries,
    molecules = input_fov$molecules,
    assay = assay_name,
    key = Key(paste0("rotated", angle), quiet = TRUE)
  )
  input[[image_name]] <- NULL
  input@images[[paste0("rotated", angle)]] <- new.fov
  return(input)
}

#' rotate
#' @keyword internal
#' 
rotate <- function(arr, angle, pos) {
  angle_rad <- angle * pi / 180
  rows <- dim(arr)[1]
  cols <- dim(arr)[2]

  center_x <- cols / 2 
  center_y <- rows / 2

  size <- max(center_x*2,center_y*2)

  rotated_arr <- array(0, dim = c(size,size,3))
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      x <- j - center_x 
      y <- i - center_y

      rotated_x <- x * cos(-angle_rad) - y * sin(-angle_rad)
      rotated_y <- x * sin(-angle_rad) + y * cos(-angle_rad)
      
      rotated_x <- rotated_x + center_x
      rotated_y <- rotated_y + center_y 

      if (rotated_x >= 1 && rotated_x <= dim(rotated_arr)[2] && rotated_y >= 1 && rotated_y <= dim(rotated_arr)[1]) {
        rotated_arr[ceiling(rotated_y), ceiling(rotated_x), ] <- arr[i, j, ]
      } 
    }
  }
  return(rotated_arr)
}
