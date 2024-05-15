#' @include internal.R
#'
NULL
#' createRasterizedObject
#' @keyword internal
#' 
createRasterizedObject <- function(input, out, name) {
  data_rast <- out$data_rast
  pos_rast_temp <- as.data.frame(out$pos_rast)
  
  if(name %in% Assays(input)){
    image <- Images(input, assay = name)[1]
  } else{
    image <- Images(input, assay = DefaultAssay(input))[1]
  }
  
  meta_rast <- out$meta_rast[,c("num_cell","type","resolution")]
  resolution <- meta_rast[["resolution"]][1]
  for (i in seq_along(out$meta_rast$cellID_list)) {
    meta_rast$cellID_list[i] <- paste(unlist(out$meta_rast$cellID_list[[i]]), collapse = ", ")
  }
  
  output <- CreateSeuratObject(
    counts = data_rast,
    assay = name,
    meta.data = meta_rast
  )
  
  class <- class(input[[image]])
  if(class == 'VisiumV1') {
    scale.use <- ScaleFactors(input[[image]])[['lowres']]
    pos_rast <- as.data.frame(cbind(tissue = rep(1, nrow(pos_rast_temp)), 
                                     row = pos_rast_temp$x, 
                                     col = pos_rast_temp$y,
                                     imagerow = ceiling(pos_rast_temp$x / scale.use),
                                     imagecol = ceiling(pos_rast_temp$y / scale.use)))
    rownames(pos_rast) <- rownames(pos_rast_temp)
    input[[image]]@scale.factors$spot <- sqrt(nrow(input[[]])/nrow(meta_rast))*Radius(input[[image]])
    
    new.fov <- new(
      Class = "VisiumV1",
      coordinates = pos_rast,
      assay = name, 
      key = "rasterized_",
      image = input[[image]]@image,
      scale.factors = input[[image]]@scale.factors,
      spot.radius = input[[image]]@scale.factors$spot
    )
  } else if(class == 'VisiumV2') {
    input[[image]]@scale.factors$spot <- sqrt(nrow(input[[]])/nrow(meta_rast))*slot(input[[image]]$centroids, name='radius')
    fov <- CreateFOV(
      pos_rast_temp[, c("x", "y")],
      type = "centroids",
      radius = input[[image]]@scale.factors[["spot"]],
      assay = name,
      key = Key("rasterized_", quiet = TRUE)
    )
    new.fov <- new(
      Class = "VisiumV2",
      boundaries = fov@boundaries,
      assay = fov@assay,
      key = fov@key,
      image = input[[image]]@image,
      scale.factors = input[[image]]@scale.factors
    )
  } else if(class == 'FOV') {
    radius <- sqrt(nrow(input[[]])/nrow(meta_rast))*slot(input[[image]]$centroids, name='radius')
    pos_rast_temp$cell <- rownames(pos_rast_temp)
    segmentations.data <- list(
      "centroids" = CreateCentroids(coords = pos_rast_temp[, c("x", "y")], nsides = 4, radius = radius),
      "segmentation" = CreateSegmentation(coords = pos_rast_temp)
    )
    new.fov <- CreateFOV(
      coords = segmentations.data,
      type = c("segmentation", "centroids"),
      radius = radius,
      molecules = input[[image]]$molecules,
      assay = name,
      key = Key("rasterized_", quiet = TRUE)
    )
  }
  output@images[[paste0(name,".rasterized")]]<- new.fov
  return(output)
}

#' rasterizeGeneExpression
#' @export
#' 
rasterizeGeneExpression <- function(input, assay_name = NULL, image = NULL, slot = "counts", resolution = 100, square = TRUE, fun = "mean", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
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
      image <- ifelse(is.null(image), Images(spe, assay=assay_name)[1], image)
      coords <- GetTissueCoordinates(spe, image = image)
      if("cell" %in% colnames(coords)){
        rownames(coords) <- coords$cell
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
    pos <- GetTissueCoordinates(input, image = image)
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
rasterizeCellType <- function(input, col_name, assay_name = NULL, image = NULL, resolution = 100, square = TRUE, fun = "sum", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
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
      image <- ifelse(is.null(image), Images(spe, assay=assay_name)[1], image)
      coords <- GetTissueCoordinates(spe, image = image)
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
    pos <- GetTissueCoordinates(input, image = image)
    if("cell" %in% colnames(pos)){
        rownames(coords) <- pos$cell
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
  ## compute rotation degrees based on the required number of permutation
  angles <- seq(0, 360, by = 360/n_perm)[1:n_perm]
  
  if (verbose) {
    message(paste0("Number of permutations = ", n_perm))
    message(paste0("Angles used for rotations are ", paste(angles, collapse = ", "), " degrees"))
  }
  
  if (is.list(input)) {
    ## combine all x,y coordinate
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
      assay_name <- DefaultAssay(spe)
      pos_orig <- data.frame(GetTissueCoordinates(spe))
      colnames(pos_orig) <- c("x","y")
      stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively. Please change column names accordingly."=colnames(pos_orig)[1:2] == c("x", "y"))
      
      output2 <- lapply(angles, function(angle) {
        ## rotate around the midrange point
        pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin = as.numeric(midrange_pt), overwrite = TRUE)
        
        pos_rotated <- as.data.frame(pos_rotated[,c("x_rotated", "y_rotated")])
        colnames(pos_rotated) <- c("x", "y")
        rownames(pos_rotated) <- rownames(pos_orig)
        
        image <- Images(spe, assay = assay_name)[[1]]
        class <- class(spe[[image]])

        if(class == 'VisiumV1') {
          scale.use <- ScaleFactors(spe[[image]])[['lowres']]
          pos_new <- slot(spe[[image]], name="coordinates")
          pos_new$row <- pos_rotated$x
          pos_new$col <- pos_rotated$y 
          pos_new$imagerow <- ceiling(pos_rotated$x / scale.use)
          pos_new$imagecol <- ceiling(pos_rotated$y / scale.use)
          rownames(pos_new) <- rownames(pos_rotated)
          visium.fov <- new(
            Class = "VisiumV1",
            coordinates = pos_new,
            assay = assay_name, 
            key = paste0("rotate",angle,"_"),
            image = spe[[image]]@image,
            scale.factors = spe[[image]]@scale.factors,
            spot.radius = spe[[image]]@spot.radius
          )
        } else if(class == 'VisiumV2') {
          pos_new <-  pos_rotated[, c("x", "y")]
          fov <- CreateFOV(
            pos_new[, c("x", "y")],
            type = "centroids",
            radius = spe[[image]]@scale.factors[["spot"]],
            assay = assay_name,
            theta = angle,
            key =paste0("rotate",angle,"_")
          )
          visium.fov <- new(
            Class = "VisiumV2",
            boundaries = fov@boundaries,
            molecules = fov@molecules,
            assay = fov@assay,
            key = fov@key,
            image = spe[[image]]@image,
            scale.factors = spe[[image]]@scale.factors
          )
        }
        spe@images[[paste0("rotate",angle)]] <- visium.fov
        return(spe)
      })
      return(output2)
    }))
    
    if (!is.null(names(input))) {
      names(output) <- paste0(rep(names(input), each = length(angles)), "_rotated_", angles)
    }
    
    return(output)
    
  } else {
    ## get original x,y coordinates
    assay_name <- DefaultAssay(input)
    pos_orig <- GetTissueCoordinates(input)
    colnames(pos_orig) <- c("x","y")
    stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively. Please change column names accordingly."=colnames(pos_orig)[1:2] == c("x", "y"))
    
    for(angle in angles) {
      ## rotate around the midrange point
      midrange_pt <- rearrr::midrange(pos_orig, cols = c("x", "y"))
      pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin = as.numeric(midrange_pt), overwrite = TRUE)
      
      pos_rotated <- as.data.frame(pos_rotated[,c("x_rotated", "y_rotated")])
      colnames(pos_rotated) <- c("x", "y")
      rownames(pos_rotated) <- rownames(pos_orig)
      
      image <- Images(input, assay = assay_name)[[1]]
      class <- class(input[[image]])
      if(class == 'VisiumV1') {
        scale.use <- ScaleFactors(input[[image]])[['lowres']]
        pos_new <- slot(input[[image]], name="coordinates")
        pos_new$row <- pos_rotated$x
        pos_new$col <- pos_rotated$y 
        pos_new$imagerow <- ceiling(pos_rotated$x / scale.use)
        pos_new$imagecol <- ceiling(pos_rotated$y / scale.use)
        rownames(pos_new) <- rownames(pos_rotated)
        
        new.fov <- new(
          Class = "VisiumV1",
          coordinates = pos_new,
          assay = assay_name, 
          key = paste0("rotate",angle,"_"),
          image = input[[image]]@image,
          scale.factors = input[[image]]@scale.factors,
          spot.radius = input[[image]]@spot.radius
        )
      } else if(class == 'VisiumV2') {
        pos_new <-  pos_rotated[, c("x", "y")]
        fov <- CreateFOV(
          pos_new[, c("x", "y")],
          type = "centroids",
          radius = input[[image]]@scale.factors[["spot"]],
          assay = assay_name,
          theta = angle,
          key = paste0("rotate",angle,"_")
        )
        new.fov <- new(
          Class = "VisiumV2",
          boundaries = fov@boundaries,
          molecules = fov@molecules,
          assay = fov@assay,
          key = fov@key,
          image = rotate(input[[image]]@image, angle),
          scale.factors = input[[image]]@scale.factors
        )
      } else if(class == 'FOV') {
        pos_new <-  pos_rotated[, c("x", "y")]
        pos_new$cell <- rownames(pos_rotated)
        radius <- slot(input[[image]]$centroids, name='radius')
        segmentations.data <- list(
          "centroids" = CreateCentroids(coords = pos_new[, c("x", "y")], radius = radius),
          "segmentation" = CreateSegmentation(coords = pos_new)
        )
        new.fov <- CreateFOV(
          coords = segmentations.data,
          type = c("segmentation", "centroids"),
          radius = radius,
          molecules = input[[image]]$molecules,
          assay = assay_name
        )
      }
      input@images[[paste0("rotate",angle)]] <- new.fov
    }
    ## return the original Seurat object with new fields of view
    return(input)
  }
}

#' rotate
#' @keyword internal
#' 
rotate <- function(arr, angle) {
  angle_rad <- angle * pi / 180
  rows <- dim(arr)[1]
  cols <- dim(arr)[2]

  midrange_pt <- rearrr::midrange(pos, cols = c("x", "y"))
  max_x <- max(pos['x'])
  max_y <- max(pos['y'])
  center_x <- cols * midrange_pt$x / max_x
  center_y <- rows * midrange_pt$y / max_y
  center_x <- 280

  size <- max(center_x*2,center_y*2)

  rotated_arr <- array(0, dim = c(size,size,3))
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      x <- j - center_x 
      y <- i - center_y

      rotated_x <- x * cos(angle_rad) - y * sin(angle_rad)
      rotated_y <- x * sin(angle_rad) + y * cos(angle_rad)
      
      rotated_x <- rotated_x + center_x
      rotated_y <- rotated_y + center_y 

      if (rotated_x >= 1 && rotated_x <= dim(rotated_arr)[2] && rotated_y >= 1 && rotated_y <= dim(rotated_arr)[1]) {
        rotated_arr[ceiling(rotated_y), ceiling(rotated_x), ] <- arr[i, j, ]
      } else {
        print(rotated_x)
      }
    }
  }
  return(rotated_arr)
}
