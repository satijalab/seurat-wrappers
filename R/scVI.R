#' @include internal.R
#'
NULL

#' scVI Integration
#' @param object A \code{StdAssay} or \code{STDAssay} instance containing
#' merged data
#' @param features Features to integrate
#' @param layers Layers to integrate
#' @param conda_env conda environment to run scVI
#' @param new.reduction Name under which to store resulting DimReduc object
#' @param ndims Dimensionality of the latent space
#' @param nlayers Number of hidden layers used for encoder and decoder NNs
#' @param gene_likelihood Distribution to use for modelling expression 
#' data: {"zinb", "nb", "poisson"}
#' @param max_epochs Number of passes through the dataset taken while
#' training the model
#' @param ... Unused - currently just capturing parameters passed in from
#' \code{Seurat::IntegrateLayers} intended for other integration methods
#'
#' @export
#'
#' @note This function requires the
#' \href{https://docs.scvi-tools.org/en/stable/installation.html}{\pkg{scvi-tools}}
#' package to be installed
#'
#' @examples
#' \dontrun{
#' # Preprocessing
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' # After preprocessing, we integrate layers, specifying a conda environment
#' obj <- IntegrateLayers(
#'   object = obj,
#'   method = scVIIntegration,
#'   new.reduction = "integrated.scvi",
#'   conda_env = "../miniconda3/envs/scvi-env",
#'   verbose = FALSE
#' )
#'
#' # Alternatively, we can integrate SCTransformed data
#' obj <- SCTransform(object = obj)
#' obj <- IntegrateLayers(
#'   object = obj, method = scVIIntegration,
#'   orig.reduction = "pca", new.reduction = "integrated.scvi",
#'   assay = "SCT", conda_env = "../miniconda3/envs/scvi-env", verbose = FALSE
#' )
#' }
#'
#' @seealso \href{https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scvi_in_R.html}{scVI}
#'
#' @return A single-element named list \code{DimReduc} elements containing
#' the integrated data
scVIIntegration <- function(
    object,
    features = NULL,
    layers = "counts",
    conda_env = NULL,
    new.reduction = "integrated.dr",
    ndims = 30,
    nlayers = 2,
    gene_likelihood = "nb",
    max_epochs = NULL,
    ...) {
  
  # import python methods from specified conda env
  reticulate::use_condaenv(conda_env, required = TRUE)
  sc <- reticulate::import("scanpy", convert = FALSE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  anndata <- reticulate::import("anndata", convert = FALSE)
  scipy <- reticulate::import("scipy", convert = FALSE)

  # if `max_epochs` is not set
  if (is.null(max_epochs)) {
    # convert `NULL` to python's `None`
    max_epochs <- reticulate::r_to_py(max_epochs)
  } else {
    # otherwise make sure it's an int
    max_epochs <- as.integer(max_epochs)
  }

  # build a meta.data-style data.frame indicating the batch for each cell
  # scVI expects a single counts matrix so we'll join our layers together
  # it also expects the raw counts matrix
  # TODO: avoid hardcoding this - users can rename their layers arbitrarily
  # so there's no gauruntee that the usual naming conventions will be followed
  if (inherits(object, what = "SCTAssay")) {
      batches <- .FindSCTBatches(object)
  } else {
      batches <- .FindBatches(object, layers = layers)
      object <- JoinLayers(object = object, layers = "counts")
  }
  
  # setup an `AnnData` python instance
  adata <- sc$AnnData(
    X = scipy$sparse$csr_matrix(
      # TODO: avoid hardcoding per comment above
      Matrix::t(LayerData(object, layer = "counts")[features, ])
    ),
    obs = batches,
    var = object[[]][features, ]
  )
  scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")

  # initialize and train the model
  model <- scvi$model$SCVI(
    adata = adata,
    n_latent = as.integer(x = ndims),
    n_layers = as.integer(x = nlayers),
    gene_likelihood = gene_likelihood
  )
  model$train(max_epochs = max_epochs)

  # extract the latent representation of the merged data
  latent <- model$get_latent_representation()
  latent <- as.matrix(latent)
  # pull the cell identifiers back out of the `AnnData` instance
  # in case anything was sorted under the hood
  rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
  # prepend the latent space dimensions with `new.reduction` to
  # give the features more readable names
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))

  # build a `DimReduc` instance
  suppressWarnings(
    latent.dr <- CreateDimReducObject(
      embeddings = latent, 
      key = new.reduction
    )
  )
  # to make it easier to add the reduction into a `Seurat` instance
  # we'll wrap it up in a named list
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction

  return(output.list)
}

attr(x = scVIIntegration, which = "Seurat.method") <- "integration"


#' Builds a data.frame with batch identifiers to use when integrating
#' \code{object}. For \code{StdAssays}, batches are split by layer.
#'
#' Internal - essentially the same as \code{Seurat:::CreateIntegrationGroups} 
#' except that it does not take in a `scale.layer` param.
#'
#' @noRd
#'
#' @param object A \code{StdAssay} instance.
#' @param layers Layers in \code{object} to integrate.
#'
#' @return A dataframe indexed on the cell identifiers from \code{object} - 
#' the dataframe contains a single column, "batch", indicating the layer/batch each cell is from
.FindBatches <- function(object, layers) {
  # build a LogMap indicating which layer each cell is from
  layer.masks <- slot(object, name = "cells")[, layers]
  # get a named vector mapping each cell to its respective layer
  layer.map <- labels(
      layer.masks,
      values = Cells(object, layer = layers)
  )
  # wrap the vector up in a data.frame
  batch.df <- as.data.frame(layer.map)
  names(batch.df) <- "batch"
  
  return(batch.df)
}


#' Builds a data.frame with batch identifiers to use when integrating
#' \code{object}. For \code{SCTAssay}s, batches are split using their
#' model identifiers. 
#'
#' Internal - essentially the same as \code{Seurat:::CreateIntegrationGroups} 
#' except that it does not take in a `scale.layer` param.
#'
#' @noRd
#'
#' @param object A \code{SCTAssay} or \code{StdAssays} instance.
#' @param layers Layers in \code{object} to integrate.
#'
#' @return A dataframe indexed on the cell identifiers from \code{object} - 
#' the dataframe contains a single column, "batch", indicating the layer/batch each cell is from
.FindSCTBatches <- function(object) {
  # build an empty data.frame indexed
  # on the cell identifiers from `object`
  batch.df <- SeuratObject::EmptyDF(n = ncol(object))
  row.names(batch.df) <- Cells(object)
  # for each
  for (sct.model in levels(object)) {
      cell.identifiers <- Cells(object, layer = sct.model)
      batch.df[cell.identifiers, "batch"] <- sct.model
  }
  return(batch.df)
}






