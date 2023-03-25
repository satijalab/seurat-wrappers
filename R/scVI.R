#' @include internal.R
#'
NULL

#' Run scVI in seurat5
#' @param object A merged Seurat object
#' @param groups A one-column data frame with grouping information
#' @param features features to use 
#' @param layers Layers to use 
#' @param conda_env conda environment to run scVI
#' @param new.reduction Name to store resulting DimReduc object as
#' @param ... Arguments passed to other methods
#'
#' @export
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
#' obj <- IntegrateLayers(object = obj, method = scVIIntegration, new.reduction = 'integrated.scvi',
#'                        conda_env = '../miniconda3/envs/scvi-env', verbose = FALSE)
#' }
#'
#' # Alternatively, we can integrate SCTransformed data 
#' obj <- SCTransform(object = obj)
#' obj <- IntegrateLayers(object = obj, method = scVIIntegration, 
#'           orig.reduction = "pca", new.reduction = 'integrated.scvi', 
#'            assay = "SCT", conda_env = '../miniconda3/envs/scvi-env', verbose = FALSE)
#'
#' @seealso \href{https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scvi_in_R.html}{scVI}
#' 
#' @return A Seurat object with embeddings and loadings

scVIIntegration <- function(
    object,
    groups = NULL,
    features = NULL,
    layers = 'counts',
    conda_env = NULL,
    new.reduction = 'integrated.dr',
    ndims = 10,
    max_epochs = 200,
    ...){
  reticulate::use_condaenv(conda_env, required = TRUE)
  sc <-  reticulate::import('scanpy', convert = FALSE)
  scvi <-  reticulate::import('scvi', convert = FALSE)
  anndata <-  reticulate::import('anndata', convert = FALSE)
  scipy <-  reticulate::import('scipy', convert = FALSE)
  adata <- sc$AnnData(
    X   = scipy$sparse$csr_matrix(Matrix::t(LayerData(object, layers = layers)[features ,]) ), #scVI requires raw counts
    obs = groups,
    var = object[][features,]
  )

  scvi$model$SCVI$setup_anndata(adata,  batch_key = 'group')
  model = scvi$model$SCVI(adata = adata, n_latent = as.integer(x = ndims))
  model$train(max_epochs = as.integer(x = max_epochs))

  latent = model$get_latent_representation()
  latent <- as.matrix(latent)
  rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))
  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent, key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  return(output.list)
}

attr(x = scVIIntegration, which = 'Seurat.method') <- 'integration'