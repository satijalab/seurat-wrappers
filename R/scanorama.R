#' @include internal.R
#'
NULL

#' Run scanorama in \link[Seurat:Assay5] through \link[Seurat:IntegrateLayers]
#'
#' @param object An \code{\link[SeuratObject:Assay5]} object
#' @param orig A \link[SeuratObject:DimReduc]{dimensional reduction} to correct
#' @param layers Layers to look for. Either a vector of layer names or a single
#' string (e.g. 'counts' or 'data')
#' @param features A vector of features to use for integration.
#' All features used by default
#' @param conda_env Path to conda environment to run scanorama (should also
#' contain the scipy python module)
#' @param new.reduction Name to store resulting
#' \link[SeuratObject:DimReduc]{dimensional reduction} object. \code{NULL} or
#' \code{FALSE} disables the dimension reduction computation
#' @param reduction.key Key for resulting
#' \link[SeuratObject:DimReduc]{dimensional reduction} object. Ignored if
#' \code{reduction.key} is \code{NULL} or \code{FALSE}
#' @param reconstructed.assay Name for the \code{assay} containing the corrected
#' expression matrix of dimension #\code{features} x #cells. \code{NULL} or
#' \code{FALSE} disables the corrected expression matrix computation
#' @param ndims.out Number of dimensions for \code{new.reduction} output.
#' Scanorama specific argument
#' @param return_dense Whether scanorama returns \code{numpy.ndarray} matrices
#' instead of \code{scipy.sparse.csr_matrix}. Scanorama specific argument
#' @param batch_size Used in the alignment vector computation. Higer values
#' increase memory burden.Scanorama specific argument
#' @param approx Whether scanorama approximates nearest neighbors computation.
#' Scanorama specific argument
#' @param sigma  Correction smoothing parameter on gaussian kernel. Scanorama
#' specific argument
#' @param alpha  Minimum alignment score. Scanorama specific argument
#' @param knn  Number of nearest neighbours used for matching. Scanorama
#' specific argument
#' @param hvg.scanorama A positive to turn scanorama's internal HVG selection on.
#' Disabled by default. Scanorama specific argument
#' @param union.features By default, scanorama uses the intersection of features
#' to perform integration. Set this parameter to \code{TRUE} to use the union.
#' Discouraged. Scanorama specific argument
#' @param sketch Turns on sketching-based acceleration by first downsampling the
#' datasets. See Hie et al., Cell Systems (2019). Disabled by default. Ignored
#' when \code{reconstructed.assay} is enabled. Scanorama specific argument
#' @param sketch_method A skething method to apply to the data. Either
#' 'geosketch' (default) or 'uniform'. Ignored when \code{reconstructed.assay}
#' is enabled or when \code{sketch} is FALSE. Scanorama specific argument
#' @param sketch_max Downsampling cutoff. Ignored when sketching disabled.
#' Scanorama specific argument
#' @param seed.use Value of scanorama's internal seed. Set to \code{NULL} to
#' disable
#' @param verbose Print messages
#' @param ... Ignored
#'
#' @return A list containing at least one of:
#' \itemize{
#'   \item a new DimReduc of name \code{reduction.name} (key set to
#'   \code{reduction.key}) with corrected embeddings matrix of \code{ndims.out}.
#'   \item a new Assay of name \code{reconstructed.assay} with corrected counts
#'   for each \code{features} (or less if \code{hvg} is set to a lower integer).
#' }
#'
#' @importFrom Seurat DefaultAssay DefaultAssay<- VariableFeatures VariableFeatures<-
#' CreateDimReducObject CreateAssayObject Cells
#' @importFrom SeuratObject Layers LayerData Features
#'
#' @export
#'@note This function requires the
#' \href{https://github.com/brianhie/scanorama}{\pkg{Scanorama}} package
#' to be installed (along with \pkg{scipy})
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
#' # After preprocessing, we integrate layers:
#' obj <- IntegrateLayers(object = obj, method = ScanoramaIntegration)
#'
#' # To disable feature expression matrix correction:
#' obj <- IntegrateLayers(object = obj, method = ScanoramaIntegration,
#'                        reconstructed.assay = NULL)
#' }
#'
#' @references Hie, B., Bryson, B. & Berger, B. Efficient integration of
#' heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol 37,
#' 685-691 (2019). \href{https://doi.org/10.1038/s41587-019-0113-3}{DOI}
#' 
#' @seealso \code{\link[Seurat]{Tool}}

ScanoramaIntegration <- function(
    object,
    orig,
    layers = NULL,
    features = NULL,
    scale.layer = 'scale.data',
    conda_env = NULL,
    new.reduction = "integrated.scanorama", 
    reduction.key = "scanorama_",
    reconstructed.assay = "scanorama.reconstructed",
    ndims.out = 100L,
    return_dense = TRUE,
    batch_size = 5e3,
    approx = TRUE,
    sigma = 15L,
    alpha = .1,
    knn = 20L,
    hvg.scanorama = NULL,
    union.features = FALSE,
    sketch = FALSE,
    sketch_method = c('geosketch', 'uniform'),
    sketch_max = 1e4,
    seed.use = 42L,
    verbose = TRUE,
    ...) {
  sketch_method <- match.arg(arg = sketch_method)
  
  reticulate::use_condaenv(conda_env, required = TRUE)
  scano <- reticulate::import('scanorama', convert = FALSE)
  scipy <- reticulate::import('scipy', convert = FALSE)
  
  layers <- layers %||% 'counts'
  layers <- suppressWarnings({Layers(object = object, search = layers)}) %||%
    {
      Layers(object = object, search = 'counts')
      warning(paste(sQuote(layers), "not found in layers. Using 'counts'"),
              call. = FALSE, immediate. = TRUE)
    }
  
  data.list <- lapply(X = layers, FUN = function(layer) {
    layer.data <- LayerData(object = object, layer = layer, features = features)
    scipy$sparse$csr_matrix(Matrix::t(layer.data))
  })
  features <- features %||% Features(object = object, layer = layers, simplify = F)
  if (! is.list(features)) {    #HVGs
    features <- list(features)[rep(1,length(data.list))]
  }
  hvg.scanorama <- hvg.scanorama %iff% as.integer(hvg.scanorama)
  
  new.reduction <- new.reduction %||% FALSE
  reconstructed.assay <- reconstructed.assay %||% FALSE
  dimred <- ! isFALSE(new.reduction)
  cntcor <- ! isFALSE(reconstructed.assay)
  if (! (dimred | cntcor)) {
    abort(message = "missing both 'new.reduction' and 'reconstructed.assay'")
  }
  if (! cntcor) {
    out <- scano$integrate(reticulate::r_to_py(data.list),
                           reticulate::r_to_py(features),
                           dimred = as.integer(ndims.out),
                           batch_size = as.integer(batch_size),
                           approx = approx,
                           sigma = sigma,
                           alpha = alpha,
                           knn = as.integer(knn),
                           hvg = hvg.scanorama,
                           union = union.features,
                           sketch = sketch,
                           sketch_method = sketch_method,
                           sketch_max = as.integer(sketch_max),
                           seed = as.integer(seed.use),
                           verbose = verbose + 0)
  } else {
    out <- scano$correct(reticulate::r_to_py(data.list),
                         reticulate::r_to_py(features),
                         return_dimred = dimred,
                         return_dense=return_dense,
                         dimred = as.integer(ndims.out),
                         batch_size = as.integer(batch_size),
                         approx = approx,
                         sigma = sigma,
                         alpha = alpha,
                         knn = as.integer(knn),
                         hvg = hvg.scanorama,
                         union = union.features,
                         seed = as.integer(seed.use),
                         verbose = verbose + 0)
    
    cnts.scano <- do.call(rbind, reticulate::py_to_r(out[0 + dimred]))
    rownames(cnts.scano) <- Cells(object, layer = layers)
    colnames(cnts.scano) <- reticulate::py_to_r(out[1 + dimred])
  }
  
  features <- unique(x = reticulate::py_to_r(out[-1]))
  output.list <- list()
  if (dimred) {
    embeddings <- do.call(rbind, reticulate::py_to_r(out[0]))
    rownames(embeddings) <- Cells(object, layer = layers)
    colnames(embeddings) <- paste0(reduction.key, 1:min(ndims.out, hvg.scanorama))
    
    output.list[[new.reduction]] <- CreateDimReducObject(
      embeddings = embeddings,
      assay = DefaultAssay(object = orig),
      key = reduction.key
    )
  }
  if (cntcor) {
    output.list[[reconstructed.assay]] <- CreateAssayObject(
      data = Matrix::t( as(cnts.scano, "sparseMatrix") )
    )
   VariableFeatures(object = output.list[[reconstructed.assay]]) <- features
  }
  return(output.list)
}

attr(x = ScanoramaIntegration, which = 'Seurat.method') <- 'integration'
