#' @include internal.R
#'
NULL

#' Run fastMNN in Seurat 5
#'
#' @param object A merged seurat object
#' @param groups A one-column data frame with grouping information
#' @param layers Layers to use 
#' @param assay Assay to use, defaults to the default assay of the first object
#' @param features Either a list of features to use when calculating batch
#' correction, or a number (2000 by default) of variable features to select.
#' @param reduction.name Name to store resulting DimReduc object as
#' @param reduction.key Key for resulting DimReduc
#' @param reconstructed.assay Name for the assay containing the low-rank
#' reconstruction of the expression matrix.
#' @param verbose Print messages
#' @param ... Extra parameters passed to \code{\link[batchelor]{fastMNN}}
#'
#' @return A Seurat object merged from the objects in \code{object.list} and a
#' new DimReduc of name \code{reduction.name} (key set to \code{reduction.key})
#' with corrected embeddings matrix as well as the rotation matrix used for the
#' PCA stored in the feature loadings slot. Also returns an expression matrix
#' reconstructed from the low-rank approximation in the
#' \code{reconstructed.assay} assay; all other metadata info
#' \code{\link[batchelor]{fastMNN}} is stored in the \code{tool} slot,
#' accessible with \code{\link[Seurat]{Tool}}
#'
#' @importFrom Seurat DefaultAssay DefaultAssay<- SelectIntegrationFeatures VariableFeatures VariableFeatures<-
#' as.SingleCellExperiment CreateDimReducObject Tool<- LogSeuratCommand
#' @importFrom rlang check_installed
#'
#' @export
#'
#' @examples 
#' \dontrun{
#' # Preprocessing
#' obj <- LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#' 
#' # After preprocessing, we integrate layers: 
#' obj <- IntegrateLayers(object = obj, method = FastMNNIntegration, 
#'   new.reduction = 'integrated.mnn', verbose = FALSE)
#'   
#' # We can also add parameters specific to FastMNN. 
#' # Here we set `k` to specify the number of nearest neighbors to use when identifying MNNs: 
#' obj <- IntegrateLayers(object = obj, method = FastMNNIntegration, 
#'   new.reduction = 'integrated.mnn', k = 15, verbose = FALSE)
#' }
#' 
#' @seealso \code{\link[batchelor]{fastMNN}} \code{\link[Seurat]{Tool}}
#' 
FastMNNIntegration <- function(
    object,
    assay = NULL,
    orig = NULL,
    groups = NULL, 
    layers = NULL,
    scale.layer = NULL,
    features = 2000,
    new.reduction = "integrated.mnn",
    reduction.key = "mnn_",
    reconstructed.assay = "mnn.reconstructed",
    verbose = TRUE,
    ...
) {
  check_installed(
    pkg = "batchelor",
    reason = "for running integration with mnnCorrect"
  )
  object <- CreateSeuratObject(object)
  if (is.numeric(x = features)) {
    if (verbose) {
      message(paste("Computing", features, "integration features"))
    }
    features <- SelectIntegrationFeatures5(object = object, features = features)
  }
  layers <- layers %||% Layers(object, search = 'data')
  if (verbose) {
    message("Converting layers to SingleCellExperiment")
  }
  objects.sce <- lapply(
    X = layers,
    FUN = function(x, f) {
      return(as.SingleCellExperiment(
        x = subset(x = object,
                   features = f,
                   cells = colnames(LayerData(object, layer = x)))
            )
        )
    },
    f = features
  )
  if (verbose) {
    message("Running fastMNN")
  }
  out <- do.call(
    what = batchelor::fastMNN,
    args = c(
      objects.sce,
      list(...)
    )
  )
  colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key, 1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
  reduction <- CreateDimReducObject(
    embeddings = SingleCellExperiment::reducedDim(x = out),
    loadings = as.matrix(SingleCellExperiment::rowData(x = out)),
    assay = DefaultAssay(object = object),
    key = reduction.key
  )
  # Add reconstructed matrix (gene x cell)
  reconstructed_assay <- CreateAssayObject(
    data = as(object = SummarizedExperiment::assay(x = out), Class = "sparseMatrix"),
  )
  # Add variable features
  VariableFeatures(object = reconstructed_assay) <- features
  #Tool(object = object) <- S4Vectors::metadata(x = out)
  #object <- LogSeuratCommand(object = object)
  output.list <- list(reduction, reconstructed_assay)
  names(output.list) <- c(new.reduction, reconstructed.assay)
  return(output.list)
}

attr(x = FastMNNIntegration, which = 'Seurat.method') <- 'integration'