#' @include internal.R
#'
NULL

#' Run GLMPCA
#'
#' @param object A Seurat object
#' @param L The number of dimensions to return (defaults to 5)
#' @param assay Assay to use, defaults to the default assay
#' @param features A list of features to use when performing GLM-PCA. If null, defaults to variable features.
#' @param reduction.name Name to store resulting DimReduc object as. Defaults to glmpca
#' @param reduction.key Key for resulting DimReduc. Defaults to GLMPC_
#' @param ... Extra parameters passed to \code{\link[glmpca]{glmpca}}
#'
#' @return A Seurat object containing the output of GLMPCA stored as a DimReduc object.
#' @importFrom Seurat DefaultAssay DefaultAssay<- CreateDimReducObject Tool<- LogSeuratCommand
#'
#' @author Will Townes
#' @references Townes, W., Hicks, SC, Aryee, MJ, Irizarry, RA. (2019). "Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model."
#' Genome Biology.
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' pbmc_small <- RunGLMPCA(pbmc_small)
#' DimPlot(pbmc_small, redunction = 'glmpca')
#' }
#'
#' @export
#'
RunGLMPCA <- function(
  object,
  L = 5,
  assay = NULL,
  features = NULL,
  reduction.name = 'glmpca',
  reduction.key = 'GLMPC_',
  verbose = TRUE,
  ...
) {
  CheckPackage(package = 'glmpca', repository = 'CRAN')
  if (!inherits(x = object,what = 'Seurat')) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object)
  data <- GetAssayData(object = object, slot = 'counts')
  features <- intersect(x = features, y = rownames(x = data))
  if (length(x = features) == 0) {
    stop("Please specify a subset of features for GLM-PCA")
  }
  data <- data[features, ]
  glmpca_results <- glmpca:::glmpca(Y = data, L = L, ...)
  glmpca_dimnames <- paste0(reduction.key, 1:L)
  factors<-as.matrix(glmpca_results$factors)
  loadings<-as.matrix(glmpca_results$loadings)
  colnames(x = factors) <- glmpca_dimnames
  colnames(x = loadings) <- glmpca_dimnames
  factors_l2_norm <- sqrt(colSums(factors^2))
  #strip S3 class "glmpca" to enable it to pass validObject()
  class(glmpca_results)<-NULL
  #save memory by removing factors and loadings since they are stored separately
  glmpca_results$factors<-glmpca_results$loadings<-NULL
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = factors,
    key = reduction.key,
    loadings = loadings,
    stdev = factors_l2_norm,
    assay = assay,
    global = TRUE,
    misc = glmpca_results
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}
