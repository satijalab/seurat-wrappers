#' @title Add Title Here
#'
#' @description Add description here.
#' 
#' @param object Describe input argument "object" here.
#'
#' @param k Describe input argument "k" here.
#'
#' @param assay Describe input argument "assay" here.
#'
#' @param features Describe input argument "features" here.
#'
#' @param reduction.name Describe input argument "reduction.name" here.
#'
#' @param reduction.key Describe input argument "reduction.key" here.
#'
#' @param verbose Describe input argument "verbose" here.
#' 
#' @param \dots Describe "..." here.
#'
#' @return Describe the return value here.
#'
#' Clarify that, unfortunately, "factors" and "loadings" in fastTopics
#' mean the opposite of what they mean in Seurat.
#' 
#' @seealso \code{\link[fastTopics]{fit_topic_model}}
#'
#' @importFrom stats prcomp
#' @importFrom Matrix t
#' 
#' @export
#'
FitTopicModel <- function (object, k = 3, assay = NULL, features = NULL,
                           reduction.name = "multinom_topic_model",
                           reduction.key = "k_", verbose = TRUE, ...) {

  # Check the input arguments, and that fastTopics is installed.
  CheckPackage(package = "stephenslab/fastTopics")
  if (!inherits(object,"Seurat"))
    stop("\"object\" must be a Seurat object",call. = FALSE)

  # Get the n x m counts matrix, where n is the number of samples
  # (cells) and m is the number of selected genes.
  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay
  features <- features %||% VariableFeatures(object)
  X <- GetAssayData(object,"counts")
  if (length(features) == 0)
    features <- rownames(X)
  else
    features <- intersect(features,rownames(X))
  X <- X[features,]
  X <- t(X)

  # TO DO: Remove any columns that are all zeros.
  
  # Fit the multinomial topic model using fastTopics.
  fit <- fit_topic_model(X,k,verbose = verbose,...)
  class(fit) <- c("list","multinom_topic_model_fit")
  
  # Retrieve the factors matrix (n x k) and loadings matrix (m x k).
  embeddings <- fit$L
  loadings <- fit$F
  colnames(embeddings) <- paste0(reduction.key,1:k)
  colnames(loadings) <- paste0(reduction.key,1:k)

  # Add the topic model fit to the Seurat object.
  object[[reduction.name]] <-
    CreateDimReducObject(embeddings,loadings,assay = assay,key = reduction.key,
                         global = TRUE,misc = fit)

  # Add a PCA dimension reduction from the mixture proportions.
  out <- prcomp(fit$L)
  colnames(out$x) <- paste0("PC_",1:k)
  colnames(out$rotation) <- paste0("PC_",1:k)
  object[["pca"]] <- 
    CreateDimReducObject(out$x[,-k],out$rotation[,-k],assay = assay,
                         key = "PC_",global = TRUE)
  
  return(LogSeuratCommand(object))
}
