#' @title Fit a Multinomial Topic Model Using fastTopics
#'
#' @description Add description here.
#' 
#' @param object A Seurat object.
#'
#' @param k The number of topics. Must be 2 or more.
#'
#' @param assay Name of assay to use; defaults to the default
#'   assay of the object.
#'
#' @param features A list of features to use for fitting the model. If
#'   \code{features = NULL}, all variable features; see
#'   \code{\link[Seurat]{VariableFeatures}}.
#'
#' @param reduction.name Name of the outputted reduction.
#'
#' @param reduction.key Key for the outputted reduction.
#'
#' @param verbose When \code{verbose = TRUE}, information about the
#'   progress of the model fitting is printed to the console. See
#'   \code{\link[fastTopics]{fit_poisson_nmf}} for an explanation of the
#'   output.
#' 
#' @param \dots \dots Additional arguments passed to
#'   \code{fit_topic_model}; see
#'   \code{\link[fastTopics]{fit_topic_model}} for details
#'
#' @return Describe the return value here.
#'
#' Clarify that, unfortunately, "factors" and "loadings" in fastTopics
#' mean the opposite of what they mean in Seurat.
#'
#' @author Peter Carbonetto
#'
#' @references
#' Dey, K. K., Hsiao, C. J. and Stephens, M. (2017). Visualizing the
#' structure of RNA-seq expression data using grade of membership
#' models. \emph{PLoS Genetics} \bold{13}, e1006599.
#'
#' @seealso \code{\link[fastTopics]{fit_topic_model}}
#'
#' @examples
#' set.seed(1)
#'
#' # Load the PBMC data.
#' data(pbmc_small)
#'
#' # Fit the multinomial topic model to the raw UMI count data; no
#' # pre-processing is needed.
#' pbmc_small <- FitTopicModel(pbmc_small,k = 3)
#'
#' # This plot shows the cells projected onto the top 2 principal
#' # components (PCs) of the topic model mixture proportions.
#' Idents(pbmc_small) <- pbmc_small$letter.idents
#' DimPlot(pbmc_small,reduction = "pca")
#'
#' # Once fitted topic model is extracted, many functions from the
#' # fastTopics package can be used. For example, the Structure plot
#' # provides an evocative visual summary of the estimated mixture
#' # proportions for each cell.
#' fit <- Misc(Reductions(pbmc_small,"multinom_topic_model"))
#' structure_plot(fit,grouping = Idents(pbmc_small),gap = 5)
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

  # Add a PCA dimension reduction calculated from the mixture
  # proportions.
  out <- prcomp(fit$L)
  colnames(out$x) <- paste0("PC_",1:k)
  colnames(out$rotation) <- paste0("PC_",1:k)
  object[["pca"]] <- 
    CreateDimReducObject(out$x[,-k],out$rotation[,-k],assay = assay,
                         key = "PC_",global = TRUE)

  # Output the updated Seurat object.
  return(LogSeuratCommand(object))
}
