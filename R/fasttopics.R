#' @title Fit a Multinomial Topic Model Using fastTopics
#'
#' @description Fits a multinomial topic model to the count data,
#'   hiding most of the complexities of model fitting. The default
#'   optimization settings used here are intended to work well in a wide
#'   range of data sets, although some fine-tuning may be needed for
#'   more difficult cases.
#'
#' @details The topic model can be understood as a dimensionality
#' reduction, in which the mixture proportions matrix, \code{L}, an n
#' x k matrix, defines a projection of the cells onto a
#' (k–1)-dimension space (n is the number of cells and k is the number
#' of topics). Therefore, we use the mixture proportions to define the
#' \dQuote{cell embedding} in the \code{DimReduc} object.
#'
#' There is a confusing bit of terminology here, because "loadings"
#' means something different depending on the context In fastTopics,
#' the \code{L} matrix is sometimes referred to as the "loadings
#' matrix" (this is the convention used in factor analysis), but it is
#' \emph{not} the same as the "feature loadings" matrix used here,
#' which adopts the convention used for principal component
#' analysis. To avoid confusion, we only refer to \code{L} as the
#' mixture proportions matrix.
#' 
#' See \code{\link[fastTopics]{fit_topic_model}} for more details.
#' 
#' @param object A Seurat object.
#'
#' @param k The number of topics. Must be 2 or more.
#'
#' @param assay Name of assay to use; defaults to the default
#'   assay of the object.
#'
#' @param features A list of features to use for fitting the model. If
#'   \code{features = NULL}, \emph{all} features are used; see
#'   \code{\link[Seurat]{VariableFeatures}}. NOTE: EXPAND ON THIS.
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
#' @param \dots Additional arguments passed to \code{fit_topic_model};
#'   see \code{\link[fastTopics]{fit_topic_model}} for details.
#'
#' @return A Seurat object with the multinomial topic model fit stored
#' as a \code{\link[Seurat]{DimReduc}} object. The cell embeddings
#' are the mixture proportions, stored as an n x k matrix in the
#' \code{cell.embeddings} slot, where n is the number of cells, and k
#' is the number of topics.  This is the same as the \code{L} matrix
#' output in \code{fit_topic_model}.
#'
#' The feature loadings (the \code{feature.loadings} slot) are an m x
#' k matrix in which column i is the vector of relative gene
#' expression levels in topic i. This is the same as the \code{F}
#' output in \code{fit_topic_model}.
#'
#' Apply \code{\link[Seurat]{Misc}} to the \code{DimReduce} object to
#' access the \dQuote{"multinom_topic_model_fit"} object outputted by
#' \code{\link[fastTopics]{fit_topic_model}}, which contains more
#' information about the model fit; see the example for an
#' illustration of how to do this, and see
#' \code{\link[fastTopics]{fit_topic_model}} more information about
#' the "multinom_topic_model_fit" object.
#'
#' An additional PCA dimension reduction is provided to, for example,
#' quickly visualize the cells using
#' \code{\link[Seurat]{DimPlot}}. The principal components are
#' computed from the topic model mixture proportions. However, this is
#' for convenience only, and we recommend extracting the
#' \dQuote{"multinom_topic_model_fit"} object and dedicated
#' visualization tools such \code{\link[fastTopics]{structure_plot}}
#' that are provided in the fastTopics package.
#' 
#' NOTE: There are only k - 1 PCs.
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
#' library(Seurat)
#' library(fastTopics)
#' set.seed(1)
#'
#' # Load the PBMC data.
#' data(pbmc_small)
#'
#' # Fit the multinomial topic model to the raw UMI count data; no
#' # pre-processing is needed. 
#' pbmc_small <- FitTopicModel(pbmc_small,k = 3)
#'
#' # This plot shows the cells projected onto the 2 principal
#' # components (PCs) of the topic model mixture proportions.
#' DimPlot(pbmc_small,reduction = "pca_topics")
#'
#' # Compare against PCA performed on the transformed count data.
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
#' @importFrom Matrix colSums
#' @importFrom Matrix t
#' @importFrom fastTopics fit_topic_model
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
  X <- GetAssayData(object,"counts")
  if (is.null(features))
    features <- rownames(X)
  else
    features <- intersect(features,rownames(X))
  X <- X[features,]
  X <- t(X)

  # Remove all-zero columns.
  i        <- which(colSums(X > 0) >= 1)
  features <- features[i]
  X        <- X[,i]
  
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
  colnames(out$x) <- paste0("TOPICPC_",1:k)
  colnames(out$rotation) <- paste0("TOPICPC_",1:k)
  object[["pca_topics"]] <- 
    CreateDimReducObject(out$x[,-k],out$rotation[,-k],assay = assay,
                         key = "TOPICPC_",global = TRUE)

  # Output the updated Seurat object.
  return(LogSeuratCommand(object))
}
