#' @title Fit a Non-negative Matrix Factorization Using fastTopics
#'
#' @description Approximate the raw count data \code{X} by the
#'   non-negative matrix factorization \code{tcrossprod(L,F)}, in which
#'   the quality of the approximation is measured by a
#'   \dQuote{divergence} criterion; equivalently, optimize the
#'   likelihood under a Poisson model of the count data, \code{X}, in
#'   which the Poisson rates are given by \code{tcrossprod(L,F)}.
#'
#' @param object A Seurat object. If this Seurat object contains a
#'   Poisson NMF dimensionality reduction object (\dQuote{DimReduc}),
#'   the model fitting will be initialized to this previously fitted
#'   Poisson NMF.
#' 
#' @param k An integer 2 or greater giving the matrix rank. This
#'   argument should only be specified if the initial fit (\code{fit0})
#'   is not provided.
#'
#' @param assay Name of assay to use; defaults to the default
#'   assay of the object.
#'
#' @param features A list of features to use for fitting the model. If
#'   \code{features = NULL}, \emph{all} features are used; in
#'   particular, \code{\link[Seurat]{VariableFeatures}} is \emph{not}
#'   used to pre-select features.
#'
#' @param reduction.name Name of the outputted reduction.
#'
#' @param reduction.key Key for the outputted reduction.
#'
#' @param reduction.pca.name Name of the outputted PCA reduction.
#'
#' @param reduction.pca.key Key for the outputted PCA reduction.
#' 
#' @param numiter The number of updates of the factors and loadings to
#'   perform.
#'
#' @param method The method to use for updating the factors and
#'   loadings. See \code{\link[fastTopics]{fit_poisson_nmf}} for
#'   details.
#'
#' @param init.method The method used to initialize the factors and
#'   loadings. See \code{\link[fastTopics]{fit_poisson_nmf}} for
#'   details.
#' 
#' @param control A list of parameters controlling the behaviour of
#'   the optimization algorithm. See
#'   \code{\link[fastTopics]{fit_poisson_nmf}} for details.
#'
#' @param verbose When \code{verbose = "progressbar"} or \code{verbose
#'   = "detailed"}, information about the progress of the model fitting
#'   is printed to the console. See \code{\link[fastTopics]{fit_poisson_nmf}}
#'   for more information.
#' 
#' @param \dots Additional arguments passed to \code{fit_poisson_nmf}.
#'
#' @return A Seurat object, in which the Poisson NMF fit is stored as
#' a Seurat \code{\link[Seurat]{DimReduc}} object. The cell embeddings
#' (stored in the \code{cell.embeddings} slot) are the loadings; this
#' is the n x k matrix \code{L} outputted by \code{fit_poisson_nmf} (n
#' is the number of cells and k is the dimension of the matrix
#' factorization).
#' 
#' The feature loadings (stored in the \code{feature.loadings} slot)
#' are the m x k factors matrix \code{F} outputted by
#' \code{fit_poisson_nmf}, where m is the number of features.
#'
#' Apply \code{\link[Seurat]{Misc}} to the \code{DimReduce} object to
#' access the \dQuote{"poisson_nmf_fit"} object outputted by
#' \code{\link[fastTopics]{fit_poisson_nmf}}; see the example for an
#' illustration.
#'
#' An additional PCA dimension reduction computed from the k loadings
#' is provided.
#' 
#' @author Peter Carbonetto
#'
#' @references
#'   Lee, D. D. and Seung, H. S. (2001). Algorithms for non-negative
#'   matrix factorization. In \emph{Advances in Neural Information
#'   Processing Systems} \bold{13}, 556–562.
#'
#' @seealso \code{\link{FitTopicModel}},
#'   \code{\link[fastTopics]{fit_poisson_nmf}}
#'
#' @examples
#' library(Seurat)
#' library(fastTopics)
#' set.seed(1)
#' data(pbmc_small)
#'
#' # Fit the non-negative matrix factorization to the raw UMI count
#' # data; no pre-processing or pre-selection of genes is needed.
#' pbmc_small <- FitPoissonNMF(pbmc_small,k = 3,numiter = 20)
#'
#' # Improve the fit by running another 20 updates.
#' pbmc_small <- FitPoissonNMF(pbmc_small,k = 3,numiter = 20)
#'
#' # This plot shows the cells projected onto the 2 principal
#' # components (PCs) of the topic mixture proportions.
#' DimPlot(pbmc_small,reduction = "pca_nmf")
#'
#' # Compare this against the top two PCs of the transformed count
#' # data.
#' DimPlot(pbmc_small,reduction = "pca")
#'
#' # Extract the non-negative matrix factorization.
#' fit <- Misc(Reductions(pbmc_small,"poisson_nmf"))
#' summary(fit)
#' 
#' @importFrom fastTopics fit_poisson_nmf
#' 
#' @export
#' 
FitPoissonNMF <- function (object, k, assay = NULL, features = NULL,
                           reduction.name = "poisson_nmf",
                           reduction.key = "k_",
                           reduction.pca.name = "pca_nmf",
                           reduction.pca.key = "NMFPC_", numiter = 100,
                           method = c("scd", "em", "mu", "ccd"),
                           init.method = c("topicscore", "random"),
                           control = list(),
                           verbose = c("progressbar", "detailed", "none"),
                           ...) {

  # Check the input arguments, and that fastTopics is installed.
  CheckPackage(package = "stephenslab/fastTopics")
  if (!inherits(object,"Seurat"))
    stop("\"object\" must be a Seurat object",call. = FALSE)

  # Check and progress input argument "verbose".
  verbose <- match.arg(verbose)
  
  # Get the n x m counts matrix, where n is the number of samples
  # (cells) and m is the number of selected features.
  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay
  X <- prepare_counts_fasttopics(object,features)
  features <- colnames(X)

  # Fit the Poisson non-negative matrix factorization using
  # fastTopics. If Seurat object has an existing "poisson_nmf"
  # reduction, use this to initialize the model fitting.
  if (is.element("poisson_nmf",Reductions(object))) {
    fit0 <- Misc(Reductions(object,"poisson_nmf"))
    fit <- fit_poisson_nmf(X,fit0 = fit0,numiter = numiter,method = method,
                           init.method = init.method,control = control,
                           verbose = verbose,...)
  } else
    fit <- fit_poisson_nmf(X,k,numiter = numiter,method = method,
                           init.method = init.method,control = control,
                           verbose = verbose,...)
  class(fit) <- c("list","poisson_nmf_fit")
  
  # Retrieve the factors matrix (n x k) and loadings matrix (m x k).
  embeddings <- fit$L
  loadings <- fit$F
  colnames(embeddings) <- paste0(reduction.key,1:k)
  colnames(loadings) <- paste0(reduction.key,1:k)

  # Add the topic model fit to the Seurat object.
  object[[reduction.name]] <-
    CreateDimReducObject(embeddings,loadings,assay = assay,
                         key = reduction.key,global = TRUE,
                         misc = fit)

  # Add a PCA dimension reduction calculated from the mixture
  # proportions.
  object[[reduction.pca.name]] <-
    pca_from_loadings_fasttopics(fit,assay,reduction.pca.key)

  # Output the updated Seurat object.
  return(LogSeuratCommand(object))
}

#' @title Fit a Multinomial Topic Model Using fastTopics
#'
#' @description Fits a multinomial topic model to the raw count data,
#'   hiding most of the complexities of model fitting. The default
#'   optimization settings used here are intended to work well in a wide
#'   range of data sets, although some fine-tuning may be needed for
#'   larger or more complex data sets.
#'
#' @details The topic model can be understood as a dimensionality
#' reduction, in which the mixture proportions matrix, \code{L}, an n
#' x k matrix, defines a projection of the cells onto a
#' (k–1)-dimension space (n is the number of cells and k is the number
#' of topics). The topic mixture proportions define the \dQuote{cell
#' embedding} in the \code{DimReduc} object.
#'
#' \emph{A warning about confusing terminology:} In fastTopics,
#' the \code{L} matrix is sometimes referred to as the "loadings
#' matrix" (because this is the convention used in factor
#' analysis). However, this is \emph{not} the same as the "feature
#' loadings" matrix used in Seurat, which adopts the convention used n
#' principal component analysis. To avoid confusion, we refer to
#' \code{L} as the mixture proportions matrix.
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
#'   \code{features = NULL}, \emph{all} features are used; in
#'   particular, \code{\link[Seurat]{VariableFeatures}} is \emph{not}
#'   used to pre-select features.
#'
#' @param reduction.name Name of the outputted reduction.
#'
#' @param reduction.key Key for the outputted reduction.
#'
#' @param reduction.pca.name Name of the outputted PCA reduction.
#'
#' @param reduction.pca.key Key for the outputted PCA reduction.
#' 
#' @param verbose When \code{verbose = "progressbar"} or \code{verbose
#'   = "detailed"}, information about the progress of the model fitting
#'   is printed to the console. See \code{\link[fastTopics]{fit_poisson_nmf}}
#'   for more information.
#' 
#' @param \dots Additional arguments passed to \code{fit_topic_model};
#'   see \code{\link[fastTopics]{fit_topic_model}} for details.
#'
#' @return A Seurat object, in which the multinomial topic model fit
#' is stored as a Seurat \code{\link[Seurat]{DimReduc}} object. The
#' cell embeddings (stored in the \code{cell.embeddings} slot) are the
#' mixture proportions; this is the n x k matrix \code{L} outputted by
#' \code{fit_topic_model} (n is the number of cells and k is the
#' number of topics).
#' 
#' The feature loadings (stored in the \code{feature.loadings} slot)
#' are an m x k matrix of relative expression levels, where m is the
#' number of features and k is the number of topics; this is the
#' matrix \code{F} outputted by \code{fit_topic_model}.
#'
#' Apply \code{\link[Seurat]{Misc}} to the \code{DimReduce} object to
#' access the \dQuote{"multinom_topic_model_fit"} object outputted by
#' \code{\link[fastTopics]{fit_topic_model}}, which contains more
#' information about the multinomial topic model fit; see the example
#' here for an illustration of how to access the model fit, and see
#' \code{\link[fastTopics]{fit_topic_model}} more information about
#' the "multinom_topic_model_fit" object.
#'
#' An additional PCA dimension reduction on \code{k-1} dimensions,
#' computed from the mixture proportions, is provided. This could be
#' useful, for example, to quickly visualize the cells using
#' \code{\link[Seurat]{DimPlot}}. The principal components are
#' computed from the topic mixture proportions. However, this
#' reduction is provided for convenience only, and we recommend to
#' extract the \dQuote{"multinom_topic_model_fit"} object and use the
#' dedicated visualization tools that are provided in the fastTopics
#' package such as \code{\link[fastTopics]{structure_plot}}.
#' 
#' @author Peter Carbonetto
#'
#' @references
#' Dey, K. K., Hsiao, C. J. and Stephens, M. (2017). Visualizing the
#' structure of RNA-seq expression data using grade of membership
#' models. \emph{PLoS Genetics} \bold{13}, e1006599.
#'
#' @seealso \code{\link{FitPoissonNMF}},
#'   \code{\link[fastTopics]{fit_topic_model}}
#'
#' @examples
#' library(Seurat)
#' library(fastTopics)
#' set.seed(1)
#' data(pbmc_small)
#'
#' # Fit the multinomial topic model to the raw UMI count data; no
#' # pre-processing or pre-selection of genes is needed.
#' pbmc_small <- FitTopicModel(pbmc_small,k = 3)
#'
#' # This plot shows the cells projected onto the 2 principal
#' # components (PCs) of the topic mixture proportions.
#' DimPlot(pbmc_small,reduction = "pca_topics")
#'
#' # Compare this against the top two PCs of the transformed count
#' # data.
#' DimPlot(pbmc_small,reduction = "pca")
#'
#' # Once fitted topic model is extracted, many functions from the
#' # fastTopics package can be used for analysis and visualization. For
#' # example, the Structure plot provides an evocative visual summary of
#' # the estimated mixture proportions for each cell.
#' fit <- Misc(Reductions(pbmc_small,"multinom_topic_model"))
#' structure_plot(fit,grouping = Idents(pbmc_small),gap = 5)
#'
#' @importFrom fastTopics fit_topic_model
#' 
#' @export
#'
FitTopicModel <- function (object, k = 3, assay = NULL, features = NULL,
                           reduction.name = "multinom_topic_model",
                           reduction.key = "k_",
                           reduction.pca.name = "pca_topics",
                           reduction.pca.key = "TOPICPC_",
                           verbose = c("progressbar", "detailed", "none"),
                           ...) {

  # Check the input arguments, and that fastTopics is installed.
  CheckPackage(package = "stephenslab/fastTopics")
  if (!inherits(object,"Seurat"))
    stop("\"object\" must be a Seurat object",call. = FALSE)

  # Check and progress input argument "verbose".
  verbose <- match.arg(verbose)
  
  # Get the n x m counts matrix, where n is the number of samples
  # (cells) and m is the number of selected features.
  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay
  X <- prepare_counts_fasttopics(object,features)
  features <- colnames(X)
  
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
    CreateDimReducObject(embeddings,loadings,assay = assay,
                         key = reduction.key,global = TRUE,
                         misc = fit)

  # Add a PCA dimension reduction calculated from the mixture
  # proportions.
  object[[reduction.pca.name]] <-
    pca_from_loadings_fasttopics(fit,assay,reduction.pca.key)

  # Output the updated Seurat object.
  return(LogSeuratCommand(object))
}

# Get the n x m counts matrix, where n is the number of samples
# (cells) and m is the number of selected features (columns). An
# additional step is taken to remove all-zero columns.
#
#' @importFrom Matrix colSums
#' @importFrom Matrix t
prepare_counts_fasttopics <- function (object, features) {
  X <- GetAssayData(object,"counts")
  if (is.null(features))
    features <- rownames(X)
  else
    features <- intersect(features,rownames(X))
  X <- X[features,]
  X <- t(X)

  # Remove all-zero columns.
  i <- which(colSums(X > 0) >= 1)
  return(X[,i])
}

# Generate a Seurat PCA dimension reduction object from the loadings
# matrix.
#
#' @importFrom stats prcomp
pca_from_loadings_fasttopics <- function (fit, assay, reduction.key,
                                          min.sdev = 1e-8) {
  k   <- ncol(fit$L)
  out <- prcomp(fit$L)
  colnames(out$x) <- paste0(reduction.key,1:k)
  colnames(out$rotation) <- paste0(reduction.key,1:k)
  cols <- which(out$sdev > min.sdev)
  return(CreateDimReducObject(out$x[,cols],out$rotation[,cols],assay = assay,
                              key = reduction.key,global = TRUE))
}
