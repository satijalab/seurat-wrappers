#' @include internal.R
#'
NULL

#' Coralysis integration in Seurat v5 
#' 
#' @description Run Coralysis integration in Seurat v5.
#' 
#' @param object A Seurat object with joint layers. 
#' @param batch Cell metadata column name representing the batch label identity 
#' to guide integration. By default uses the default assay. 
#' @param threads Threads to pass on to \code{Coralysis}. By default it is not 
#' performed parallelization. 
#' @param ndims A positive integer denoting the number of principal components 
#' to calculate and select. Default is \code{50}.
#' @param assay Seurat assay to retrieve layers from. By default uses the default 
#' assay. 
#' @param layers Seurat layers to retrieve the normalized expression matrix from. 
#' By default uses the layer \code{data}, i.e., the normalized expression data. 
#' @param features A vector of features to use for integration. Otherwise 2,000
#' highly variable features are internally selected by default. Set the argument 
#' of \code{scale.layer} to \code{NULL} in the function \code{IntegrateLayers}  
#' to use a set of \code{features} other than default. 
#' @param new.reduction Name of the new integrated dimensional reduction.
#' @param reduction.key Reduction key for the new (integrated) reduction. 
#' @param reconstructed.assay Name of the integrated assay comprising the cell 
#' cluster probabilities from Coralysis. 
#' @param save.coralysis.metadata Filename to save Coralysis metadata locally 
#' as an RDS object. Useful to re-use the ICP models generated during integration. 
#' By default it is not saved. 
#' @param verbose Logical value indicating whether to print progress messages.
#' @param verbose.icp Logical value indicating whether to print progress messages 
#' during the ICP modeling step. Defaults to \code{FALSE}, meaning no messages 
#' are printed. Useful to debug ICP. 
#' @param ... Arguments passed on to \code{Coralysis}. Check possible arguments 
#' with \code{?Coralysis::RunParallelDivisiveICP}.
#' 
#' @importFrom Seurat Assays as.SingleCellExperiment DefaultAssay 
#' SelectIntegrationFeatures5 CreateDimReducObject
#' @importFrom SeuratObject Layers CreateSeuratObject
#' @importFrom SingleCellExperiment cbind reducedDim reducedDim<-
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom rlang check_installed
#' 
#' @return A list with the new integrated reduction (\code{new.reduction}) and 
#' the integrated assay comprising the cell cluster probabilities from Coralysis
#' (\code{reconstructed.assay}).  
#' 
#' @export 
#' 
#' @note This function requires the \href{https://github.com/elolab/Coralysis}{\pkg{Coralysis}} 
#' package to be installed.
#' 
#' @examples
#' \dontrun{
#' # Basic Seurat workflow
#' seu <- SeuratData::LoadData("ifnb")
#' seu[["RNA"]] <- split(seu[["RNA"]], f = seu$stim)
#' seu <- FindVariableFeatures(seu)
#' seu <- ScaleData(seu)
#' seu <- RunPCA(seu)
#' 
#' # Coralysis integration with default parameters 
#' seu <- IntegrateLayers(object = seu, 
#'   method = CoralysisIntegration,
#'   new.reduction = "integrated.coralysis",
#'   batch = "stim",
#'   threads = 4)
#'   
#'  # Coralysis integration with other than default parameters 
#'  # (see other options `?Coralysis::RunParallelDivisiveICP`)
#' seu <- IntegrateLayers(object = seu, 
#'   method = CoralysisIntegration,
#'   new.reduction = "integrated.coralysis",
#'   batch = "stim",
#'   threads = 4, 
#'   k = 8, 
#'   L = 25)
#' }
#' 
#' @author António GG Sousa
#' 
#' @references Sousa A, Smolander J, Junttila S, Elo L (2025). “Coralysis enables 
#' sensitive identification of imbalanced cell types and states in single-cell data 
#' via multi-level integration.” bioRxiv. doi:10.1101/2025.02.07.637023
#' 
#' @seealso \code{\link[Coralysis]{RunParallelDivisiveICP}} \code{\link[Coralysis]{RunPCA}}
#' 
CoralysisIntegration <- function(object, 
                                 batch,
                                 threads = 1,
                                 ndims = 50,
                                 assay = NULL, 
                                 layers = NULL, 
                                 features = NULL,
                                 new.reduction = "integrated.dr",
                                 reduction.key = "coralysis_",
                                 reconstructed.assay = "coralysis",
                                 save.coralysis.metadata = NULL,
                                 verbose = TRUE,
                                 verbose.icp = FALSE,
                                 ...) {
  # Check Coralysis installation
  check_installed(
    pkg = "Coralysis",
    reason = "for running integration with Coralysis"
  )
  
  # Check if input is a Seurat object; otherwise assay
  if (inherits(x = object, what = "Seurat")) { # If Seurat object, check input params
    if (!is.null(assay) && !(assay %in% Assays(object))) {
      stop("'assay' is not available in 'object'")
    } else if (!is.null(layers) && !(layers %in% Layers(object))) {
      stop("'layers' are not available in 'object'") 
    }
    # Define layers & assay
    assay <- assay %||% DefaultAssay(object)
    layers <- layers %||% Layers(object, search = "data")
  } else { # Create Seurat object from Assay5 & add batch label identity
    object <- CreateSeuratObject(object, assay = assay)
    cells <- unlist(lapply(X = layers, FUN = function(l) {
      ncol(LayerData(object = object, layer = l)) 
    }))
    object[[batch]] <- rep(x = layers, times = cells)   
  }
  
  # Features for integration: provide features or 2,000 highly variable are selected
  if (is.null(features)) {
    if (verbose) {
      message("\nSelecting 2000 highly variable integration features.")
    }
    features <- features %||% SelectIntegrationFeatures5(object = object)
  }
  
  # Create SingleCellExperiment object
  if (verbose) {
    cells <- unlist(lapply(X = layers, FUN = function(l) {
      ncol(LayerData(object = object, layer = l)) 
    }))
    message(paste("\nCreating SingleCellExperiment object:\n", 
                  "Assay:", assay, "\n", 
                  "Layers:", paste(layers, collapse = ", "), "\n",
                  "Features:", length(features), 
                  paste0("(", paste(head(features), collapse=", "), " ... ", paste(tail(features), collapse=", "), ")"), "\n", 
                  "Cells:", paste0(paste(cells, collapse = ", "), " (total: ", sum(cells), ")")))
  }
  object.sce.list <- lapply(X = layers, FUN = function(l, f) {
    as.SingleCellExperiment(
      x = subset(x = object, 
                 features = f, 
                 cells = colnames(LayerData(object = object, layer = l)))
    )
  }, f = features)
  object.sce <- do.call(cbind, object.sce.list)
  rm(object)
  invisible(gc())
  
  # Coralysis integration
  if (verbose) {
    message("\nIntegrating data with Coralysis (v.", packageVersion(pkg = "Coralysis"), ").")
  }
  # Parse input params for ICP
  icp.params <- names(as.list(args(Coralysis::RunParallelDivisiveICP)))
  icp.params <- icp.params[icp.params != ""]
  icp.args <- c(list(object = object.sce, batch.label = batch, threads = threads, verbose = verbose.icp), 
                list(...))
  icp.args <- icp.args[icp.params]
  icp.args <- Filter(Negate(is.null), icp.args)
  object.sce <- do.call(what = Coralysis::RunParallelDivisiveICP, args = icp.args)
  
  # Calculating integrated PCA
  if (verbose) {
    message("\nPerforming integrated reduction...")
  }
  object.sce <- Coralysis::RunPCA(object = object.sce, 
                                  assay.name = "joint.probability", 
                                  p = ndims,
                                  scale = TRUE,
                                  center = TRUE,
                                  threshold = 0,
                                  pca.method = "irlba",
                                  return.model = FALSE,
                                  select.icp.tables = NULL,
                                  features = NULL,
                                  dimred.name = new.reduction)
  
  if (verbose) {
    message("\nRetrieving integrated reduction and assay...")
  }
  # Parsing integrated dimred
  colnames(reducedDim(x = object.sce, type = new.reduction)) <- paste0(reduction.key, seq_len(ndims))
  reduction <- CreateDimReducObject(embeddings = reducedDim(x = object.sce, type = new.reduction), 
                                    key = reduction.key, 
                                    assay = assay)
  
  # Creating new assay
  k <- metadata(object.sce)$coralysis$k
  Ks <- log2(k)
  joint.probability <- Coralysis::GetCellClusterProbability(object = object.sce, icp.round = Ks)
  joint.probability <- t(joint.probability)
  colnames(joint.probability) <- colnames(object.sce)
  row.names(joint.probability) <- paste0("k", seq_len(nrow(joint.probability)))
  reconstructed_assay <- CreateAssayObject(data = as(object = joint.probability, Class = "sparseMatrix"))
  
  # Add cells & features used for integration to metadata
  if (!is.null(save.coralysis.metadata)) {
    if (verbose) {
      message("\nSaving Coralysis metadata: ", save.coralysis.metadata)
    }
    metadata(object.sce)$coralysis$cells <- colnames(object.sce)
    metadata(object.sce)$coralysis$features <- row.names(object.sce)
    saveRDS(object = metadata(object.sce), file = save.coralysis.metadata)
  }

  # Return
  res.list <- list(reduction, reconstructed_assay)
  names(res.list) <- c(new.reduction, reconstructed.assay)
  return(res.list)
}

attr(x = CoralysisIntegration, which = "Seurat.method") <- "integration"
