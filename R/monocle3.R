#' @include internal.R
#' @importFrom Seurat DefaultAssay Idents<-
#' @importFrom methods as slot<- slot
#'
NULL

#' Convert objects to Monocle3 \code{cell_data_set} objects
#'
#' @param x An object
#' @param ... Arguments passed to other methods
#'
#' @return A \code{cell_data_set} object
#'
#' @name as.cell_data_set
#' @rdname as.cell_data_set
#'
#' @aliases as.CellDataSet
#'
#' @export
#'
as.cell_data_set <- function(x, ...) {
  CheckPackage(package = 'cole-trapnell-lab/monocle3', repository = 'github')
  UseMethod(generic = 'as.cell_data_set', object = x)
}

#' @inheritParams Seurat::as.SingleCellExperiment
#' @param reductions A vector of dimensional reductions add to the
#' \code{cell_data_set} object; defaults to all dimensional reductions
#' calculated from \code{assay} and all \link[Seurat:IsGlobal]{global}
#' dimensional reductions
#' @param default.reduction Name of dimensional reduction to use for clustering
#' name
#' @param graph Name of graph to be used for clustering results
#' @param group.by Name of cell-level metadata column to use as identites; pass
# \code{NULL} to use the active identites
#'
#' @importFrom Seurat as.SingleCellExperiment GetAssayData Loadings
#' Embeddings Stdev Idents
#'
#' @rdname as.cell_data_set
#' @method as.cell_data_set Seurat
#' @export
#'
as.cell_data_set.Seurat <- function(
  x,
  assay = DefaultAssay(object = x),
  reductions = AssociatedDimReducs(object = x, assay = assay),
  default.reduction = DefaultDimReduc(object = x, assay = assay),
  graph = paste0(assay, '_snn'),
  group.by = NULL,
  ...
) {
  # Add assay data
  # Cheat and use as.SingleCellExperiment
  cds <- as(
    object = as.SingleCellExperiment(x = x, assay = assay),
    Class = 'cell_data_set'
  )
  # Add DimReducs: Embeddings become a reduced dim, Loadings go to
  # preprocess_aux$gene_loadings, Stdev goes go preprocess@aux$prop_var_expl
  # First, reset the ones from as.SingleCellExperiment
  SingleCellExperiment::reducedDims(x = cds)[SingleCellExperiment::reducedDimNames(x = cds)] <- NULL
  reductions <- intersect(
    x = reductions,
    y = AssociatedDimReducs(object = x, assay = assay)
  )
  for (reduc in reductions) {
    SingleCellExperiment::reducedDims(x = cds)[[reduc]] <- Embeddings(object = x[[reduc]])
    loadings <- Loadings(object = x[[reduc]])
    if (!IsMatrixEmpty(x = loadings)) {
      slot(object = cds, name = 'preprocess_aux')[['gene_loadings']] <- loadings
    }
    stdev <- Stdev(object = x[[reduc]])
    if (length(x = stdev)) {
      slot(object = cds, name = 'preprocess_aux')[['prop_var_expl']] <- stdev
    }
  }
  # Add clustering information
  # TODO: Figure out if I need to add relations, distMatrix, or clusters/partitions
  if (!is.null(x = group.by)) {
    Idents(object = x) <- group.by
  }
  if (graph %in% names(x = x)) {
    g <- igraph::graph_from_adjacency_matrix(
      adjmatrix = x[[graph]],
      weighted = TRUE
    )
    partitions <- igraph::components(graph = g)$membership[colnames(x = x)]
    slot(object = cds, name = 'clusters')[[default.reduction]] <- list(
      cluster_result = list(
        g = g,
        relations = NULL,
        distMatrix = 'matrix',
        coord = NULL,
        edge_links = NULL,
        optim_res = list(
          membership = as.integer(x = Idents(object = x)),
          modularity = NA_real_
        )
      ),
      partitions = as.factor(x = partitions),
      clusters = Idents(object = x)
    )
  }
  # TODO: Add translated resutls from learn_graph
  return(cds)
}

#' @param loadings Name of dimensional reduction to save loadings to, if present;
#' defaults to first dimensional reduction present (eg.
#' \code{SingleCellExperiment::reducedDimNames(x)[1]})
#' @param clusters Name of clustering method to use for setting identity classes
#'
#' @importFrom Seurat as.Seurat Loadings<- as.Graph DefaultAssay<-
#'
#' @rdname as.Seurat.extras
#' @method as.Seurat cell_data_set
#' @export
#'
as.Seurat.cell_data_set <- function(
  x,
  counts = 'counts',
  data = NULL,
  assay = 'RNA',
  project = 'cell_data_set',
  loadings = NULL,
  clusters = NULL,
  ...
) {
  CheckPackage(package = 'cole-trapnell-lab/monocle3', repository = 'github')
  # Cheat and pull most information using as.SingleCellExperiment
  # cell_data_set objects inherit SingleCellExperiment
  object <- suppressWarnings(expr = as.Seurat(
    x = as(object = x, Class = 'SingleCellExperiment'),
    assay = assay,
    counts = counts,
    data = data,
    project = project
  ))
  # Pull feature loading
  lds.reduc <- loadings %||% SingleCellExperiment::reducedDimNames(x = x)[1]
  if (!is.na(x = lds.reduc)) {
    loadings <- slot(object = x, name = 'preprocess_aux')[['gene_loadings']]
    if (!is.null(x = loadings)) {
      Loadings(object = object[[lds.reduc]], projected = FALSE) <- loadings
    }
  }
  # Pull cluster information
  if (length(x = slot(object = x, name = 'clusters'))) {
    clusters <- clusters %||% DefaultDimReduc(object = object)
    Idents(object = object) <- monocle3::clusters(
      x = x,
      reduction_method = clusters
    )
    graph <- slot(object = x, name = 'clusters')[[clusters]]$cluster_result$g[]
    graph <- as.Graph(x = graph)
    DefaultAssay(object = graph) <- DefaultAssay(object = object)
    object[[paste0(DefaultAssay(object = graph), '_monocle3_graph')]] <- graph
  }
  # TODO: Pull trajectory information
  return(object)
}

#' @importFrom Seurat Reductions
#'
LearnGraph <- function(object, reduction = DefaultDimReduc(object = object), ...) {
  CheckPackage(package = 'cole-trapnell-lab/monocle3', repository = 'github')
  if (reduction != 'UMAP') {
    if ('UMAP' %in% Reductions(object = object)) {
      ''
    }
    reduc <- object[[reduction]]
    suppressWarnings(expr = object[['UMAP']] <- reduc)
  }
  cds <- as.cell_data_set.Seurat(
    x = object,
    assay = DefaultAssay(object = object[['UMAP']]),
    reductions = 'UMAP',
    default.reduction = 'UMAP',
    ...
  )
  cds <- monocle3::learn_graph(cds = cds, ...)
  return(cds)
  # if (reduction != 'UMAP') {
  #   object[['UMAP']] <- NULL
  # }
  # return(object)
}
