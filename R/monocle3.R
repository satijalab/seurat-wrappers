#' @include internal.R
#' @importFrom Seurat DefaultAssay Idents<-
#' @importFrom methods as slot<- slot
#'
NULL

clusters.key <- 'monocle3_clusters'
partitions.key <- 'monocle3_partitions'

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
#' @details The \code{\link[Seurat]{Seurat}} method utilizes
#' \code{\link[Seurat]{as.SingleCellExperiment}} to transfer over expression
#' and cell-level metadata. The following additional information is also
#' transferred over:
#' \itemize{
#'  \item Cell emebeddings are transferred over to the
#'  \code{\link[SingleCellExperiment]{reducedDims}} slot. Dimensional reduction
#'  names are converted to upper-case (eg. \dQuote{umap} to \dQuote{UMAP}) to
#'  match Monocle 3 style
#'  \item Feature loadings are transfered to
#'  \code{cds@reduce_dim_aux$gene_loadings} if present. \strong{NOTE}: only the
#'  feature loadings of the last dimensional reduction are transferred over
#'  \item Standard deviations are added to
#'  \code{cds@reduce_dim_aux$prop_var_expl} if present. \strong{NOTE}: only the
#'  standard deviations of the last dimensional reduction are transferred over
#'  \item Clustering information is transferred over in the following manner: if
#'  cell-level metadata entries \dQuote{monocle3_clusters} and
#'  \dQuote{monocle3_partitions} exist, then these will be set as the clusters
#'  and partitions, with no nearest neighbor graph being added to the object;
#'  otherwise, Seurat's nearest-neighbor graph will be converted to an
#'  \code{\link[igraph]{igraph}} object and added to the \code{cell_data_set}
#'  object along with Seurat's clusters. No partition information is added when
#'  using Seurat's clsuters
#' }
#'
#' @seealso \code{\link[Seurat]{as.SingleCellExperiment}}
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
  # Ensure we have a counts assay
  if (is.null(x = SummarizedExperiment::assays(x = cds)$counts)) {
    SummarizedExperiment::assays(x = cds)$counts <- SummarizedExperiment::assays(x = cds)[[1]]
  }
  # Add Size_factor
  if (!"Size_Factor" %in% colnames(x = SummarizedExperiment::colData(x = cds))) {
    size.factor <- paste0('nCount_', assay)
    if (size.factor %in% colnames(x = x[[]])) {
      SummarizedExperiment::colData(x = cds)$Size_Factor <- x[[size.factor, drop = TRUE]]
    }
  }
  # Add DimReducs: Embeddings become a reduced dim, Loadings go to
  # reduce_dim_aux$gene_loadings, Stdev goes go reduce_dim_aux$prop_var_expl
  # First, reset the ones from as.SingleCellExperiment
  SingleCellExperiment::reducedDims(x = cds)[SingleCellExperiment::reducedDimNames(x = cds)] <- NULL
  reductions <- intersect(
    x = reductions,
    y = AssociatedDimReducs(object = x, assay = assay)
  )
  for (reduc in reductions) {
    SingleCellExperiment::reducedDims(x = cds)[[toupper(x = reduc)]] <- Embeddings(object = x[[reduc]])
    loadings <- Loadings(object = x[[reduc]])
    if (!IsMatrixEmpty(x = loadings)) {
      slot(object = cds, name = 'reduce_dim_aux')[['gene_loadings']] <- loadings
    }
    stdev <- Stdev(object = x[[reduc]])
    if (length(x = stdev)) {
      slot(object = cds, name = 'reduce_dim_aux')[['prop_var_expl']] <- stdev
    }
  }
  # Add clustering information
  # TODO: Figure out if I need to add relations, distMatrix, or clusters/partitions
  if (!is.null(x = group.by)) {
    Idents(object = x) <- group.by
  }
  # if (clusters.key %in% colnames(x = x[[]])) {
  clusters.list <- if (is.null(x = group.by) && all(c(clusters.key, partitions.key) %in% colnames(x = x[[]]))) {
    message("Using existing Monocle 3 cluster membership and partitions")
    list(
      partitions = factor(x = x[[partitions.key, drop = TRUE]]),
      clusters = factor(x = x[[clusters.key, drop = TRUE]])
    )
  } else if (graph %in% names(x = x)) {
    g <- igraph::graph_from_adjacency_matrix(
      adjmatrix = x[[graph]],
      weighted = TRUE
    )
    # TODO: figure out proper partitioning scheme
    # partitions <- igraph::components(graph = g)$membership[colnames(x = x)]
    warning(
      "Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object",
      call. = FALSE,
      immediate. = TRUE
    )
    partitions <- rep_len(x = 1, length.out = ncol(x = x))
    list(
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
      partitions = factor(x = partitions),
      clusters = Idents(object = x)
    )
  } else {
    list()
  }
  if (length(x = clusters.list)) {
    slot(object = cds, name = 'clusters')[[toupper(x = default.reduction)]] <- clusters.list
  }
  # TODO: Add translated results from learn_graph
  return(cds)
}

#' @param loadings Name of dimensional reduction to save loadings to, if present;
#' defaults to first dimensional reduction present (eg.
#' \code{SingleCellExperiment::reducedDimNames(x)[1]}); pass \code{NA} to
#' suppress transfer of loadings
#' @param clusters Name of clustering method to use for setting identity classes
#'
#' @importFrom Seurat as.Seurat Loadings<- as.Graph DefaultAssay<-
#'
#' @details The \code{cell_data_set} method for \code{\link[Seurat]{as.Seurat}}
#' utilizes the \code{\link[Seurat::as.Seurat]{SingleCellExperiment}} method of
#' \code{\link[Seurat]{as.Seurat}} to handle moving over expression data, cell
#' embeddings, and cell-level metadata. The following additional information
#' will also be transfered over:
#' \itemize{
#'  \item Feature loadings from \code{cds@reduce_dim_aux$gene_loadings} will be
#'  added to the dimensional reduction specified by \code{loadings} or the name
#'  of the first dimensional reduction that contains "pca" (case-insensitive) if
#'  \code{loadings} is not set
#'  \item Monocle 3 clustering will be set as the default identity class. In
#'  addition, the Monocle 3 clustering will be added to cell-level metadata as
#'  \dQuote{monocle3_clusters}, if present
#'  \item Monocle 3 partitions will be added to cell-level metadata as
#'  \dQuote{monocle3_partitions}, if present
#'  \item Monocle 3 pseudotime calculations will be added to
#'  \dQuote{monocle3_pseudotime}, if present
#'  \item The nearest-neighbor graph, if present, will be converted to a
#'  \code{\link[Seurat]{Graph}} object, and stored as
#'  \dQuote{\code{assay}_monocle3_graph}
#' }
#'
#' @seealso \code{\link[Seurat]{as.Seurat.SingleCellExperiment}}
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
  # Pull feature loadings
  lds.reduc <- ifelse(
    test = is.null(x = loadings),
    yes = grep(
      pattern = 'pca',
      x = SingleCellExperiment::reducedDimNames(x = x),
      ignore.case = TRUE,
      value = TRUE
    ),
    no = loadings
  )
  if (length(x = lds.reduc) && !is.na(x = lds.reduc)) {
    loadings <- slot(object = x, name = 'reduce_dim_aux')[['gene_loadings']]
    if (!is.null(x = loadings)) {
      Loadings(object = object[[lds.reduc]], projected = FALSE) <- loadings
    }
  }
  # Pull cluster information and pseudotime
  if (length(x = slot(object = x, name = 'clusters'))) {
    clusters <- clusters %||% DefaultDimReduc(object = object)
    object[[clusters.key]] <- Idents(object = object) <- monocle3::clusters(
      x = x,
      reduction_method = clusters
    )
    object[[partitions.key]] <- monocle3::partitions(
      x = x,
      reduction_method = clusters
    )
    graph <- slot(object = x, name = 'clusters')[[clusters]]$cluster_result$g[]
    try(
      expr = {
        graph <- as.Graph(x = graph)
        DefaultAssay(object = graph) <- DefaultAssay(object = object)
        object[[paste0(DefaultAssay(object = graph), '_monocle3_graph')]] <- graph
      },
      silent = TRUE
    )
    try(
      expr = object[['monocle3_pseudotime']] <- monocle3::pseudotime(
        x = cds,
        reduction_method = clusters
      ),
      silent = TRUE
    )
  }
  # TODO: Pull trajectory information
  return(object)
}

#' Run \code{link[monocle3]{learn_graph}} on a \code{\link[Seurat]{Seurat}} object
#'
#' @param object A \code{\link[Seurat]{Seurat}} object
#' @param reduction Name of reduction to use for learning the pseudotime graph
#' @param ... Arguments passed to \code{\link[monocle3]{learn_graph}}
#'
#' @return A \code{\link[monocle3]{cell_data_set}} object with the pseudotime graph
#'
#' @importFrom Seurat Reductions
#'
#' @seealso \code{\link[monocle3]{learn_graph}} \code{\link[monocle3]{cell_data_set}}
#'
# @export
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
  cds <- as.cell_data_set(
    x = object,
    assay = DefaultAssay(object = object[['UMAP']]),
    reductions = 'UMAP',
    default.reduction = 'UMAP'
  )
  cds <- monocle3::learn_graph(cds = cds, ...)
  return(cds)
  # if (reduction != 'UMAP') {
  #   object[['UMAP']] <- NULL
  # }
  # return(object)
}
