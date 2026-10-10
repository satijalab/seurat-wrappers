#' Find Clean, Cell-Defining Markers via U-method
#'
#' Evaluates expression probability dynamics across single-cell subsets to pinpoint
#' and highlight cell-defining genes. Unlike standard differential expression tools
#' that rely on average log-fold changes (which can be skewed by expression outliers
#' or localized inflation within minor sub-populations), the U-method evaluates
#' maximum-outside cluster contrast. This structural optimization excels at identifying
#' robust markers across highly sparse datasets and advanced spatial transcriptomic platforms,
#' operating with remarkable speed and completely bypassing the need for downstream normalization, 
#' batch correction, or data smoothing.
#' 
#' # Pro tip:
#' In high-resolution datasets, omitting mixed or exceptionally small micro-clusters 
#' from the background contrast yields a significantly cleaner and more reliable analysis.
#'
#' @param object A Seurat object.
#' @param group.by Category to group cells by (e.g., "Ident", "Celltype"). 
#'                 Defaults to the active cell identity if NULL.
#' @param Uscore.threshold Minimum U-score cutoff. Setting a baseline of 0 allows 
#'                         for a broad exploratory evaluation of the entire landscape. 
#'                         Once clusters are well-established, a value of 0.2 serves 
#'                         as a highly robust, empirically tested heuristic. Default is 0.
#' @param n Number of top-ranked markers to export per cluster for the discrete gene set output. Default is 3.
#' @param omit.clusters Character vector of specific clusters (e.g., mixed, transitional, 
#'                         or trace clusters) to exclude from background contrast computations 
#'                         to prevent baseline inflation.
#' @param assay Assay to pull data from (defaults to the Seurat object's DefaultAssay).
#' @param slot Assay slot to use for computing cell expression probabilities (default "counts").
#' @param ... Additional arguments passed directly to the underlying Umethod package execution.
#'
#' @return A named list containing two core components:
#' \itemize{
#'   \item{\code{gene_list_long}}{: A comprehensive data frame ranking every analyzed gene per cluster, detailing calculated U-scores, cluster expression probabilities (P_in), background expression probabilities (P_out), and associated p-values.}
#'   \item{\code{gene_set}}{: A structured matrix or data frame consolidating the top \code{n} highly ranked cell-defining markers structured column-wise by cluster for downstream signature evaluation.}
#' }
#'
#' @references Stein et al. (2026). The U-method: Leveraging expression probability 
#' for robust biological marker detection. bioRxiv preprint. 
#' \url{https://www.biorxiv.org/content/10.64898/2026.03.31.715470v1.full}
#'
#' @export
FindUniqueMarkers.Seurat <- function(
    object,
    group.by = NULL,
    Uscore.threshold = 0,
    n = 3,
    omit.clusters = NULL,
    assay = NULL,
    slot = "counts",
    ...
) {
  # 1. Dependency check
  if (!requireNamespace("Umethod", quietly = TRUE)) {
    stop("Please install 'Umethod' to use this wrapper: devtools::install_github('YanuvS-Dev/Umethod')")
  }
  
  # 2. Extract grouping vector
  assay <- assay %||% Seurat::DefaultAssay(object)
  
  if (is.null(group.by)) {
    group_vector <- Seurat::Idents(object)
  } else {
    group_vector <- object[[group.by, drop = TRUE]]
  }
  
  # 3. CRITICAL CHECK 1: Ensure there are at least 2 distinct clusters
  unique_labels <- unique(na.omit(group_vector))
  if (length(unique_labels) < 2) {
    stop(paste0(
      "The grouping variable contains only ", length(unique_labels), " active cluster(s). ",
      "U-method requires at least 2 distinct clusters to perform background contrast computations."
    ))
  }
  
  # 4. CRITICAL CHECK 2: Ensure all factor levels have > 0 cells
  cell_counts <- table(group_vector)
  zero_cell_clusters <- names(cell_counts)[cell_counts == 0]
  
  if (length(zero_cell_clusters) > 0) {
    missing_labels <- paste(zero_cell_clusters, collapse = ", ")
    
    stop(paste0(
      "label ", missing_labels, " had 0 cells. ",
      "please run again as.factor on the object to ensure clusters have positive cell counts."
    ))
  }
  
  # 5. Invoke core algorithm directly without modifying underlying object metadata
  umethod_results <- Umethod::FindUniqueMarkers(
    obj = object,
    group_by = group.by,
    Uscore = Uscore.threshold,
    n = n,
    omitCluster = omit.clusters,
    ...
  )
  
  return(umethod_results)
}