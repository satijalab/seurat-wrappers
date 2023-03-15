#' @include internal.R
#'
NULL

# Runs Wilcoxon Rank Sum using the Presto package
#
# @param data.use Data matrix to test
# @param cells.1 Group 1 cells
# @param cells.2 Group 2 cells
# @param verbose Print a progress bar
# @param ... Extra parameters passed to wilcox.test
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# features
#
#' @importFrom stats wilcox.test
#
PrestoDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  # NOTE: do not use logfc from presto
  group.info <- factor(
    c(rep(x = "Group1", length = length(x = cells.1)),
      rep(x = "Group2", length = length(x = cells.2))),
    levels = c("Group1", "Group2"))
  names(x = group.info) <- c(cells.1, cells.2)
  data.use <- data.use[, names(x = group.info), drop = FALSE]
  res <- presto::wilcoxauc(X = data.use, y = group.info)
  res <- res[1:(nrow(x = res)/2), c('pval','auc')]
  colnames(x = res)[1] <- 'p_val'
  return(as.data.frame(x = res, row.names = rownames(x = data.use)))
}

#' A Presto-based implementation of FindMarkers that runs Wilcoxon tests for the given identity classes
#'
#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param reduction Reduction to use in differential expression testing - will test for DE on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping. Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#' @param mean.fxn Function to use for fold change or average difference calculation.
#' If NULL, the appropriate function will be chose according to the slot used
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame. If NULL, the fold change column will be named
#' according to the logarithm base (eg, "avg_log2FC"), or if using the scale.data
#' slot "avg_diff".
#' @param base The base with respect to which logarithms are computed.
#'
#' @importFrom rlang duplicate
#' @importFrom utils assignInNamespace
#' @importFrom Seurat FindMarkers
#'
#' @export
#' @seealso https://github.com/immunogenomics/presto
RunPresto <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = 'data',
  reduction = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  ...
) {
  if (test.use != 'wilcox') {
    stop("Differential expression test must be `wilcox`")
  }

  CheckPackage(package = 'immunogenomics/presto', repository = 'github')
  orig.fxn <- rlang::duplicate(x = Seurat:::WilcoxDETest)
  assignInNamespace(
    x = "WilcoxDETest",
    value = PrestoDETest,
    ns = "Seurat")

  tryCatch(
    expr = res <- FindMarkers(
      object = object,
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = group.by,
      subset.ident = subset.ident,
      assay = assay,
      slot = slot,
      reduction = reduction,
      features = features,
      logfc.threshold = logfc.threshold,
      test.use = "wilcox",
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      verbose = verbose,
      only.pos = only.pos,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed,
      latent.vars = latent.vars,
      min.cells.feature = min.cells.feature,
      min.cells.group = min.cells.group,
      mean.fxn = mean.fxn,
      fc.name = fc.name,
      base = base,
      ...
    ),
    finally = assignInNamespace(
      x = "WilcoxDETest",
      value = orig.fxn,
      ns = "Seurat")
  )

  return(res)
}

#' A Presto-based implementation of FindAllMarkers that runs Wilcoxon tests for all identity classes
#'
#' Finds markers (Wilcoxon-differentially expressed genes) for each of the identity classes in a dataset
#'
#' @inheritParams RunPresto
#' @param node A node to find markers for and all its children; requires
#' \code{\link{BuildClusterTree}} to have been run previously; replaces \code{FindAllMarkersNode}
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, logFC, etc.)
#'
#' @importFrom stats setNames
#' @importFrom rlang duplicate
#' @importFrom utils assignInNamespace
#' @importFrom Seurat FindAllMarkers
#'
#' @export
#'
#' @aliases RunPrestoAllNode
#' @seealso https://github.com/immunogenomics/presto
RunPrestoAll <- function(
  object,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = 'wilcox',
  slot = 'data',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  return.thresh = 1e-2,
  ...
) {
  if (test.use != 'wilcox') {
    stop("Differential expression test must be `wilcox`")
  }

  CheckPackage(package = 'immunogenomics/presto', repository = 'github')
  orig.fxn <- rlang::duplicate(x = Seurat:::WilcoxDETest)
  assignInNamespace(
    x = "WilcoxDETest",
    value = PrestoDETest,
    ns = "Seurat")

  tryCatch(
    expr = res <- FindAllMarkers(
      object = object,
      assay = assay,
      features = features,
      logfc.threshold = logfc.threshold,
      test.use = "wilcox",
      slot = slot,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      node = node,
      verbose = verbose,
      only.pos = only.pos,
      max.cells.per.ident = max.cells.per.ident,
      random.seed = random.seed,
      latent.vars = latent.vars,
      min.cells.feature = min.cells.feature,
      min.cells.group = min.cells.group,
      mean.fxn = mean.fxn,
      fc.name = fc.name,
      base = base,
      return.thresh = return.thresh,
      ...
    ),
    finally = assignInNamespace(
      x = "WilcoxDETest",
      value = orig.fxn,
      ns = "Seurat")
  )

  return(res)
}
