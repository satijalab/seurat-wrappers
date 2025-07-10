#' Run a decoupleR method on a Seurat object
#'
#' Generic wrapper to apply any decoupleR method (e.g., "ulm", "viper", "wmean", "mlm")
#' to the expression matrix from a Seurat object, then store scores in the tools slot.
#'
#' @param object    A Seurat object
#' @param network   A prior-knowledge data.frame with columns: source, target, weight
#' @param method    The method to apply (e.g., "ulm", "viper", "wmean", "mlm")
#' @param assay     The assay in the Seurat object to use (default is "RNA")
#' @param layer     The data layer to use (default is "data")
#' @param store_as  The name to store the results under in the Seurat tools slot (optional)
#' @param verbose   Logical, whether to display messages during the process
#' @param ...       Additional arguments passed to the decoupleR function
#'
#' @return The Seurat object with the decoupleR results stored in the tools slot
#'
#' @references Badia-i-Mompel, P. et al. Bioinformatics Advances, 2022. https://doi.org/10.1093/bioadv/vbac016
#'
#' @export
RunDecoupleRSeurat <- function(
  object,
  network,
  method = "ulm",
  assay = "RNA",
  layer = "data",
  store_as = NULL,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Please install the 'decoupleR' package.")
  }

  expr_mat <- GetAssayData(object, assay = assay, slot = layer) |> as.matrix()
  keep <- apply(expr_mat, 1, function(x) all(is.finite(x)))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  if (nrow(expr_mat) == 0) {
    stop("All rows in the expression matrix are NA or Inf.")
  }

  if (!inherits(expr_mat, "dgCMatrix")) {
    expr_mat <- Matrix::Matrix(expr_mat, sparse = TRUE)
  }

  store_as <- store_as %||% paste0("decoupleR_", method)
  network <- dplyr::distinct(network, source, target, .keep_all = TRUE)

  if (verbose) {
    message("→ Running decoupleR method: ", method)
  }

  result <- tryCatch({
    decoupleR::decouple(expr_mat, network = network, method = method, ...)
  }, error = function(e) {
    stop("Error during decoupleR execution: ", e$message)
  })

  Tool(object) <- list2(!!!Tool(object), !!store_as := result)
  object <- LogSeuratCommand(object = object)
  if (verbose) message("→ Stored in tools[['", store_as, "']]")
  return(object)
}