#' Run a decoupleR method on a Seurat object
#'
#' Generic wrapper to apply any decoupleR method (e.g., "ulm", "viper", "wmean", "mlm")
#' to the expression matrix from a Seurat object, then store scores in the tools slot.
#'
#' @param object      A Seurat object
#' @param network     A prior-knowledge data.frame with columns: source, target, weight
#' @param method      The method to apply (e.g., "ulm", "viper", "wmean", "mlm")
#' @param assay       The assay in the Seurat object to use (default is "RNA")
#' @param layer       The data layer to extract (default is "data")
#' @param store_as    The name under which results are stored in the tools slot (optional)
#' @param verbose     Logical, whether to display progress messages
#' @param ...         Additional arguments passed to the decoupleR function
#'
#' @return The Seurat object with decoupleR results stored in the tools slot
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
  # Check if decoupleR is installed
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Please install the 'decoupleR' package.")
  }

  # Extract the expression matrix from the specified assay and layer
  expression_matrix <- GetAssayData(object, assay = assay, slot = layer) |> as.matrix()

  # Remove rows (genes) that contain only non-finite values
  valid_rows <- apply(expression_matrix, 1, function(gene_expr) all(is.finite(gene_expr)))
  expression_matrix <- expression_matrix[valid_rows, , drop = FALSE]

  if (nrow(expression_matrix) == 0) {
    stop("All rows in the expression matrix are non-finite (NA or Inf).")
  }

  # Ensure matrix is in sparse format (dgCMatrix)
  if (!inherits(expression_matrix, "dgCMatrix")) {
    expression_matrix <- Matrix::Matrix(expression_matrix, sparse = TRUE)
  }

  # Generate a default name for storing results if not provided
  result_name <- store_as %||% paste0("decoupleR_", method)

  # Deduplicate network interactions (based on source/target)
  network <- dplyr::distinct(network, source, target, .keep_all = TRUE)

  if (verbose) {
    message("→ Running decoupleR method: ", method)
  }

  # Run the decoupleR method on the expression matrix
  result <- tryCatch({
    decoupleR::decouple(expression_matrix, network = network, method = method, ...)
  }, error = function(e) {
    stop("Error during decoupleR execution: ", e$message)
  })

  # Store the results in the Seurat object tools slot
  Tool(object) <- list2(!!!Tool(object), !!result_name := result)

  # Log the command in the Seurat object
  object <- LogSeuratCommand(object = object)

  if (verbose) {
    message("→ Results stored in object@tools[['", result_name, "']]")
  }

  return(object)
}
