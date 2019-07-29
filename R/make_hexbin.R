#' Bivariate binning of single cell data into hexagon cells.
#'
#' \code{make_hexbin} returns a
#'   \code{\link[Seurat]{Seurat}} object of binned hexagon cells.
#'
#' @param sce A \code{\link[Seurat]{Seurat}} object.
#' @param nbins The number of bins partitioning the range of the first
#'   component of the chosen dimension reduction.
#' @param dimension_reduction A string indicating the reduced dimension
#'   result to calculate hexagon cell representation of.
#'
#' @details This function bins observations with computed reduced dimension
#'   results into hexagon cells. For a \code{\link[Seurat]{Seurat}} object the
#'   results from this function are stored in \code{slot(object, "misc")}.
#'   The list contains two items. The first item stores a vector specifying the
#'   hexagon ID for each observation. The second item stores a matrix with the
#'   x and y positions of the hexagon cells and the number of observations in
#'   each of them.
#'
#' @return A \code{\link[Seurat]{Seurat}} object.
#' @importFrom hexbin hexbin hcell2xy
#' @import Seurat
#' @rdname make_hexbin
#' @export make_hexbin
#'
#' @author Saskia Freytag
#' @seealso \code{\link[schex]{make_hexbin}}
#' @references Delile, Julien, et al.
#'    "Single cell transcriptomics reveals spatial and temporal dynamics of gene
#'     expression in the developing mouse spinal cord." Development 146.12
#'     (2019): dev173807.
#'
#'  @examples
#' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
make_hexbin <- function(sce,
                                            nbins = 80,
                                            dimension_reduction = "UMAP") {

  if (class(sce) != "Seurat") {
    stop("Please provide a Seurat object.")
  }

  if(!is.element(tolower(dimension_reduction), names(slot(sce, "reductions")))) {
    stop("Specify existing dimension reduction.")
  }

  func <- paste0("dr <- slot(slot(sce, 'reductions')$",
                 tolower(dimension_reduction), ", 'cell.embeddings')")

  eval(parse(text = func))

  res <- .make_hexbin_helper(dr, nbins)
  slot(sce, "misc")$hexbin <- res

  return(sce)

}

.make_hexbin_helper <- function(dr, nbins = 80) {

  xbnds <- range(c(dr[, 1]))
  ybnds <- range(c(dr[, 2]))

  drhex <- hexbin(dr[, 1],
                  dr[, 2],
                  nbins,
                  xbnds = xbnds,
                  ybnds = ybnds,
                  IDs = TRUE
  )
  cID <- slot(drhex, "cID")
  drhex <- cbind(as.numeric(hcell2xy(drhex)$x),
                 as.numeric(hcell2xy(drhex)$y),
                 as.numeric(slot(drhex, "count")))

  colnames(drhex) <- c("x", "y", "number_of_cells")

  res <- list(cID=cID, hexbin.matrix=drhex)

  return(res)

}




