#' Plot of density of observations from single cell data in bivariate hexagon
#'   cells.
#'
#' @param sce A \code{\link[Seurat]{Seurat}} object.
#' @param title A string containing the title of the plot.
#' @param xlab A string containing the title of the x axis.
#' @param ylab A string containing the title of the y axis.
#'
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @import Seurat
#' @import ggplot2
#' @importFrom dplyr as_tibble
#' @rdname plot_hexbin_density
#' @export plot_hexbin_density
#'
#' @author Saskia Freytag
#' @seealso \code{\link[schex]{plot_hexbin_density}}
#' @references Delile, Julien, et al.
#'    "Single cell transcriptomics reveals spatial and temporal dynamics of gene
#'     expression in the developing mouse spinal cord." Development 146.12
#'     (2019): dev173807.
#'
#' @examples
#' #' # For Seurat object
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
#' plot_hexbin_density(pbmc_small)
plot_hexbin_density <- function(sce,
                                                    title=NULL,
                                                    xlab=NULL,
                                                    ylab=NULL) {

  if (class(sce) != "Seurat") {
    stop("Please provide a Seurat object.")
  }

  out <- slot(sce, "misc")$hexbin[[2]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  if(is.null(title)) {
    title <- "Density"
  }

  if(is.null(xlab)) {
    xlab <- "x"
  }

  if(is.null(ylab)) {
    ylab <- "y"
  }

  out <- as_tibble(out)

  ggplot(out, aes_string("x", "y", fill="number_of_cells")) +
    geom_hex(stat = "identity") + scale_fill_viridis_c() +
    theme_classic() + theme(legend.position="bottom") + ggtitle(title) +
    labs(x=xlab, y=ylab) + theme(legend.title=element_blank())

}
