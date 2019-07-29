#' Group label plot position.
#'
#' @param sce A \code{\link[Seurat]{Seurat}} object.
#' @param col The name referring to one column in meta data for which
#'   the label position on the plot is calculated for every level. The chosen
#'   column needs to be a factor.
#'
#' @return A tibble.
#' @importFrom tidyr nest
#' @importFrom cluster pam
#' @import Seurat
#' @import dplyr
#' @rdname make_hexbin_label
#' @export make_hexbin_label
#'
#' @author Saskia Freytag
#' @seealso \code{\link[schex]{make_hexbin_label}}
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
#' make_hexbin_label(pbmc_small, col="RNA_snn_res.1")
make_hexbin_label <- function(sce, col){

  if (class(sce) != "Seurat") {
    stop("Please provide a Seurat object.")
  }

  func <- paste0("!(is.factor(sce$", col, ")|is.character(sce$", col, "))")

  if (any(!col %in% colnames(slot(sce, "meta.data")))) {
    stop("Column cannot be found in slot(sce, 'meta.data').")
  }

  if(eval(parse(text=func))){
    stop("The specified column must be a character or a factor.")
  }

  out <- slot(sce, "misc")$hexbin[[2]]
  cID <- slot(sce, "misc")$hexbin[[1]]

  if(is.null(out)){
    stop("Compute hexbin representation before plotting.")
  }

  action <- "majority"

  name_s <- paste0("sce$", col)
  func <- paste0("x <- ", name_s)

  eval(parse(text = func))

  hh <- .make_hexbin_function(x, action, cID)
  out <- as_tibble(out)

  label <- paste0(col, "_majority")

  func1 <- paste0("out$", label, " <- hh")
  eval(parse(text=func1))

  label.df_2 <- out %>%
    dplyr::select_("x", "y", label) %>%
    dplyr::group_by_(label) %>%
    nest()

  label.df_3 <- lapply(label.df_2$data, function(x) cluster::pam(x, 1)$medoids)
  names(label.df_3) <- unlist(label.df_2 %>% select_(label))

  label.df_3 <- Reduce(rbind, label.df_3)
  label.df_3 <- data.frame(
    x = label.df_3[, 1],
    y = label.df_3[, 2],
    label = unlist(label.df_2 %>% select_(label))
  )
  rownames(label.df_3) <- NULL

  label.df_3
}

