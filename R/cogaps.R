#' @include internal.R
#'
NULL

#' Run CoGAPs on a Seurat object
#'
#' @param object Seurat object
#' @param assay Assay to pull data from
#' @param slot Slot to pull data from.
#' @param params \code{\link[CoGAPS]{CogapsParams}} object for specifying parameter settings
#' @param temp.file Name of temporary data matrix file to create if running in a distributed mode.
#' Setting to TRUE will generate the file name using \code{tempfile}.
#' @param reduction.name Name of the CoGAPS reduction returned
#' @param reduction.key Key for the CoGAPS reduction returned
#'
#' @return Returns a Seurat object with the CoGAPS results stored as a \code{\link{DimReduc}} object
#' @seealso \code{\link[CoGAPS]{CoGAPS}}
#' @references E.J. Fertig, J. Ding, A.V. Favorov, G. Parmigiani, and M.F. Ochs (2010) CoGAPS: an
#' integrated R/C++ package to identify overlapping patterns of activation of biological processes
#' from expression data. Bioinformatics 26:2792-2793.
#'
#' @importFrom Matrix writeMM
#'
#' @export
RunCoGAPS <- function(
  object,
  assay = NULL,
  slot = "counts",
  params = NULL,
  temp.file = NULL,
  reduction.name = "CoGAPS",
  reduction.key = "CoGAPS_",
  ...
) {
  SeuratWrappers:::CheckPackage(package = 'CoGAPS', repository = 'bioconductor')
  assay <- assay %||% DefaultAssay(object = object)
  dat <- GetAssayData(object = object, assay = assay, slot = slot)
  dat <- log2(x = dat + 1)
  geneNames <- rownames(x = dat)
  sampleNames <- colnames(x = dat)
  if (!is.null(temp.file)) {
    if (isTRUE(x = temp.file)) {
      temp.file <- paste0(tempfile(), ".mtx")
    } else if (file.exists(temp.file)) {
      stop("temp.file already exists and would be overwritten. Please either remove or specify a new name.")
    }
    dat <- as(object = dat, Class = "dgCMatrix")
    Matrix::writeMM(obj = dat, file = temp.file)
    dat <- temp.file
  } else {
    dat <- as.matrix(x = dat)
  }
  if (!is.null(x = params)) {
    CoGAPS_results <- CoGAPS::CoGAPS(
      data = dat,
      params = params,
      geneNames = geneNames,
      sampleNames = sampleNames,
      ...
    )
  } else {
    CoGAPS_results <- CoGAPS::CoGAPS(
      data = dat,
      geneNames = geneNames,
      sampleNames = sampleNames,
      ...
    )
  }
  object[["CoGAPS"]] <- CreateDimReducObject(
    embeddings = slot(object = CoGAPS_results, name = "sampleFactors"),
    loadings = slot(object = CoGAPS_results, name = "featureLoadings"),
    key = "CoGAPS_",
    assay = assay
  )
  return(object)
}
