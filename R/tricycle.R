#' @include internal.R
#'
NULL

#' Run estimate_cycle_position on a Seurat object
#'
#' This function run estimate_cycle_position function on Seurat object. It uses
#' the tricycle internal reference projection matrix.
#'
#' @param object Seurat object
#' @param assay Assay to use, defaults to the default assay
#' @param slot Slot to use. It should be library size adjusted **log-expression** values.
#' Note that it is convention that we rename "logcounts" to "data" when converting SingleCellExperiment to Seurat object.
#' See also \code{\link[Seurat]{as.Seurat}}. Defaults to "data"
#' @param reduction.name Name of the cell cycle projection returned
#' @param reduction.key Key for the cell cycle projection returned
#' @param gname Alternative rownames of \code{object}. If provided, this will be used to map genes within \code{object} with genes in reference.
#' If not provided, the rownames of \code{object} will be used instead. Default: NULL
#' @param gname.type The type of gene names as in \code{gname} or rownames of \code{object}. It can be either 'ENSEMBL' or 'SYMBOL'. Default: 'ENSEMBL'
#' @param species The type of species in \code{object}. It can be either 'mouse' or 'human'. Default: 'mouse'
#' @param AnnotationDb An AnnotationDb objects. If the user provides rownames in the format of Ensembl IDs and project human data,
#'  this object will be used to map Ensembl IDs to gene SYMBOLs. If no AnnotationDb object being given, the function will use \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}.
#' @param center.pc1 The center of PC1 when defining the angle. Default: 0
#' @param center.pc2 The center of PC2 when defining the angle. Default: 0
#'
#' @export
Runtricycle <- function(
  object,
  assay = NULL,
  slot = "data",
  reduction.name = "tricycleEmbedding",
  reduction.key = "tricycleEmbedding_",
  gname = NULL,
  gname.type = c("ENSEMBL", "SYMBOL"),
  species = c("mouse", "human"),
  AnnotationDb = NULL,
  center.pc1 = 0,
  center.pc2 = 0) {
  SeuratWrappers:::CheckPackage(package = 'tricycle', repository = 'bioconductor')
  assay <- assay %||% DefaultAssay(object = object)
  data.m <- GetAssayData(object = object, assay = assay, slot = slot)

  projection.m <- tricycle:::.project_cycle_space(data.m, ref.m = NULL, gname = gname, gname.type = gname.type, species = species, AnnotationDb = AnnotationDb)
  object$tricyclePosition <- tricycle:::.getTheta(projection.m, center.pc1 = center.pc1, center.pc2 = center.pc2)

  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = projection.m,
    key = reduction.key,
    assay = assay
  )
  return(object)

}
