#' @importFrom BiocManager install
#' @importFrom remotes install_github
#' @importFrom Seurat IsGlobal Reductions
#'
NULL

#' @docType package
#' @name SeuratWrappers-package
#' @rdname SeuratWrappers-package
#'
"_PACKAGE"

# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

# Get dimensional reduction information associated with an assay
#
# @param object A \code{Seurat} object
# @param assay Name of assay that dimensional reduction objects should be
# associated with
# @param global Include global dimensional reductions
#
# @return A vector of dimensional reduction names
#
# @keywords internal
#
AssociatedDimReducs <- function(
  object,
  assay = DefaultAssay(object = object),
  global = TRUE
) {
  return(Filter(
    f = function(x) {
      check <- DefaultAssay(object = object[[x]]) == assay
      if (global) {
        check <- c(check, IsGlobal(object = object[[x]]))
      }
      return(any(check))
    },
    x = Reductions(object = object)
  ))
}

# Find the default DimReduc
#
# Searches for DimReducs matching 'umap', 'tsne', or 'pca', case-insensitive, and
# in that order. Priority given to DimReducs matching the DefaultAssay or assay specified
# (eg. 'pca' for the default assay weights higher than 'umap' for a non-default assay)
#
# @param object A Seurat object
# @param assay Name of assay to use; defaults to the default assay of the object
#
# @return The default DimReduc, if possible
#
#
DefaultDimReduc <- function(object, assay = NULL) {
  assay <- assay %||% DefaultAssay(object = object)
  drs.use <- c('umap', 'tsne', 'pca')
  dim.reducs <- Reductions(object = object)
  drs.assay <- Filter(
    f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    },
    x = dim.reducs
  )
  if (length(x = drs.assay) > 0) {
    index <- lapply(
      X = drs.use,
      FUN = grep,
      x = drs.assay,
      ignore.case = TRUE
    )
    index <- Filter(f = length, x = index)
    if (length(x = index) > 0) {
      return(drs.assay[min(index[[1]])])
    }
  }
  index <- lapply(
    X = drs.use,
    FUN = grep,
    x = dim.reducs,
    ignore.case = TRUE
  )
  index <- Filter(f = length, x = index)
  if (length(x = index) < 1) {
    stop(
      "Unable to find a DimReduc matching one of '",
      paste(drs.use[1:(length(x = drs.use) - 1)], collapse = "', '"),
      "', or '",
      drs.use[length(x = drs.use)],
      "', please specify a dimensional reduction to use",
      call. = FALSE
    )
  }
  return(dim.reducs[min(index[[1]])])
}

# Check to ensure a package is installed
#
# @param package Name of pacakge to check
# @param repository Repository that package is available on;
# choose from 'bioconductor', 'github', or 'cran'
# @param ... Extra parameters passed to BiocManager::install, remotes::install_github, or install.packages, depending on \code{repository}
#
#' @importFrom utils menu install.packages
#
CheckPackage <- function(package, repository, ...) {
  if (!requireNamespace(package = basename(path = package), quietly = TRUE)) {
    if (interactive()) {
      message("Package ", package, " is not yet installed")
      message("Install now?")
      choice <- menu(choices = c('yes', 'no'))
      if (choice == 1) {
        repository <- match.arg(
          arg = tolower(x = repository),
          choices = c('github', 'bioconductor', 'cran')
        )
        switch(
          EXPR = repository,
          'github' = remotes::install_github(repo = package, ...),
          'bioconductor' = BiocManager::install(pkgs = package, ...),
          'cran' = install.packages(pkgs = package, ...),
          stop("Unknown repository ", repository, call. = FALSE)
        )
        return(invisible(x = NULL))
      }
    }
    stop("Unable to find package ", package, ", please install", call. = FALSE)
  }
}

# Check if a matrix is empty
#
# Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#
# @param x A matrix
#
# @return Whether or not \code{x} is empty
#
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}
