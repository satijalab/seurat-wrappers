#' @importFrom BiocManager install
#' @importFrom remotes install_github
#'
NULL

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

# Check to ensure a package is installed
#
# @param package Name of pacakge to check
# @param repository Repository that package is available on;
# choose from 'bioconductor', 'github', or 'cran'
# @param ... Extra parameters passed to BiocManager::install, remotes::install_github, or install.packages, depending on \code{repository}
#
#' @importFrom utils install.packages
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
