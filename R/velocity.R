#' @include internal.R
#'
NULL

#' @inheritParams Seurat::CreateSeuratObject
#' @param default.assay Name or index of matrix to use as default assay;
#' defaults to name of first matrix in list
#' @param slot Name of slot to store matrix in; choose from 'counts' or 'data'
#'
#' @importFrom methods new
#' @importFrom utils txtProgressBar packageVersion setTxtProgressBar
#' @importFrom Seurat as.Seurat CreateAssayObject Key<- CreateSeuratObject
#'
#' @details
#' The \code{list} method for \code{\link[Seurat]{as.Seurat}} takes a named list
#' of matrices (dense or sparse) and creates a single \code{Seurat} object where
#' each matrix is its own assay. The names of the list are taken to be the names
#' of the assays. If not present, assays will be named as "Assay#" where "#" is
#' the index number in the list of matrices. Objects will be constructed as follows:
#' \itemize{
#'   \item By default, all matrices are assumed to be raw counts and will be stored
#'   in the \code{counts} slot. This can be changed to store in the matrix in the
#'   \code{data} slot instead. The \code{slot} parameter is vectorized, so different
#'   matrices can be stored in either \code{counts} or \code{data}
#'   \item For any and all matrices designated as \code{counts}, the \code{min.cells}
#'   and \code{min.features} filtering will be applied. These parameters are vectorized,
#'   so different filterings can be applied to different matrices
#'   \item No extra information (eg. \code{project}) can be provided to
#'   \code{\link[Seurat]{CreateSeuratObject}}
#' }
#'
#' @rdname as.Seurat.extras
#' @export
#' @method as.Seurat list
#'
as.Seurat.list <- function(
  x,
  default.assay = 1,
  slot = 'counts',
  min.cells = 0,
  min.features = 0,
  verbose = TRUE,
  ...
) {
  if (!all(sapply(X = x, FUN = inherits, what = c('matrix', 'dgCMatrix')))) {
    stop("All values must be either a matrix or dgCMatrix", call. = FALSE)
  }
  names(x = x) <- names(x = x) %||% rep_len(x = '', length.out = length(x = x))
  names(x = x)[nchar(x = names(x = x)) == 0] <- paste0('Assay', which(x = nchar(x = names(x = x)) == 0))
  if (is.numeric(x = default.assay)) {
    default.assay <- names(x = x)[default.assay]
  }
  if (!default.assay %in% names(x = x)) {
    stop(
      "Cannot find specified default assay '",
      default.assay,
      "' in the list of matrices",
      call. = FALSE
    )
  }
  slot <- rep_len(x = slot, length.out = length(x = x))
  min.cells <- rep_len(x = min.cells, length.out = length(x = x))
  min.features <- rep_len(x = min.features, length.out = length(x = x))
  if (!all(slot %in% c('counts', 'data'))) {
    stop("'slot' must be either 'counts' or 'data'")
  }
  names(x = slot) <- names(x = min.cells) <- names(x = min.features) <- names(x = x)
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(x = x), style = 3, file = stderr())
  }
  if (slot[[default.assay]] == 'data') {
    assays <- list(CreateAssayObject(data = x[[default.assay]]))
    names(x = assays) <- default.assay
    suppressWarnings(expr = Key(object = assays[[default.assay]]) <- tolower(x = default.assay))
    object <- new(
      Class = 'Seurat',
      assays = assays,
      meta.data = data.frame(row.names = colnames(x = assays[[default.assay]])),
      version = packageVersion(pkg = 'Seurat'),
      project.name = 'SeuratProject'
    )
    DefaultAssay(object = object) <- default.assay
  } else {
    object <- CreateSeuratObject(
      counts = x[[default.assay]],
      assay = default.assay,
      min.cells = min.cells[[default.assay]],
      min.features = min.features[[default.assay]]
    )
  }
  if (verbose) {
    setTxtProgressBar(pb = pb, value = 1 + pb$getVal())
  }
  for (i in names(x = x)) {
    if (i == default.assay) {
      next
    }
    if (slot[[i]] == 'data') {
      suppressWarnings(expr = object[[i]] <- CreateAssayObject(data = x[[i]]))
    } else {
      suppressWarnings(
        expr = object[[i]] <- CreateAssayObject(
          counts = x[[i]],
          min.cells = min.cells[[i]],
          min.features = min.features[[i]]
        )
      )
    }
    if (verbose) {
      setTxtProgressBar(pb = pb, value = 1 + pb$getVal())
    }
  }
  if (verbose) {
    close(con = pb)
  }
  return(object)
}

#' Load RNA Velocity data from a loom file
#'
#' This is a wrapper around \code{\link[velocyto.R]{read.loom.matrices}}, but sends
#' messages to \code{stderr} instead of \code{stdout} (or silences messages with
#' \code{verbose = FALSE})
#'
#' @param file Path to loom file
#' @param engine Method to load data data, choose from 'hdf5r' or 'h5'
#' @param verbose Display progress updates
#'
#' @importFrom utils capture.output
#'
#' @export
#'
#' @seealso \code{\link[velocyto.R]{read.loom.matrices}}
#'
ReadVelocity <- function(file, engine = 'hdf5r', verbose = TRUE) {
  CheckPackage(package = 'velocyto-team/velocyto.R', repository = 'github')
  if (verbose) {
    sink(file = stderr(), type = 'output')
    on.exit(expr = sink())
    ldat <- velocyto.R::read.loom.matrices(file = file, engine = engine)
  } else {
    invisible(x = capture.output(ldat <- velocyto.R::read.loom.matrices(
      file = file,
      engine = engine
    )))
  }
  return(ldat)
}

#' Run RNA Velocty
#'
#' @param object A \code{Seurat} object
#' @param spliced Name of spliced assay
#' @param unspliced Name of unspliced assay
#' @param ambiguous Optional name of ambiguous assay
#' @param spliced.average,unspliced.average Required minimum average expression count for the spliced and unspliced expression matrices
#' @param reduction Name of reduction to use
#' @param group.by Factor to group cells by
#' @param cells Vector of cells to use; defaults to all cells
#' (see \code{\link[velocyto.R]{gene.relative.velocity.estimates}:steady.state.cells})
#' @param graph Optional name of nearest neighbor graph to use
#' @param ncores Number of cores to use
#' @param verbose Display progress updates
#' @param ... Extra parameters passed to \code{\link[velocyto.R]{gene.relative.velocity.estimates}}
#'
#' @return ...
#'
#' @importFrom stats as.dist
#' @importFrom Seurat FetchData GetAssayData
#'
#' @export
#'
#' @seealso \code{\link[velocyto.R]{gene.relative.velocity.estimates}} \code{\link[Seurat]{Tool}}
#'
RunVelocity <- function(
  object,
  spliced = 'spliced',
  unspliced = 'unspliced',
  ambiguous = NULL,
  spliced.average = 0.2,
  unspliced.average = 0.05,
  reduction = 'pca',
  group.by = 'ident',
  cells = NULL,
  graph = NULL,
  ncores = 1,
  verbose = TRUE,
  ...
) {
  CheckPackage(package = 'velocyto-team/velocyto.R', repository = 'github')
  # return(invisible(x = NULL))
  # Collect data from Seurat object
  clusters <- FetchData(object = object, vars = group.by)[, , drop = TRUE]
  names(x = clusters) <- colnames(x = object)
  if (!is.factor(x = clusters)) {
    clusters <- as.factor(x = clusters)
  }
  if (verbose) {
    message("Filtering genes in the spliced matrix")
  }
  spliced.matrix <- velocyto.R::filter.genes.by.cluster.expression(
    emat = GetAssayData(object = object, assay = spliced),
    clusters = clusters,
    min.max.cluster.average = spliced.average
  )
  if (verbose) {
    message("Filtering genes in the unspliced matrix")
  }
  unspliced.matrix <- velocyto.R::filter.genes.by.cluster.expression(
    emat = GetAssayData(object = object, assay = unspliced),
    clusters = clusters,
    min.max.cluster.average = unspliced.average
  )
  if (verbose) {
    message("Calculating embedding distance matrix")
  }
  cell.dist <- as.dist(
    m = 1 - velocyto.R::armaCor(
      mat = t(x = Embeddings(object = object, reduction = reduction))
    )
  )
  # Set arguments
  args <- list(...)
  defaults <- as.list(x = formals(fun = velocyto.R::gene.relative.velocity.estimates))
  args <- args[intersect(x = names(x = args), y = names(x = defaults))]
  defaults.use <- setdiff(x = names(x = defaults), y = names(x = args))
  args[defaults.use] <- defaults[defaults.use]
  args$emat <- spliced.matrix
  args$nmat <- unspliced.matrix
  args$smat <- ambiguous %iff% GetAssayData(object = object, assay = ambiguous)
  args$steady.state.cells <- cells %||% colnames(x = object)
  args$cell.dist <- cell.dist
  args$cellKNN <- graph %iff% object[[graph]]
  args$n.cores <- ncores
  args$verbose <- verbose
  # Run velocity
  sink(file = stderr(), type = 'output')
  on.exit(expr = sink())
  cd <- do.call(what = velocyto.R::gene.relative.velocity.estimates, args = args)
  Tool(object = object) <- cd
  return(object)
}

#' RNA Velocity Plot
#'
#' @inheritParams Seurat::DimPlot
#' @param ... Extra parameters passed on to \code{\link[velocyto.R]{show.velocity.on.embedding.cor}}
#'
#' @return Nothing, shows plot
#'
#' @importFrom Seurat Tool Embeddings
#'
# @export
#'
#' @keywords internal
#'
#' @seealso \code{\link[velocyto.R]{show.velocity.on.embedding.cor}}
#'
VeloPlot <- function(
  object,
  reduction = NULL,
  ...
) {
  .NotYetImplemented()
  CheckPackage(package = 'velocyto-team/velocyto.R', repository = 'github')
  velocity <- Tool(object = object, slot = 'RunVelocity')
  if (is.null(x = velocity)) {
    stop("Please run RunVelocity on this Seurat object")
  }
  reduction <- reduction %||% {
    default.reductions <- c("umap", "tsne", "pca")
    object.reductions <- Filter(
      f = function(x) {
        return(inherits(x = object[[x]], what = 'DimReduc'))
      },
      x = names(x = object)
    )
    reduc.use <- min(which(x = default.reductions %in% object.reductions))
    default.reductions[reduc.use]
  }
  embeddings <- Embeddings(object = object, reduction = reduction)
  velocyto.R::show.velocity.on.embedding.cor(
    emb = embeddings,
    vel = velocity,
    ...
  )
  return(invisible(x = NULL))
}
