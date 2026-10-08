\name{RunLocalMAP}
\alias{RunLocalMAP}
\alias{RunLocalMAP.Seurat}
\alias{RunLocalMAP.default}
\title{Run LocalMAP (Local Manifold Approximation Projection)}
\usage{
RunLocalMAP(object, ...)

\method{RunLocalMAP}{Seurat}(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  slot = "data",
  n_components = 2,
  n_neighbors = NULL,
  lr = 1,
  num_iters = 250L,
  init = "random",
  reduction.name = "localmap",
  reduction.key = "LocalMAP_",
  verbose = TRUE,
  seed.use = 42L,
  ...
)

\method{RunLocalMAP}{default}(
  object,
  assay = NULL,
  n_components = 2,
  n_neighbors = NULL,
  lr = 1,
  num_iters = 250L,
  init = "random",
  reduction.key = "LocalMAP_",
  verbose = TRUE,
  seed.use = 42L,
  ...
)
}
\arguments{
  \item{object}{A Seurat object or matrix-like object (e.g., cells × genes).}
  \item{\dots}{Additional arguments passed to the `localmap.LocalMAP` Python constructor.}
  \item{reduction}{Name of input reduction (e.g., “pca”) when using a Seurat object.}
  \item{dims}{Integer vector specifying dimensions of input reduction.}
  \item{features}{Character vector of feature names to use instead of dims.}
  \item{assay}{Assay name in Seurat (default: \code{NULL}, uses default assay).}
  \item{slot}{Slot in Seurat assay (default: "data").}
  \item{n_components}{Number of output embedding dimensions (default: 2).}
  \item{n_neighbors}{Number of neighbors for graph construction; if \code{NULL}, defaults to 10 (or scaled for large datasets).}
  \item{lr}{Learning rate for optimizer (default: 1).}
  \item{num_iters}{Number of optimization iterations (default: 250).}
  \item{init}{Initialization method: "random" or "pca" (default: "random").}
  \item{reduction.name}{Name under which embedding is stored in Seurat object (default: "localmap").}
  \item{reduction.key}{Prefix for embedding column names (default: "LocalMAP_").}
  \item{verbose}{Logical: print progress messages (default: TRUE).}
  \item{seed.use}{Integer seed for reproducibility (default: 42).}
}
\description{
Performs LocalMAP, an algorithm for dimension Reduction with Locally Adjusted Graphs.
Constructs a local-neighbor graph and optimizes a low-dimensional embedding to preserve neighborhood structure.
Intended to be used on high-dimensional data such as single-cell RNA-seq.
}
\details{
This function wraps the Python \code{pacmap.LocalMAP} algorithm via \pkg{reticulate}. Inputs can be raw feature matrices or precomputed reductions (e.g., PCA). The Seurat method adds support for \code{features} or \code{dims} selection and stores the resulting embedding as a reduction in the Seurat object.
}
\examples{
## Using Seurat object
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
seurat_obj <- RunLocalMAP(
  object = seurat_obj,
  features = Seurat::VariableFeatures(seurat_obj),
  n_components = 2,
  seed.use = 123
)
DimPlot(seurat_obj, reduction = "localmap")

## Using matrix directly
mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
dr <- RunLocalMAP.default(mat, n_components = 3, init = "pca")
head(dr@cell.embeddings)
}
\author{
Yiyang Sun
}
