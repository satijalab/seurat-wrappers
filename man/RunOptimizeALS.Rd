% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liger.R
\name{RunOptimizeALS}
\alias{RunOptimizeALS}
\alias{optimizeALS}
\title{Run optimizeALS on a Seurat object}
\usage{
RunOptimizeALS(
  object,
  k,
  assay = NULL,
  split.by = "orig.ident",
  lambda = 5,
  thresh = 1e-06,
  max.iters = 30,
  reduction.name = "iNMF_raw",
  reduction.key = "riNMF_",
  nrep = 1,
  H.init = NULL,
  W.init = NULL,
  V.init = NULL,
  rand.seed = 1,
  print.obj = FALSE,
  ...
)
}
\arguments{
\item{object}{A merged Seurat object}

\item{k}{Inner dimension of factorization (number of factors). Run suggestK to determine
appropriate value; a general rule of thumb is that a higher k will be needed for datasets with
more sub-structure.}

\item{assay}{Assay to use, defaults to the default assay of the first object}

\item{split.by}{Attribute for splitting, defaults to "orig.ident"}

\item{lambda}{Regularization parameter. Larger values penalize dataset-specific effects more
strongly (ie. alignment should increase as lambda increases). Run suggestLambda to determine
most appropriate value for balancing dataset alignment and agreement (default 5.0).}

\item{thresh}{Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh.
(default 1e-6)}

\item{max.iters}{Maximum number of block coordinate descent iterations to perform (default 30).}

\item{reduction.name}{Name to store resulting DimReduc object as}

\item{reduction.key}{Key for resulting DimReduc}

\item{nrep}{Number of restarts to perform (iNMF objective function is non-convex, so taking the
best objective from multiple successive initializations is recommended). For easier
reproducibility, this increments the random seed by 1 for each consecutive restart, so future
factorizations of the same dataset can be run with one rep if necessary. (default 1)}

\item{H.init}{Initial values to use for H matrices. (default NULL)}

\item{W.init}{Initial values to use for W matrix (default NULL)}

\item{V.init}{Initial values to use for V matrices (default NULL)}

\item{rand.seed}{Random seed to allow reproducible results (default 1).}

\item{print.obj}{Print objective function values after convergence (default FALSE).}

\item{...}{Arguments passed to other methods}
}
\value{
A Seurat object with embeddings and loadings from \code{\link[liger]{optimizeALS}}
stored as a DimReduc object with name \code{reduction.name} (key set to \code{reduction.key});
per-dataset feature loadings matrices stored in the \code{tool} slot, accessible with
\code{\link[Seurat]{Tool}}
}
\description{
Run optimizeALS on a Seurat object
}
\seealso{
\code{\link[rliger]{optimizeALS}} \code{\link[Seurat]{Tool}}
}
