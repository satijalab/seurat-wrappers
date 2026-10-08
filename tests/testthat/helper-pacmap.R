# Shared helpers for the RunPaCMAP tests.

skip_without_pacmapr <- function() {
  testthat::skip_if_not_installed("pacmapr")
}

skip_without_seurat <- function() {
  testthat::skip_if_not_installed("Seurat")
}

# Small synthetic matrix (rows = cells, cols = features) with an obvious
# two-cluster structure, so a working PaCMAP call has a chance of separating
# the clusters even at low n.
make_toy_matrix <- function(n_per_cluster = 40L, p = 12L, seed = 1L) {
  set.seed(seed)
  n <- 2L * n_per_cluster
  offset <- matrix(0, nrow = n, ncol = p)
  offset[seq_len(n_per_cluster), ] <- 5
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) + offset
  rownames(X) <- paste0("cell", seq_len(n))
  colnames(X) <- paste0("gene", seq_len(p))
  X
}

make_toy_seurat <- function(n_per_cluster = 40L, p = 30L, seed = 1L) {
  testthat::skip_if_not_installed("Seurat")
  X <- make_toy_matrix(n_per_cluster = n_per_cluster, p = p, seed = seed)
  # Seurat wants features x cells
  counts <- t(X) - min(t(X)) + 1
  storage.mode(counts) <- "integer"
  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE, nfeatures = p)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = 5L, verbose = FALSE,
                        features = Seurat::VariableFeatures(obj))
  obj
}
