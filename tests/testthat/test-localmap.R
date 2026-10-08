# Tests for RunLocalMAP (R/localmap.R), which now dispatches to
# pacmapr::localmap() instead of the Python pacmap module.

test_that("RunLocalMAP.default returns a DimReduc with expected shape", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  red <- RunLocalMAP(
    object       = X,
    n_components = 2L,
    num_iters    = 50L,
    verbose      = FALSE,
    seed.use     = 11L
  )

  expect_s4_class(red, "DimReduc")
  emb <- Seurat::Embeddings(red)
  expect_equal(dim(emb), c(nrow(X), 2L))
  expect_equal(rownames(emb), rownames(X))
  expect_equal(colnames(emb), c("LocalMAP_1", "LocalMAP_2"))
  expect_true(all(is.finite(emb)))
})

test_that("RunLocalMAP.default honors n_components", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  red <- RunLocalMAP(
    object       = X,
    n_components = 3L,
    num_iters    = 50L,
    verbose      = FALSE,
    seed.use     = 11L
  )
  expect_equal(ncol(Seurat::Embeddings(red)), 3L)
  expect_equal(colnames(Seurat::Embeddings(red)),
               c("LocalMAP_1", "LocalMAP_2", "LocalMAP_3"))
})

test_that("RunLocalMAP.default is reproducible under the same seed", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  r1 <- RunLocalMAP(X, num_iters = 50L, verbose = FALSE, seed.use = 42L)
  r2 <- RunLocalMAP(X, num_iters = 50L, verbose = FALSE, seed.use = 42L)
  expect_equal(Seurat::Embeddings(r1), Seurat::Embeddings(r2))
})

test_that("RunLocalMAP.default accepts a custom reduction.key", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  red <- RunLocalMAP(
    object        = X,
    num_iters     = 50L,
    verbose       = FALSE,
    seed.use      = 11L,
    reduction.key = "LM_"
  )
  expect_equal(colnames(Seurat::Embeddings(red)), c("LM_1", "LM_2"))
})

test_that("RunLocalMAP.default preserves cluster structure on a toy input", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix(n_per_cluster = 40L)
  labels <- rep(c("A", "B"), each = 40L)

  red <- RunLocalMAP(X, num_iters = 100L, verbose = FALSE, seed.use = 11L)
  emb <- Seurat::Embeddings(red)

  ca <- colMeans(emb[labels == "A", , drop = FALSE])
  cb <- colMeans(emb[labels == "B", , drop = FALSE])
  wa <- mean(apply(emb[labels == "A", , drop = FALSE], 1,
                   function(r) sqrt(sum((r - ca)^2))))
  wb <- mean(apply(emb[labels == "B", , drop = FALSE], 1,
                   function(r) sqrt(sum((r - cb)^2))))
  between <- sqrt(sum((ca - cb)^2))
  expect_gt(between, max(wa, wb))
})

test_that("RunLocalMAP.default forwards low_dist_thres to pacmapr::localmap", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  # Two very different thresholds should yield non-identical embeddings.
  r1 <- RunLocalMAP(X, num_iters = 50L, verbose = FALSE, seed.use = 11L,
                    low_dist_thres = 1)
  r2 <- RunLocalMAP(X, num_iters = 50L, verbose = FALSE, seed.use = 11L,
                    low_dist_thres = 100)
  expect_false(isTRUE(all.equal(Seurat::Embeddings(r1),
                                Seurat::Embeddings(r2))))
})

test_that("RunLocalMAP.Seurat wires features -> matrix path correctly", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  feats <- Seurat::VariableFeatures(obj)

  obj <- RunLocalMAP(
    object    = obj,
    features  = feats,
    num_iters = 50L,
    verbose   = FALSE,
    seed.use  = 11L
  )

  expect_true("localmap" %in% Seurat::Reductions(obj))
  emb <- Seurat::Embeddings(obj[["localmap"]])
  expect_equal(nrow(emb), ncol(obj))
  expect_equal(rownames(emb), colnames(obj))
  expect_equal(ncol(emb), 2L)
})

test_that("RunLocalMAP.Seurat wires dims -> matrix path correctly", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  obj <- RunLocalMAP(
    object    = obj,
    reduction = "pca",
    dims      = 1:5,
    num_iters = 50L,
    verbose   = FALSE,
    seed.use  = 11L
  )
  emb <- Seurat::Embeddings(obj[["localmap"]])
  expect_equal(dim(emb), c(ncol(obj), 2L))
  expect_equal(rownames(emb), colnames(obj))
})

test_that("RunLocalMAP.Seurat errors when neither dims nor features supplied", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  expect_error(
    RunLocalMAP(object = obj, num_iters = 50L, verbose = FALSE),
    regexp = "dims|features"
  )
})

test_that("RunLocalMAP.Seurat errors when fewer features than n_components", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  feats <- Seurat::VariableFeatures(obj)[1]
  expect_error(
    RunLocalMAP(
      object       = obj,
      features     = feats,
      n_components = 2L,
      num_iters    = 50L,
      verbose      = FALSE
    ),
    regexp = "more features than n_components"
  )
})

test_that("Missing pacmapr yields an informative error", {
  if (requireNamespace("pacmapr", quietly = TRUE)) {
    skip("pacmapr is installed; the missing-package guard is not exercised here")
  }
  X <- make_toy_matrix(n_per_cluster = 10L, p = 4L)
  expect_error(
    RunLocalMAP(X, num_iters = 20L, verbose = FALSE),
    regexp = "pacmapr"
  )
})
