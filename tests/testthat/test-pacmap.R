# Tests for RunPaCMAP (R/pacmap.R), which now dispatches to the native
# pacmapr package instead of the Python pacmap module.

test_that("RunPaCMAP.default returns a DimReduc with expected shape", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  red <- RunPaCMAP(
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
  expect_equal(colnames(emb), c("PaCMAP_1", "PaCMAP_2"))
  expect_true(all(is.finite(emb)))
})

test_that("RunPaCMAP.default honors n_components", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  red <- RunPaCMAP(
    object       = X,
    n_components = 3L,
    num_iters    = 50L,
    verbose      = FALSE,
    seed.use     = 11L
  )
  expect_equal(ncol(Seurat::Embeddings(red)), 3L)
  expect_equal(colnames(Seurat::Embeddings(red)),
               c("PaCMAP_1", "PaCMAP_2", "PaCMAP_3"))
})

test_that("RunPaCMAP.default is reproducible under the same seed", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  r1 <- RunPaCMAP(X, num_iters = 50L, verbose = FALSE, seed.use = 42L)
  r2 <- RunPaCMAP(X, num_iters = 50L, verbose = FALSE, seed.use = 42L)
  expect_equal(Seurat::Embeddings(r1), Seurat::Embeddings(r2))
})

test_that("RunPaCMAP.default accepts a custom reduction.key", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix()
  red <- RunPaCMAP(
    object        = X,
    num_iters     = 50L,
    verbose       = FALSE,
    seed.use      = 11L,
    reduction.key = "PC_"
  )
  expect_equal(colnames(Seurat::Embeddings(red)), c("PC_1", "PC_2"))
})

test_that("RunPaCMAP.default preserves cluster structure on a toy input", {
  skip_without_pacmapr()
  skip_without_seurat()

  X <- make_toy_matrix(n_per_cluster = 40L)
  labels <- rep(c("A", "B"), each = 40L)

  red <- RunPaCMAP(X, num_iters = 100L, verbose = FALSE, seed.use = 11L)
  emb <- Seurat::Embeddings(red)

  ca <- colMeans(emb[labels == "A", , drop = FALSE])
  cb <- colMeans(emb[labels == "B", , drop = FALSE])
  wa <- mean(apply(emb[labels == "A", , drop = FALSE], 1,
                   function(r) sqrt(sum((r - ca)^2))))
  wb <- mean(apply(emb[labels == "B", , drop = FALSE], 1,
                   function(r) sqrt(sum((r - cb)^2))))
  between <- sqrt(sum((ca - cb)^2))

  # Between-cluster centroid distance should be well above the average
  # within-cluster spread if PaCMAP separated the two blobs at all.
  expect_gt(between, max(wa, wb))
})

test_that("RunPaCMAP.Seurat wires features -> matrix path correctly", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  feats <- Seurat::VariableFeatures(obj)

  obj <- RunPaCMAP(
    object    = obj,
    features  = feats,
    num_iters = 50L,
    verbose   = FALSE,
    seed.use  = 11L
  )

  expect_true("pacmap" %in% Seurat::Reductions(obj))
  emb <- Seurat::Embeddings(obj[["pacmap"]])
  expect_equal(nrow(emb), ncol(obj))
  expect_equal(rownames(emb), colnames(obj))
  expect_equal(ncol(emb), 2L)
})

test_that("RunPaCMAP.Seurat wires dims -> matrix path correctly", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  obj <- RunPaCMAP(
    object    = obj,
    reduction = "pca",
    dims      = 1:5,
    num_iters = 50L,
    verbose   = FALSE,
    seed.use  = 11L
  )

  emb <- Seurat::Embeddings(obj[["pacmap"]])
  expect_equal(dim(emb), c(ncol(obj), 2L))
  expect_equal(rownames(emb), colnames(obj))
})

test_that("RunPaCMAP.Seurat errors when neither dims nor features supplied", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  expect_error(
    RunPaCMAP(object = obj, num_iters = 50L, verbose = FALSE),
    regexp = "dims|features"
  )
})

test_that("RunPaCMAP.Seurat errors when fewer features than n_components", {
  skip_without_pacmapr()
  skip_without_seurat()

  obj <- make_toy_seurat()
  feats <- Seurat::VariableFeatures(obj)[1]  # only one feature
  expect_error(
    RunPaCMAP(
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
  # Simulate absence by shadowing requireNamespace within this test's env.
  # We only exercise the guard when pacmapr is really missing, otherwise
  # the wrapper would run for real.
  if (requireNamespace("pacmapr", quietly = TRUE)) {
    skip("pacmapr is installed; the missing-package guard is not exercised here")
  }
  X <- make_toy_matrix(n_per_cluster = 10L, p = 4L)
  expect_error(
    RunPaCMAP(X, num_iters = 20L, verbose = FALSE),
    regexp = "pacmapr"
  )
})
