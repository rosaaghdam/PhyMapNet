test_that("kernel_from_dist works", {
  dist <- matrix(c(0,1,1,0), 2, 2, dimnames = list(c("a","b"),c("a","b")))
  K1 <- kernel_from_dist(dist, alpha = 1, kernel = "gaussian")
  K2 <- kernel_from_dist(dist, alpha = 1, kernel = "laplacian")
  expect_equal(dim(K1), c(2,2))
  expect_equal(dim(K2), c(2,2))
  expect_true(all(diag(K1) == 1))
  expect_true(all(diag(K2) == 1))
})

test_that("public normalization defaults contain revised methods only", {
  expect_identical(eval(formals(phymapnet_fit)$normalization), c("log", "clr", "tss"))
  expect_identical(eval(formals(phymapnet_reliability)$normalizations), c("log", "clr", "tss"))
})

test_that("fit and reliability accept a supplied distance matrix", {
  otu <- matrix(
    c(1, 4, 2, 3, 2, 1, 4, 3, 5, 2, 3, 1),
    nrow = 4,
    dimnames = list(paste0("s", 1:4), c("a", "b", "c"))
  )
  dist <- matrix(
    c(0, 1, 2, 1, 0, 1, 2, 1, 0),
    nrow = 3,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )

  fit <- phymapnet_fit(
    otu, dist,
    alpha = 1, k = 2, epsilon1 = 0.1, epsilon2 = 0.1,
    kernel = "gaussian", normalization = "log"
  )
  res <- phymapnet_reliability(
    otu, dist,
    alpha_range = 1, k_range = 2,
    epsilon1_range = 0.1, epsilon2_range = 0.1,
    kernels = "gaussian", normalizations = "log", progress = FALSE
  )

  expect_identical(fit$taxa, c("a", "b", "c"))
  expect_equal(fit$dist, dist)
  expect_equal(res$N_models, 1)
  expect_identical(rownames(res$rel_mat), c("a", "b", "c"))
})

test_that("GMPR is not a supported revised-analysis normalization", {
  otu <- matrix(1:12, nrow = 4, dimnames = list(paste0("s", 1:4), c("a", "b", "c")))
  dist <- matrix(
    c(0, 1, 2, 1, 0, 1, 2, 1, 0),
    nrow = 3,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )

  expect_error(
    phymapnet_fit(otu, dist, normalization = "gmpr"),
    "'arg' should be one of"
  )
  expect_error(
    phymapnet_reliability(
      otu, dist,
      alpha_range = 1, k_range = 2,
      epsilon1_range = 0.1, epsilon2_range = 0.1,
      normalizations = "gmpr", progress = FALSE
    ),
    "subset of"
  )
})
