test_that("phymapnet_prepare_inputs aligns taxa", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(a:1,b:1,c:1);")
  otu <- matrix(1, nrow=3, ncol=2, dimnames=list(paste0("s",1:3), c("a","b")))
  expect_warning(
    prep <- phymapnet_prepare_inputs(otu, tree, prune = TRUE),
    "dropped 0 OTU-only taxa and 1 tree-only taxa"
  )
  expect_equal(colnames(prep$otu), rownames(prep$dist))
  expect_equal(colnames(prep$otu), prep$taxa)
})

test_that("tree pruning drops nonshared taxa and follows OTU column order", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(a:1,b:1,c:1);")
  otu <- matrix(
    1,
    nrow = 3,
    ncol = 3,
    dimnames = list(paste0("s", 1:3), c("b", "a", "x"))
  )

  expect_warning(
    prep <- phymapnet_prepare_inputs(otu, tree, prune = TRUE),
    "dropped 1 OTU-only taxa and 1 tree-only taxa"
  )

  expect_identical(prep$taxa, c("b", "a"))
  expect_identical(colnames(prep$otu), c("b", "a"))
  expect_identical(rownames(prep$dist), c("b", "a"))
  expect_identical(colnames(prep$dist), c("b", "a"))
})

test_that("named phylogenetic distance matrices are supported and aligned", {
  otu <- matrix(
    1,
    nrow = 3,
    ncol = 3,
    dimnames = list(paste0("s", 1:3), c("b", "a", "x"))
  )
  dist <- matrix(
    c(0, 2, 3, 2, 0, 1, 3, 1, 0),
    nrow = 3,
    dimnames = list(c("c", "a", "b"), c("c", "a", "b"))
  )

  expect_warning(
    prep <- phymapnet_prepare_inputs(otu, dist),
    "dropped 1 OTU-only taxa and 1 tree/distance-only taxa"
  )

  expect_identical(prep$taxa, c("b", "a"))
  expect_identical(rownames(prep$dist), c("b", "a"))
  expect_equal(unname(prep$dist), matrix(c(0, 1, 1, 0), nrow = 2))
})

test_that("invalid distance matrices fail clearly", {
  otu <- matrix(1, nrow = 3, ncol = 2, dimnames = list(paste0("s", 1:3), c("a", "b")))
  nonsymmetric <- matrix(c(0, 1, 2, 0), nrow = 2, dimnames = list(c("a", "b"), c("a", "b")))
  unnamed <- matrix(c(0, 1, 1, 0), nrow = 2)

  expect_error(phymapnet_prepare_inputs(otu, nonsymmetric), "must be symmetric")
  expect_error(phymapnet_prepare_inputs(otu, unnamed), "taxa names")
})
