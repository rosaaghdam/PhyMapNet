################################################################################
# Local installation and smoke test for phymapnet 0.1.3
#
# Run from the repository root:
#   Rscript test_local_phymapnet_0.1.3.R
#
# This script installs the source archive into a temporary R library so it does
# not overwrite an existing personal installation. It then verifies the revised
# public normalization options and runs small tree- and distance-input tests.
################################################################################

rm(list = ls())

package_tarball <- "phymapnet_0.1.3.tar.gz"
expected_version <- "0.1.3"

if (!file.exists(package_tarball)) {
  stop("Package tarball not found: ", package_tarball, call. = FALSE)
}

test_lib <- tempfile("phymapnet_test_library_")
dir.create(test_lib, recursive = TRUE)

cat("Installing", package_tarball, "into temporary library:\n", test_lib, "\n\n")

install.packages(
  package_tarball,
  repos = NULL,
  type = "source",
  lib = test_lib,
  quiet = TRUE
)

library(phymapnet, lib.loc = test_lib)
library(ape)

installed_version <- as.character(packageVersion("phymapnet", lib.loc = test_lib))
if (!identical(installed_version, expected_version)) {
  stop(
    "Installed version is ", installed_version,
    "; expected ", expected_version, ".",
    call. = FALSE
  )
}

normalization_defaults <- eval(formals(phymapnet_reliability)$normalizations)
expected_normalizations <- c("log", "clr", "tss")
if (!identical(normalization_defaults, expected_normalizations)) {
  stop(
    "Unexpected reliability normalization defaults: ",
    paste(normalization_defaults, collapse = ", "),
    call. = FALSE
  )
}

cat("Installed package version:", installed_version, "\n")
cat("Reliability normalizations:", paste(normalization_defaults, collapse = ", "), "\n\n")

# A small deterministic example with one OTU-only and one tree-only taxon.
otu <- matrix(
  c(
    2, 0, 4, 1,
    3, 5, 0, 2,
    1, 4, 2, 3
  ),
  nrow = 4,
  dimnames = list(paste0("sample_", 1:4), c("b", "a", "x"))
)

tree <- ape::read.tree(text = "(a:1,b:1,c:1);")

cat("Testing tree input and pruning/alignment...\n")
prep_tree <- withCallingHandlers(
  phymapnet_prepare_inputs(otu, tree, prune = TRUE),
  warning = function(w) {
    cat("  Expected alignment warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }
)

stopifnot(identical(prep_tree$taxa, c("b", "a")))
stopifnot(identical(colnames(prep_tree$otu), rownames(prep_tree$dist)))
stopifnot(identical(rownames(prep_tree$dist), colnames(prep_tree$dist)))

fit_tree <- suppressWarnings(
  phymapnet_fit(
    otu = otu,
    tree = tree,
    alpha = 0.5,
    k = 2,
    epsilon1 = 0.1,
    epsilon2 = 0.1,
    kernel = "gaussian",
    th_sparsity = 0.90,
    normalization = "log",
    prune_tree = TRUE
  )
)

cat("  Tree workflow retained taxa:", paste(fit_tree$taxa, collapse = ", "), "\n\n")

# Named phylogenetic distance matrix with the same shared taxa plus an unused taxon.
DIS <- matrix(
  c(
    0, 2, 3,
    2, 0, 1,
    3, 1, 0
  ),
  nrow = 3,
  dimnames = list(c("c", "a", "b"), c("c", "a", "b"))
)

cat("Testing distance-matrix input and alignment...\n")
prep_dist <- withCallingHandlers(
  phymapnet_prepare_inputs(otu, DIS),
  warning = function(w) {
    cat("  Expected alignment warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }
)

stopifnot(identical(prep_dist$taxa, c("b", "a")))
stopifnot(identical(rownames(prep_dist$dist), c("b", "a")))

fit_dist <- suppressWarnings(
  phymapnet_fit(
    otu = otu,
    tree = DIS,
    alpha = 0.5,
    k = 2,
    epsilon1 = 0.1,
    epsilon2 = 0.1,
    kernel = "gaussian",
    th_sparsity = 0.90,
    normalization = "log"
  )
)

res_dist <- suppressWarnings(
  phymapnet_reliability(
    otu = otu,
    tree = DIS,
    th_fixed = 0.90,
    alpha_range = 0.5,
    k_range = 2,
    epsilon1_range = 0.1,
    epsilon2_range = 0.1,
    kernels = "gaussian",
    normalizations = "log",
    consensus_cut = 0.50,
    progress = FALSE
  )
)

stopifnot(identical(fit_dist$taxa, c("b", "a")))
stopifnot(identical(rownames(res_dist$rel_mat), c("b", "a")))
stopifnot(res_dist$N_models == 1)

cat("  Distance workflow retained taxa:", paste(fit_dist$taxa, collapse = ", "), "\n")
cat("  Reliability models fitted:", res_dist$N_models, "\n\n")

cat("All local installation and smoke tests passed for phymapnet ", installed_version, ".\n", sep = "")
