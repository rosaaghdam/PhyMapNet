################################################################################
# Real-data local test for phymapnet 0.1.3
#
# Run from the repository root:
#   Rscript test_local_phymapnet_0.1.3_real_data.R
#
# This installs the package archive into a temporary library and evaluates the
# real filtered study datasets using a small parameter grid. It does NOT run
# the full paper ensemble.
#
# Tested datasets:
#   - smoking:   60 samples x 174 taxa
#   - caffeine:  98 samples x 299 taxa
#   - hmp_stool: 187 samples x 312 taxa
#
# For each dataset the script:
#   1. aligns OTU data to the source tree;
#   2. aligns OTU data to the stored distance matrix;
#   3. verifies the tree-derived and stored distances agree;
#   4. runs phymapnet_fit() with tree and distance inputs;
#   5. runs phymapnet_reliability() with distance input on four models.
################################################################################

rm(list = ls())

package_tarball <- "phymapnet_0.1.3.tar.gz"
expected_version <- "0.1.3"

if (!file.exists(package_tarball)) {
  stop("Package tarball not found: ", package_tarball, call. = FALSE)
}

required_files <- c(
  "data/filter/smoking_otu_filtered.csv",
  "data/filter/smoking_dist_matrix.csv",
  "data/filter/caffeine_otu_filtered.csv",
  "data/filter/caffeine_dist_matrix.csv",
  "data/filter/hmp_stool_otu_filtered.csv",
  "data/filter/hmp_stool_dist_matrix.csv",
  "data/original_data/data_smoking.RData",
  "data/original_data/data_caff.RData",
  "data/original_data/data_hmp.RData"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required study files: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

test_lib <- tempfile("phymapnet_real_data_library_")
dir.create(test_lib, recursive = TRUE)

cat("Installing", package_tarball, "into temporary library...\n")
install.packages(
  package_tarball,
  repos = NULL,
  type = "source",
  lib = test_lib,
  quiet = TRUE
)

library(phymapnet, lib.loc = test_lib)

installed_version <- as.character(packageVersion("phymapnet", lib.loc = test_lib))
stopifnot(identical(installed_version, expected_version))

cat("Installed package version:", installed_version, "\n")
cat("Running small real-data tests; full paper ensemble is not run.\n\n")

dataset_info <- list(
  smoking = list(rdata = "data_smoking.RData", object = "data_smoking"),
  caffeine = list(rdata = "data_caff.RData", object = "data_caff"),
  hmp_stool = list(rdata = "data_hmp.RData", object = "data_hmp")
)

# One normalization and one kernel; four combinations of alpha and k.
fit_parameters <- list(
  alpha = 0.05,
  k = 3,
  epsilon1 = 0.1,
  epsilon2 = 0.1,
  kernel = "gaussian",
  th_sparsity = 0.95,
  normalization = "log"
)

reliability_parameters <- list(
  th_fixed = 0.95,
  alpha_range = c(0.05, 0.10),
  k_range = c(2, 3),
  epsilon1_range = 0.1,
  epsilon2_range = 0.1,
  kernels = "gaussian",
  normalizations = "log",
  consensus_cut = 0.50,
  progress = FALSE
)

summary_rows <- list()

for (dataset in names(dataset_info)) {
  cat("Dataset:", dataset, "\n")

  otu <- as.matrix(
    read.csv(
      file.path("data", "filter", paste0(dataset, "_otu_filtered.csv")),
      row.names = 1,
      check.names = FALSE
    )
  )
  colnames(otu) <- gsub("^X", "", colnames(otu))

  data_env <- new.env()
  load(
    file.path("data", "original_data", dataset_info[[dataset]]$rdata),
    envir = data_env
  )
  tree <- get(dataset_info[[dataset]]$object, envir = data_env)$tree

  DIS <- as.matrix(
    read.csv(
      file.path("data", "filter", paste0(dataset, "_dist_matrix.csv")),
      row.names = 1,
      check.names = FALSE
    )
  )

  prep_tree <- suppressWarnings(phymapnet_prepare_inputs(otu, tree, prune = TRUE))
  prep_dist <- suppressWarnings(phymapnet_prepare_inputs(otu, DIS))

  stopifnot(identical(prep_tree$taxa, prep_dist$taxa))
  max_distance_difference <- max(abs(prep_tree$dist - prep_dist$dist))
  if (max_distance_difference > 1e-10) {
    stop(
      dataset, ": stored distance matrix does not match tree-derived distances.",
      call. = FALSE
    )
  }

  fit_tree_elapsed <- system.time({
    fit_tree <- do.call(
      phymapnet_fit,
      c(list(otu = otu, tree = tree, prune_tree = TRUE), fit_parameters)
    )
  })[["elapsed"]]

  fit_dist_elapsed <- system.time({
    fit_dist <- do.call(
      phymapnet_fit,
      c(list(otu = otu, tree = DIS), fit_parameters)
    )
  })[["elapsed"]]

  stopifnot(identical(fit_tree$taxa, fit_dist$taxa))
  stopifnot(isTRUE(all.equal(fit_tree$kernel_mat, fit_dist$kernel_mat, tolerance = 1e-10)))
  stopifnot(identical(fit_tree$adjacency, fit_dist$adjacency))

  reliability_elapsed <- system.time({
    res <- do.call(
      phymapnet_reliability,
      c(list(otu = otu, tree = DIS), reliability_parameters)
    )
  })[["elapsed"]]

  stopifnot(res$N_models == 4)
  stopifnot(identical(rownames(res$rel_mat), prep_dist$taxa))
  stopifnot(all(is.finite(res$rel_mat)))

  summary_rows[[dataset]] <- data.frame(
    dataset = dataset,
    samples = nrow(otu),
    taxa = length(prep_dist$taxa),
    max_distance_difference = max_distance_difference,
    tree_fit_seconds = fit_tree_elapsed,
    distance_fit_seconds = fit_dist_elapsed,
    reliability_models = res$N_models,
    reliability_seconds = reliability_elapsed,
    max_reliability = max(res$rel_mat),
    stringsAsFactors = FALSE
  )

  cat("  samples x taxa:", nrow(otu), "x", length(prep_dist$taxa), "\n")
  cat("  tree and stored distance input agree: yes\n")
  cat("  reliability models fitted:", res$N_models, "\n")
  cat("  elapsed reliability seconds:", round(reliability_elapsed, 3), "\n\n")
}

summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL

cat("Summary:\n")
print(summary_table)
cat("\nAll real-data local tests passed for phymapnet ", installed_version, ".\n", sep = "")
