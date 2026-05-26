################################################################################
# Step 1: Generate parameter-specific reliability networks for sensitivity study
#
# IMPORTANT:
#   The complete calculation is time-consuming. Precomputed outputs used for
#   the manuscript are included under result/reliable_score_all/. Do not rerun
#   the full calculation unless regeneration is required.
#
# Full regeneration for all datasets, from the repository root:
#   Rscript 2-sensitivity_analysis/generate_sensitivity_reliability_networks.R
#
# Run one dataset only:
#   PHYMAPNET_DATASETS=caffeine \
#     Rscript 2-sensitivity_analysis/generate_sensitivity_reliability_networks.R
#
# Fast validation run (one dataset and one model; does not overwrite results):
#   PHYMAPNET_SMOKE_TEST=true PHYMAPNET_DATASETS=caffeine \
#     Rscript 2-sensitivity_analysis/generate_sensitivity_reliability_networks.R
#
# Outputs for full regeneration:
#   result/reliable_score_all/<dataset>/*.csv
################################################################################

rm(list = ls())

if (!requireNamespace("phymapnet", quietly = TRUE)) {
  stop("Install the phymapnet package before running this script.", call. = FALSE)
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The ape package is required to run this script.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

repo_root <- "."
available_datasets <- c("smoking", "caffeine", "hmp_stool")
requested <- Sys.getenv("PHYMAPNET_DATASETS", unset = paste(available_datasets, collapse = ","))
datasets <- trimws(strsplit(requested, ",", fixed = TRUE)[[1]])
if (!length(datasets) || any(!datasets %in% available_datasets)) {
  stop("PHYMAPNET_DATASETS must contain: smoking, caffeine, and/or hmp_stool.", call. = FALSE)
}

smoke_test <- identical(tolower(Sys.getenv("PHYMAPNET_SMOKE_TEST", unset = "false")), "true")
output_root <- if (smoke_test) {
  Sys.getenv("PHYMAPNET_SMOKE_OUTPUT_DIR", unset = file.path(tempdir(), "phymapnet_sensitivity_smoke"))
} else {
  file.path(repo_root, "result", "reliable_score_all")
}

# These fixed settings reproduce the parameter groups used for manuscript
# Figure 2 and Supplementary Tables S3/S4. Ten alpha settings result in
# 3 + 2 + 10 + 9 + 10 + 10 = 44 parameter-specific comparisons per dataset.
alpha_range <- seq(0.01, 0.10, by = 0.01)
k_range <- 2:10
epsilon_full <- seq(0, 1, by = 0.1)
kernels <- c("gaussian", "laplacian")
normalizations <- c("log", "clr", "tss")
epsilon_windows <- lapply(0:9, function(i) c(i / 10, (i + 1) / 10))

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

load_dataset <- function(dataset) {
  otu <- as.matrix(
    read.csv(
      file.path(repo_root, "data", "filter", paste0(dataset, "_otu_filtered.csv")),
      row.names = 1,
      check.names = FALSE
    )
  )
  colnames(otu) <- gsub("^X", "", colnames(otu))

  input_info <- switch(
    dataset,
    smoking = list(file = "data_smoking.RData", object = "data_smoking"),
    caffeine = list(file = "data_caff.RData", object = "data_caff"),
    hmp_stool = list(file = "data_hmp.RData", object = "data_hmp")
  )

  env <- new.env()
  load(file.path(repo_root, "data", "original_data", input_info$file), envir = env)
  tree <- get(input_info$object, envir = env)$tree

  shared <- intersect(colnames(otu), tree$tip.label)
  if (length(shared) < 2) stop("Fewer than two OTU taxa match the tree for ", dataset, ".", call. = FALSE)
  tree <- ape::keep.tip(tree, shared)
  otu <- otu[, tree$tip.label, drop = FALSE]

  stopifnot(identical(colnames(otu), tree$tip.label))
  list(otu = otu, tree = tree)
}

save_edges <- function(result, dataset, name) {
  edges <- result$edge_list
  colnames(edges)[1:3] <- c("taxa1", "taxa2", "value")
  dataset_dir <- file.path(output_root, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(edges, file.path(dataset_dir, paste0(name, ".csv")), row.names = FALSE)
}

fit_ensemble <- function(dat, alpha, k, epsilon1, epsilon2, kernel, normalization) {
  phymapnet::phymapnet_reliability(
    otu = dat$otu,
    tree = dat$tree,
    th_fixed = 0.95,
    alpha_range = alpha,
    k_range = k,
    epsilon1_range = epsilon1,
    epsilon2_range = epsilon2,
    kernels = kernel,
    normalizations = normalization,
    consensus_cut = 0.65,
    prune_tree = TRUE,
    progress = !smoke_test
  )
}

# -------------------------------------------------------------------------
# Fast code validation run: one complete dataset, one model
# -------------------------------------------------------------------------

if (smoke_test) {
  dataset <- datasets[1]
  dat <- load_dataset(dataset)
  result <- fit_ensemble(
    dat,
    alpha = 0.05,
    k = 2,
    epsilon1 = 0.1,
    epsilon2 = 0.1,
    kernel = "gaussian",
    normalization = "log"
  )
  save_edges(result, dataset, "smoke_test")
  cat("Smoke test completed:", dataset, nrow(dat$otu), "samples x", ncol(dat$otu), "taxa\n")
  cat("Models fitted:", result$N_models, "\n")
  cat("Output:", file.path(output_root, dataset, "smoke_test.csv"), "\n")
  quit(save = "no")
}

# -------------------------------------------------------------------------
# Full generation workflow
# -------------------------------------------------------------------------

for (dataset in datasets) {
  dat <- load_dataset(dataset)
  cat("\nDataset:", dataset, "|", nrow(dat$otu), "samples x", ncol(dat$otu), "taxa\n")

  cat("  Baseline ensemble\n")
  result <- fit_ensemble(dat, alpha_range, k_range, epsilon_full, epsilon_full,
                         kernels, normalizations)
  save_edges(result, dataset, "baseline")

  for (normalization in normalizations) {
    cat("  Normalization:", normalization, "\n")
    result <- fit_ensemble(dat, alpha_range, k_range, epsilon_full, epsilon_full,
                           kernels, normalization)
    save_edges(result, dataset, paste0("norm_", normalization))
  }

  for (kernel in kernels) {
    cat("  Kernel:", kernel, "\n")
    result <- fit_ensemble(dat, alpha_range, k_range, epsilon_full, epsilon_full,
                           kernel, normalizations)
    save_edges(result, dataset, paste0("kernel_", kernel))
  }

  for (alpha in alpha_range) {
    cat("  Alpha:", alpha, "\n")
    result <- fit_ensemble(dat, alpha, k_range, epsilon_full, epsilon_full,
                           kernels, normalizations)
    save_edges(result, dataset, paste0("alpha_", alpha))
  }

  for (k in k_range) {
    cat("  k:", k, "\n")
    result <- fit_ensemble(dat, alpha_range, k, epsilon_full, epsilon_full,
                           kernels, normalizations)
    save_edges(result, dataset, paste0("k_", k))
  }

  for (i in seq_along(epsilon_windows)) {
    window <- epsilon_windows[[i]]
    epsilon_range <- seq(window[1], window[2], by = 0.01)
    cat("  epsilon1 window:", i, "\n")
    result <- fit_ensemble(dat, alpha_range, k_range, epsilon_range, epsilon_full,
                           kernels, normalizations)
    save_edges(result, dataset, paste0("eps1_", i))
  }

  for (i in seq_along(epsilon_windows)) {
    window <- epsilon_windows[[i]]
    epsilon_range <- seq(window[1], window[2], by = 0.01)
    cat("  epsilon2 window:", i, "\n")
    result <- fit_ensemble(dat, alpha_range, k_range, epsilon_full, epsilon_range,
                           kernels, normalizations)
    save_edges(result, dataset, paste0("eps2_", i))
  }
}

cat("\nAll sensitivity network files generated in:", normalizePath(output_root), "\n")
