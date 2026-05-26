################################################################################
# Filtering and phylogenetic-distance preparation for the three study datasets
#
# Run from the repository root:
#   Rscript 1-filtering/1_run_filtering.R
#
# Inputs:
#   data/original_data/data_smoking.RData
#   data/original_data/data_caff.RData
#   data/original_data/data_hmp.RData
#
# Outputs:
#   data/filter/<dataset>_otu_filtered.csv
#   data/filter/<dataset>_dist_matrix.csv
#
# To test without overwriting published filtered inputs:
#   PHYMAPNET_FILTER_OUTPUT_DIR=/tmp/phymapnet_filter_test \
#     Rscript 1-filtering/1_run_filtering.R
################################################################################

rm(list = ls())

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The ape package is required to run filtering.", call. = FALSE)
}

out_dir <- Sys.getenv("PHYMAPNET_FILTER_OUTPUT_DIR", unset = "data/filter")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# This reproduces the filtering calculation used to prepare the analysis data.
filter_dataset <- function(data_list, dataset, zero_threshold, out_dir) {
  otu <- as.matrix(data_list$otu.tab)
  tree <- data_list$tree

  n_samples_raw <- nrow(otu)
  n_taxa_raw <- ncol(otu)

  # Retain samples containing at least 15 observed taxa.
  otu <- otu[rowSums(otu > 0) >= 15, , drop = FALSE]

  # Preserve the original analysis definition of the zero proportion:
  # the denominator is the number of samples in the raw input table.
  zero_proportion <- colSums(otu == 0) / n_samples_raw
  keep_by_zero <- names(zero_proportion)[zero_proportion < zero_threshold]

  nonzero_median <- apply(otu, 2, function(x) median(x[x != 0]))
  nonzero_median[is.na(nonzero_median)] <- 0
  keep_by_median <- names(nonzero_median)[nonzero_median > 0]

  keep_taxa <- intersect(keep_by_zero, keep_by_median)
  shared_taxa <- intersect(tree$tip.label, keep_taxa)
  if (length(shared_taxa) < 2) {
    stop("Fewer than two filtered taxa match the tree for ", dataset, ".", call. = FALSE)
  }

  tree_filtered <- ape::keep.tip(tree, shared_taxa)
  otu_filtered <- otu[, tree_filtered$tip.label, drop = FALSE]
  dist_filtered <- as.matrix(ape::cophenetic.phylo(tree_filtered))

  stopifnot(identical(colnames(otu_filtered), rownames(dist_filtered)))
  stopifnot(identical(rownames(dist_filtered), colnames(dist_filtered)))

  write.csv(
    otu_filtered,
    file.path(out_dir, paste0(dataset, "_otu_filtered.csv")),
    row.names = TRUE
  )
  write.csv(
    dist_filtered,
    file.path(out_dir, paste0(dataset, "_dist_matrix.csv")),
    row.names = TRUE
  )

  data.frame(
    dataset = dataset,
    zero_threshold = zero_threshold,
    samples_raw = n_samples_raw,
    taxa_raw = n_taxa_raw,
    samples_filtered = nrow(otu_filtered),
    taxa_filtered = ncol(otu_filtered),
    stringsAsFactors = FALSE
  )
}

load("data/original_data/data_smoking.RData")
load("data/original_data/data_caff.RData")
load("data/original_data/data_hmp.RData")

results <- rbind(
  filter_dataset(data_smoking, "smoking", 0.90, out_dir),
  filter_dataset(data_caff, "caffeine", 0.75, out_dir),
  filter_dataset(data_hmp, "hmp_stool", 0.65, out_dir)
)

cat("Filtered data written to:", normalizePath(out_dir), "\n\n")
print(results, row.names = FALSE)
