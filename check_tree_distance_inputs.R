# Check that tree and stored phylogenetic distance inputs agree for an OTU table.
#
# Usage:
#   library(phymapnet)
#   source("check_tree_distance_inputs.R")
#   check_tree_distance_inputs(otu, tree, DIS)

check_tree_distance_inputs <- function(otu, tree, DIS, tolerance = 1e-10) {
  from_tree <- suppressWarnings(phymapnet_prepare_inputs(as.matrix(otu), tree, prune = TRUE))
  from_dist <- suppressWarnings(phymapnet_prepare_inputs(as.matrix(otu), as.matrix(DIS)))

  same_taxa <- identical(from_tree$taxa, from_dist$taxa)
  max_diff <- if (same_taxa) max(abs(from_tree$dist - from_dist$dist)) else NA_real_

  data.frame(
    n_taxa = length(from_tree$taxa),
    same_taxa = same_taxa,
    max_distance_difference = max_diff,
    same_distance = same_taxa && is.finite(max_diff) && max_diff <= tolerance
  )
}
