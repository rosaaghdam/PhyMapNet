#' Prepare inputs for PhyMapNet
#'
#' Computes a patristic distance matrix from a phylogenetic tree, or validates
#' a supplied phylogenetic distance matrix, then aligns taxa to the OTU table.
#'
#' @param otu A samples x taxa count/abundance matrix (rownames=samples, colnames=taxa).
#' @param tree A `phylo` object with tip labels as taxa names, or a named,
#'   symmetric taxa x taxa phylogenetic distance matrix.
#' @param prune Logical; if TRUE and `tree` is a `phylo` object, prunes tree
#'   tips not found in OTU before distances are computed. Input taxa are always
#'   aligned to their shared set.
#' @return A list with otu (aligned), dist (aligned), taxa.
#' @export
phymapnet_prepare_inputs <- function(otu, tree, prune = TRUE) {
  .check_matrix(otu, "otu")

  if (inherits(tree, "phylo")) {
    if (prune) {
      keep <- intersect(tree$tip.label, colnames(otu))
      if (length(keep) < 2) stop("After pruning, fewer than 2 taxa remain.", call. = FALSE)
      otu_only <- setdiff(colnames(otu), keep)
      tree_only <- setdiff(tree$tip.label, keep)
      if (length(otu_only) || length(tree_only)) {
        warning(
          sprintf(
            "Retaining %d shared taxa; dropped %d OTU-only taxa and %d tree-only taxa.",
            length(keep), length(otu_only), length(tree_only)
          ),
          call. = FALSE
        )
      }
      tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
    }

    dist <- as.matrix(ape::cophenetic.phylo(tree))
  } else {
    dist <- .check_dist(tree)
  }

  .align_otu_dist(otu, dist, warn = !(inherits(tree, "phylo") && prune))
}
