#' Prepare inputs for PhyMapNet
#'
#' Computes the patristic distance matrix from a phylogenetic tree and aligns
#' taxa between the OTU table and the tree.
#'
#' @param otu A samples x taxa count/abundance matrix (rownames=samples, colnames=taxa).
#' @param tree A phylo object with tip labels as taxa names.
#' @param prune Logical; if TRUE, prunes tree tips not found in OTU.
#' @return A list with otu (aligned), dist (aligned), taxa.
#' @export
phymapnet_prepare_inputs <- function(otu, tree, prune = TRUE) {
  .check_matrix(otu, "otu")
  .check_tree(tree)

  if (prune) {
    keep <- intersect(tree$tip.label, colnames(otu))
    if (length(keep) < 2) stop("After pruning, fewer than 2 taxa remain.", call. = FALSE)
    tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  }

  dist <- ape::cophenetic.phylo(tree)
  dist <- as.matrix(dist)

  .align_otu_dist(otu, dist)
}
