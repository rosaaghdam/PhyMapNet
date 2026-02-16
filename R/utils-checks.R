.check_matrix <- function(x, name) {
  if (!is.matrix(x)) stop(name, " must be a matrix.", call. = FALSE)
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(name, " must have rownames and colnames.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_tree <- function(tree) {
  if (!inherits(tree, "phylo")) stop("tree must be an object of class 'phylo'.", call. = FALSE)
  invisible(TRUE)
}

.align_otu_dist <- function(otu, dist) {
  common <- intersect(colnames(otu), rownames(dist))
  if (length(common) < 2) stop("Not enough overlap between OTU taxa and distance taxa.", call. = FALSE)
  otu2  <- otu[, common, drop = FALSE]
  dist2 <- dist[common, common, drop = FALSE]
  list(otu = otu2, dist = dist2, taxa = common)
}
