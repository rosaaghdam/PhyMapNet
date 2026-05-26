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

.check_dist <- function(dist) {
  if (!is.matrix(dist) || !is.numeric(dist)) {
    stop("Distance input must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(dist) != ncol(dist)) {
    stop("Distance matrix must be square.", call. = FALSE)
  }
  if (is.null(rownames(dist)) || is.null(colnames(dist))) {
    stop("Distance matrix must have taxa names in rownames and colnames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(dist)) || anyDuplicated(colnames(dist))) {
    stop("Distance matrix taxa names must be unique.", call. = FALSE)
  }
  if (!setequal(rownames(dist), colnames(dist))) {
    stop("Distance matrix rownames and colnames must contain the same taxa.", call. = FALSE)
  }

  dist <- dist[rownames(dist), rownames(dist), drop = FALSE]
  if (any(!is.finite(dist))) stop("Distance matrix must contain only finite values.", call. = FALSE)
  if (any(dist < 0)) stop("Distance matrix cannot contain negative values.", call. = FALSE)
  if (!isTRUE(all.equal(dist, t(dist), tolerance = 1e-8))) {
    stop("Distance matrix must be symmetric.", call. = FALSE)
  }
  if (any(abs(diag(dist)) > 1e-8)) {
    stop("Distance matrix diagonal must be zero.", call. = FALSE)
  }
  dist
}

.align_otu_dist <- function(otu, dist, warn = TRUE) {
  common <- intersect(colnames(otu), rownames(dist))
  if (length(common) < 2) stop("Not enough overlap between OTU taxa and distance taxa.", call. = FALSE)
  otu_only <- setdiff(colnames(otu), common)
  dist_only <- setdiff(rownames(dist), common)
  if (warn && (length(otu_only) || length(dist_only))) {
    warning(
      sprintf(
        "Retaining %d shared taxa; dropped %d OTU-only taxa and %d tree/distance-only taxa.",
        length(common), length(otu_only), length(dist_only)
      ),
      call. = FALSE
    )
  }
  otu2  <- otu[, common, drop = FALSE]
  dist2 <- dist[common, common, drop = FALSE]
  list(otu = otu2, dist = dist2, taxa = common)
}
