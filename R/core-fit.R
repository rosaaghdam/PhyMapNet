#' Internal PhyMapNet precision MAP estimator
#' @keywords internal
phymapnet_precision <- function(Y, nue, epsilon1, epsilon2, C) {
  p <- ncol(Y)
  n <- nrow(Y)

  S1 <- stats::cov(Y)
  V <- sqrt(diag(S1))
  V_diag <- diag(V)

  omega_hat <- V_diag %*% C %*% V_diag
  omega_hat <- omega_hat + diag(epsilon1, p)

  omega_star <- ((n/(n+nue))*S1) + ((nue/(nue+n))*omega_hat)
  nue_star <- nue + n

  IB <- as.matrix((nue_star) * omega_star) + epsilon2
  B <- solve(IB)

  mode_star <- (nue_star - p - 1) * B
  mode_star
}

#' Quantile sparsification
#' @keywords internal
sparse_quantile <- function(Isigma_star, quantile_level) {
  upper_vals <- abs(Isigma_star[upper.tri(Isigma_star)])
  th <- stats::quantile(upper_vals, quantile_level, names = FALSE, na.rm = TRUE)
  Isigma_sparse <- ifelse(abs(Isigma_star) < th, 0, Isigma_star)
  theta <- ifelse(Isigma_sparse != 0, 1L, 0L)
  list(Isigma_sparse = Isigma_sparse, theta = theta, threshold = th)
}

#' Fit a single PhyMapNet model
#'
#' @param otu samples x taxa matrix.
#' @param tree phylo tree with tips matching taxa.
#' @param alpha kernel bandwidth (>0).
#' @param k neighborhood scaling (integer >= 1). Uses K_neighbors = k * p internally.
#' @param epsilon1 diagonal jitter for omega_hat.
#' @param epsilon2 jitter for IB.
#' @param kernel "gaussian" or "laplacian".
#' @param th_sparsity quantile level for sparsification (e.g., 0.95).
#' @param normalization "log","gmpr","clr","tss".
#' @param prune_tree prune tree tips not in OTU.
#' @return A list with precision_map, adjacency, threshold, taxa, dist, kernel_mat.
#' @export
phymapnet_fit <- function(
  otu, tree,
  alpha = 0.05,
  k = 5,
  epsilon1 = 0,
  epsilon2 = 0,
  kernel = c("gaussian","laplacian"),
  th_sparsity = 0.95,
  normalization = c("log","gmpr","clr","tss"),
  prune_tree = TRUE
) {
  kernel <- match.arg(kernel)
  normalization <- match.arg(normalization)

  prep <- phymapnet_prepare_inputs(otu, tree, prune = prune_tree)
  otuA <- prep$otu
  dist <- prep$dist
  taxa <- prep$taxa

  norm_fun <- get_normalizer(normalization)
  Y <- norm_fun(otuA)

  p <- ncol(Y)
  K_neighbors <- as.integer(k * p)
  C <- kernel_from_dist(dist, alpha, kernel)

  precision <- phymapnet_precision(Y, nue = K_neighbors, epsilon1 = epsilon1, epsilon2 = epsilon2, C = C)
  sp <- sparse_quantile(precision, th_sparsity)

  adj <- sp$theta
  diag(adj) <- 0L

  list(
    precision_map = precision,
    adjacency = adj,
    threshold = sp$threshold,
    taxa = taxa,
    dist = dist,
    kernel_mat = C,
    params = list(alpha=alpha,k=k,epsilon1=epsilon1,epsilon2=epsilon2,kernel=kernel,th_sparsity=th_sparsity,normalization=normalization)
  )
}
