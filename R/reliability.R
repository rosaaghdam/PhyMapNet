#' Edge reliability via hyperparameter ensemble
#'
#' Runs an ensemble over (alpha, k, epsilon1, epsilon2, kernel, normalization) and
#' returns edge reliability as selection frequency under fixed sparsification threshold.
#'
#' @param otu samples x taxa matrix.
#' @param tree phylo tree.
#' @param th_fixed fixed quantile threshold for sparsification across all models (e.g., 0.95).
#' @param alpha_range numeric vector.
#' @param k_range integer vector.
#' @param epsilon1_range numeric vector.
#' @param epsilon2_range numeric vector.
#' @param kernels character vector: "gaussian" and/or "laplacian".
#' @param normalizations character vector: subset of c("log","gmpr","clr","tss").
#' @param consensus_cut reliability cutoff for binary consensus (default 0.5).
#' @param prune_tree prune tree tips not in OTU.
#' @param progress print progress every `progress_every` models.
#' @param progress_every integer.
#' @return A list with rel_mat, consensus_mat, edge_list, N_models, grid.
#' @export
phymapnet_reliability <- function(
  otu, tree,
  th_fixed = 0.95,
  alpha_range = seq(0.01, 0.12, by = 0.01),
  k_range = 2:10,
  epsilon1_range = seq(0, 1, by = 0.1),
  epsilon2_range = seq(0, 1, by = 0.1),
  kernels = c("gaussian"),
  normalizations = c("log","gmpr","clr","tss"),
  consensus_cut = 0.5,
  prune_tree = TRUE,
  progress = TRUE,
  progress_every = 500
) {
  kernels <- intersect(kernels, c("gaussian","laplacian"))
  if (length(kernels) == 0) stop("kernels must include 'gaussian' and/or 'laplacian'.", call. = FALSE)

  prep <- phymapnet_prepare_inputs(otu, tree, prune = prune_tree)
  otuA <- prep$otu
  dist <- prep$dist
  taxa <- prep$taxa
  p <- length(taxa)

  count_mat <- matrix(0L, nrow = p, ncol = p, dimnames = list(taxa, taxa))

  N_models <- length(normalizations) * length(alpha_range) * length(k_range) *
    length(epsilon1_range) * length(epsilon2_range) * length(kernels)

  model_idx <- 0L
  t0 <- proc.time()

  for (norm_nm in normalizations) {
    norm_fun <- get_normalizer(norm_nm)
    Y <- norm_fun(otuA)

    for (alpha in alpha_range) {
      KER <- list()
      if ("gaussian" %in% kernels)  KER$gaussian  <- exp(-(dist^2) / (2 * alpha^2))
      if ("laplacian" %in% kernels) KER$laplacian <- exp(-dist / alpha)

      for (k in k_range) {
        K_neighbors <- as.integer(k * p)

        for (epsilon1 in epsilon1_range) {
          for (epsilon2 in epsilon2_range) {
            for (ker_nm in names(KER)) {
              C <- KER[[ker_nm]]

              precision <- phymapnet_precision(Y, nue = K_neighbors, epsilon1 = epsilon1, epsilon2 = epsilon2, C = C)
              sp <- sparse_quantile(precision, th_fixed)
              B <- sp$theta
              diag(B) <- 0L
              B[lower.tri(B)] <- 0L

              count_mat <- count_mat + (B > 0)
              model_idx <- model_idx + 1L

              if (progress && (model_idx %% progress_every == 0L)) {
                elapsed <- (proc.time() - t0)[3]
                message(sprintf("%d / %d models (%.1f%%)   elapsed: %.1fs",
                                model_idx, N_models, 100*model_idx/N_models, elapsed))
              }
            }
          }
        }
      }
    }
  }

  rel_upper <- count_mat / N_models
  rel_mat <- rel_upper + t(rel_upper)
  diag(rel_mat) <- 0
  rel_mat[rel_mat > 1] <- 1

  consensus_mat <- (rel_mat >= consensus_cut) * 1L
  diag(consensus_mat) <- 0L

  ut <- upper.tri(rel_upper, diag = FALSE)
  edge_list <- data.frame(
    node_i = taxa[row(rel_upper)[ut]],
    node_j = taxa[col(rel_upper)[ut]],
    reliability = as.numeric(rel_upper[ut]),
    stringsAsFactors = FALSE
  )
  edge_list <- edge_list[order(-edge_list$reliability), , drop = FALSE]

  list(
    rel_mat = rel_mat,
    consensus_mat = consensus_mat,
    edge_list = edge_list,
    N_models = N_models,
    grid = list(
      alpha_range=alpha_range, k_range=k_range,
      epsilon1_range=epsilon1_range, epsilon2_range=epsilon2_range,
      kernels=kernels, normalizations=normalizations, th_fixed=th_fixed,
      consensus_cut=consensus_cut
    )
  )
}
