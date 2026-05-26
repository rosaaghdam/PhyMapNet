# Functions used by the bootstrap comparison for PhyMapNet estimation and
# binary-network evaluation. These preserve the calculations used in the
# original analysis script; they are kept local to make this workflow portable.

phymapnet <- function(Y, nue, epsilon1, epsilon2, C) {
  p <- ncol(Y)
  n <- nrow(Y)
  S1 <- cov(Y)
  V <- sqrt(diag(S1))
  V_diag <- diag(V)
  omega_hat <- V_diag %*% C %*% V_diag
  omega_hat <- omega_hat + diag(epsilon1, p)
  omega_star <- ((n / (n + nue)) * S1) + ((nue / (nue + n)) * omega_hat)
  nue_star <- nue + n
  IB <- as.matrix(nue_star * omega_star) + epsilon2
  B <- solve(IB)
  (nue_star - p - 1) * B
}

sparse_quantile <- function(Y, Isigma_star, quantile_level) {
  upper_tri_values <- abs(Isigma_star[upper.tri(Isigma_star)])
  ths <- quantile(upper_tri_values, quantile_level)
  Isigma_sparse <- ifelse(abs(Isigma_star) < ths, 0, Isigma_star)
  theta <- ifelse(Isigma_sparse != 0, 1, 0)
  list(Isigma_sparse = Isigma_sparse, theta = theta, threshold = ths)
}

evaluate <- function(number_of_genes, AD, TrueNet) {
  matches <- matrix(0, nrow(TrueNet), number_of_genes)
  for (i in seq_len(nrow(TrueNet))) {
    for (j in i:number_of_genes) {
      if (AD[i, j] == 1 && TrueNet[i, j] == 1) {
        matches[i, j] <- 1
      } else if (AD[i, j] == 0 && TrueNet[i, j] == 0) {
        matches[i, j] <- 0
      } else if (AD[i, j] == 1 && TrueNet[i, j] == 0) {
        matches[i, j] <- 2
      } else {
        matches[i, j] <- 3
      }
    }
  }

  matches_up <- matches[upper.tri(matches, diag = FALSE)]
  TP <- length(which(matches_up == 1))
  TN <- length(which(matches_up == 0))
  FP <- length(which(matches_up == 2))
  FN <- length(which(matches_up == 3))
  diagnostic_measures_undir <- cbind(TP, TN, FP, FN)

  Recall <- TP / (TP + FN)
  Precision <- TP / (TP + FP)
  Specificity <- TN / (TN + FP)
  Accuracy <- (TP + TN) / (TP + FP + TN + FN)
  F_score <- 2 * (Recall * Precision) / (Recall + Precision)
  accuracy_measures_undir <- cbind(Recall, Specificity, Precision, Accuracy, F_score)

  list(diagnostic_measures_undir, accuracy_measures_undir)
}
