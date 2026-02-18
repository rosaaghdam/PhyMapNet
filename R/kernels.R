#' Build a phylogenetic kernel from a distance matrix
#'
#' @param dist A taxa x taxa distance matrix.
#' @param alpha Positive bandwidth parameter.
#' @param kernel Either "gaussian" or "laplacian".
#' @return A taxa x taxa kernel matrix.
#' @keywords internal
kernel_from_dist <- function(dist, alpha, kernel = c("gaussian","laplacian")) {
  kernel <- match.arg(kernel)
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("alpha must be a single positive number.", call. = FALSE)
  }
  if (kernel == "gaussian") {
    exp(-(dist^2) / (2 * alpha^2))
  } else {
    exp(-dist / alpha)
  }
}
