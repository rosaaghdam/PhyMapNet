#' @keywords internal
norm_log1p <- function(x) log(x + 1)

#' @keywords internal
norm_gmpr <- function(x) {
  sf <- GUniFrac::GMPR(as.data.frame(x))  # rows=samples
  sweep(x, 1, sf, "/")
}

#' @keywords internal
norm_clr <- function(x, pseudo = 0.01) {
  y <- x + pseudo
  as.matrix(compositions::clr(y))
}

#' @keywords internal
norm_tss <- function(x, scale = 1) {
  x <- as.matrix(x)
  rs <- rowSums(x)
  rs_safe <- ifelse(rs == 0, 1, rs)
  y <- sweep(x, 1, rs_safe, "/")
  if (!is.null(scale) && is.finite(scale) && scale != 1) y <- y * scale
  y
}

#' @keywords internal
get_normalizer <- function(normalization = c("log","gmpr","clr","tss")) {
  normalization <- match.arg(normalization)
  switch(
    normalization,
    log  = norm_log1p,
    gmpr = norm_gmpr,
    clr  = norm_clr,
    tss  = norm_tss
  )
}
