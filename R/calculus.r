#' @title Statistical Inference with Persistence Landscapes
#' @description Compute test statistics and conduct null hypothesis tests for
#'   persistence data using persistence landscapes. See Section 3 of Bubenik
#'   (2015).
#'
#' @name statistical-inference
#' @include PersistenceLandscape.r
#' @inheritParams landscape
#' @inheritParams arithmetic-operations
#' @param supports List of support intervals for landscape levels.
#' @param r Non-negative number; the power of the coefficient \eqn{1/k} in the
#'   indicator linear form.
#' @param p Positive integer or infinity; the power used to compute an integral.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape').
#' @seealso PersistenceLandscape-methods
#' @example inst/examples/ex-inference.r
NULL

#' @rdname statistical-inference
#' @export
pl_integral <- function(pl, p = 1) {
  p <- ensure_p(p)
  pl$integral(p)
}

#' @rdname arithmetic-operations
#' @export
pl_distance <- function(pl1, pl2, p = 2) {
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  pl1$distance(pl2, p)
}

#' @rdname arithmetic-operations
#' @export
pl_dist <- function(pl_list, p = 2) {
  p <- ensure_p(p)
  PLdist(pl_list, p)
}

#' @rdname arithmetic-operations
#' @export
pl_norm <- function(pl, p = 2) {
  p <- ensure_p(p)
  pl$norm(p)
}

#' @rdname statistical-inference
#' @export
pl_indicator <- function(pl, supports, r = 0) {
  pl$indicator(supports, r = r)
}

#' @rdname statistical-inference
#' @export
pl_indicator_form <- function(pl, supports, r = 0, p = 1) {
  p <- ensure_p(p)
  pl$indicator_form(supports, r = r, p = p)
}
