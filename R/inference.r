#' @title Statistical Inference with Persistence Landscapes
#' @description Compute test statistics and conduct null hypothesis tests for
#'   persistence data using persistence landscapes. See Section 3 of Bubenik
#'   (2015).
#'
#' @name statistical-inference
#' @include PersistenceLandscape.r
#' @param pl A persistent landscape.
#' @param r Non-negative number; the power of the coefficient \eqn{1/k} in the
#'   indicator linear form.
#' @param p Positive integer or infinity; the power used to compute an integral.
#' @param center Double; where to center the moment.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape').
#' @seealso PersistenceLandscape-methods
#' @example inst/examples/ex-inference.r
NULL

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
