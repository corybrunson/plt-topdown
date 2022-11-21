#' @title Hilbert Space Operations on Persistence Landscapes
#' @description Calculate sums, scalar multiples, means, and inner products of
#'   persistent landscapes.
#'
#' @name landscape-operations
#' @aliases hilbert-space-operations
#' @include PersistenceLandscape.r
#' @param pl,pl1,pl2 Persistent landscapes.
#' @param mult Double; a real-valued scale factor.
#' @param pl_list A list of persistent landscapes.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape').
#' @seealso PersistenceLandscape-methods
#' @example inst/examples/ex-landscape-operations.r
NULL

#' @rdname landscape-operations
#' @export
pl_sum <- function(pl1, pl2) {
  # PLsum(pl1, pl2)
  pl1$add(pl2)
}

#' @rdname landscape-operations
#' @export
pl_scale <- function(pl, mult = 1) {
  # PLscale(mult, pl)
  pl$scale(mult)
}

#' @rdname landscape-operations
#' @export
pl_mean <- function(pl_list) {
  PLaverage(pl_list)
}

#' @rdname landscape-operations
#' @export
pl_inner_product <- function(pl1, pl2) {
  # PLinner(pl1, pl2)
  pl1$inner(pl2)
}