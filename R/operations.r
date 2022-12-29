#' @title Arithmetic Operations on Persistence Landscapes
#' @description Calculate sums, scalar multiples, absolute values, means, inner
#'   products, extrema, moments, distances, norms, and products with indicator
#'   functions of persistent landscapes. These operations arise from the Hilbert
#'   space structure on persistence landscapes ().
#'
#' @name arithmetic-operations
#' @include PersistenceLandscape.r
#' @param pl,pl1,pl2 Persistent landscapes.
#' @param pl_list A list of persistent landscapes.
#' @param mult Double; a real-valued scale factor.
#' @param level Positive integer; the envelope of the persistence landscape (up
#'   to) whose moment to calculate.
#' @param p Positive integer or infinity; the power used to compute a norm or
#'   moment.
#' @param center Double; where to center the moment.
#' @return A persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape').
#' @seealso PersistenceLandscape-methods
#' @example inst/examples/ex-operations.r
NULL

#' @rdname arithmetic-operations
#' @export
pl_sum <- function(pl1, pl2) {
  # PLsum(pl1, pl2)
  pl1$add(pl2)
}

#' @rdname arithmetic-operations
#' @export
pl_scale <- function(pl, mult = 1) {
  # PLscale(mult, pl)
  pl$scale(mult)
}

#' @rdname arithmetic-operations
#' @export
pl_abs <- function(pl) {
  # PLabs(pl)
  pl$abs()
}

#' @rdname arithmetic-operations
#' @export
pl_mean <- function(pl_list) {
  PLaverage(pl_list)
}

#' @rdname arithmetic-operations
#' @export
pl_inner <- function(pl1, pl2) {
  # PLinner(pl1, pl2)
  pl1$inner(pl2)
}

#' @rdname arithmetic-operations
#' @export
pl_min <- function(pl, level = 1L) {
  pl$infimum(level)
}

#' @rdname arithmetic-operations
#' @export
pl_max <- function(pl, level = 1L) {
  pl$supremum(level)
}

#' @rdname arithmetic-operations
#' @export
pl_range <- function(pl, level = 1L) {
  c(pl$infimum(level), pl$supremum(level))
}

#' @rdname arithmetic-operations
#' @export
pl_vmin <- function(pl, level = 1L) {
  vapply(seq(level), function(l) pl$infimum(level = l), double())
}

#' @rdname arithmetic-operations
#' @export
pl_vmax <- function(pl, level = 1L) {
  vapply(seq(level), function(l) pl$supremum(level = l), double())
}

#' @rdname arithmetic-operations
#' @export
pl_vrange <- function(pl, level = 1L) {
  cbind(
    vapply(seq(level), function(l) pl$infimum(level = l), double()),
    vapply(seq(level), function(l) pl$supremum(level = l), double())
  )
}

#' @rdname arithmetic-operations
#' @export
pl_moment <- function(pl, p = 1L, center = 0, level = 1L) {
  p <- ensure_p(p)
  pl$moment(p, center, level)
}

#' @rdname arithmetic-operations
#' @export
pl_vmoment <- function(pl, p = 1L, center = 0, level = NULL) {
  p <- ensure_p(p)
  if (is.null(level)) level <- pl_num_envelopes(pl)
  vapply(
    seq(level),
    function(l) pl$moment(p = p, center = center, level = l),
    double()
  )
}

#' @rdname arithmetic-operations
#' @export
pl_distance <- function(pl1, pl2, p = 2) {
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  # PLdistance(pl1, pl2)
  pl1$distance(pl2, p)
}

#' @rdname arithmetic-operations
#' @export
pl_norm <- function(pl, p = 2) {
  p <- ensure_p(p)
  pl$norm(p)
}

#' @rdname arithmetic-operations
#' @export
pl_indicate <- function(pl, indicator) {
  pl$indicator_form(indicator)
}

# pre-process power
ensure_p <- function(p) {
  # only allow positive integer powers
  if (p < 1 || (p != Inf && p %% 1 != 0))
    stop("`p` must be a positive integer or infinity.")
  p
}
