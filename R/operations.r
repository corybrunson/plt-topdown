#' @title Arithmetic and Statistical Operations on Persistence Landscapes
#' @description Calculate sums, scalar multiples, absolute values, means, inner
#'   products, extrema, and moments of persistent landscapes. These operations
#'   arise from the Hilbert space structure on persistence landscapes (Bubenik,
#'   2015).
#'
#' @details These functions are prefixed `pl_*()` to help users access them via
#'   tab-completion. Some take their names from the underlying S4 class methods
#'   and are only provided to enable composition via pipes: `add`, `scale`,
#'   `abs`, `inner`, `min` (`infimum`), `max` (`supremum`), and `moment`;
#'   `range` combines `min` and `max`. Others mimic classic R functions to
#'   handle lists of persistence landscapes: `sum`, `diff`, `mean`, `var`, and
#'   `sd`. Finally, some are vectorizations of the preceding: `vmin`, `vmax`,
#'   `vrange`, and `vmoment`.
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
#'   'Rcpp_PersistenceLandscape'), a real number, or a vector of real numbers.
#' @seealso PersistenceLandscape-methods
#' @example inst/examples/ex-operations.r
NULL

#' @rdname arithmetic-operations
#' @export
pl_add <- function(pl1, pl2) {
  pl1$add(pl2)
}

#' @rdname arithmetic-operations
#' @export
pl_sum <- function(pl_list) {
  PLsum(pl_list)
}

#' @rdname arithmetic-operations
#' @export
pl_scale <- function(pl, mult) {
  pl$scale(mult)
}

#' @rdname arithmetic-operations
#' @export
pl_abs <- function(pl) {
  pl$abs()
}

#' @rdname arithmetic-operations
#' @export
pl_diff <- function(pl_list) {
  PLdiff(pl_list)
}

#' @rdname arithmetic-operations
#' @export
pl_mean <- function(pl_list) {
  PLaverage(pl_list)
}

#' @rdname arithmetic-operations
#' @export
pl_var <- function(pl_list, p = 2) {
  if (length(pl_list) <= 1L) return(NA_real_)
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  PLvar(pl_list, p)
}

#' @rdname arithmetic-operations
#' @export
pl_sd <- function(pl_list, p = 2) {
  if (length(pl_list) <= 1L) return(NA_real_)
  p <- ensure_p(p)
  if (p == Inf) p <- 0
  PLsd(pl_list, p)
}

#' @rdname arithmetic-operations
#' @export
pl_inner <- function(pl1, pl2) {
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
  vapply(seq(level), function(l) pl$infimum(l), double())
}

#' @rdname arithmetic-operations
#' @export
pl_vmax <- function(pl, level = 1L) {
  vapply(seq(level), function(l) pl$supremum(l), double())
}

#' @rdname arithmetic-operations
#' @export
pl_vrange <- function(pl, level = 1L) {
  cbind(
    vapply(seq(level), function(l) pl$infimum(l), double()),
    vapply(seq(level), function(l) pl$supremum(l), double())
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
    function(l) pl$moment(p, center, l),
    double()
  )
}
