#' @title Persistence Landscapes
#' @description Compute persistence landscapes from persistence data.
#'
#' @name landscape
#' @include PersistenceLandscape.r
#' @param pd Persistence data (or diagram), stored as a 2-column matrix or as a
#'   persistence diagram object with element `pairs` being a list of 2-column
#'   matrices (birth-death pairs).
#' @param degree Non-negative integer; if input is a persistence diagram object,
#'   then the dimension for which to compute a landscape. (For degree $d$, the
#'   $(d+1)$th matrix in the list will be selected.)
#' @param exact Set to `TRUE` for exact computation, `FALSE` (default) for
#'   discrete.
#' @param min_x,max_x Domain thresholds for discrete PL; if not specified, then
#'   taken to be the support of the PL.
#' @param by Domain grid diameter for discrete PL; if not specified, then set to
#'   the power of 10 that yields between 100 and 1000 intervals.
#' @param max_y Numeric; the threshold used to compute the persistence
#'   diagram (could be infinite).
#' @param pl A persistence landscape as returned by `landscape()`.
#' @return `landscape()` returns a persistence landscape (an object of S4 class
#'   'Rcpp_PersistenceLandscape'). Other functions return summary information
#'   about such an object.
#' @example inst/examples/ex-landscape-exact.r
#' @example inst/examples/ex-landscape-discrete.r
#' @export
landscape <- function(
    pd, degree = NULL,
    exact = FALSE,
    min_x = NULL, max_x = NULL, by = NULL, max_y = NULL
) {
  
  # birth-death pairs matrix `diagram` with upper bound `max_y`
  if (inherits(pd, "persistence")) {
    if (is.null(degree))
      stop("`landscape()` requires a homological degree (`degree = <int>`).")
    diagram <- pd$pairs[[degree + 1L]]
    if (! is.na(pd$threshold)) max_y <- pd$threshold
  } else if (is.atomic(pd)) {
    diagram <- pd
    stopifnot(ncol(diagram) >= 2L, is.numeric(diagram))
  }
  max_y <- max_y %||% max(diagram[, 2L])
  
  # content check
  if (is.null(diagram) || all(is.na(diagram))) {
    stop("`landscape()` requires non-empty, non-missing persistence data.")
  }
  
  # infer any missing parameters from the diagram
  min_x <- min_x %||% min(diagram)
  max_x <- max_x %||% max(diagram)
  if (! min_x < max_x) stop("Must have `min_x < max_x`.")
  # grid of between 100 and 1000 intervals of length a power of 10
  by <- by %||% 10 ^ (floor(log(max_x - min_x, 10)) - 2L)
  
  # construct persistence landscape
  new(PersistenceLandscape, diagram, exact, min_x, max_x, by, max_y)
}

# pl_str <- function(pl) {
#   stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
#   if (inherits(internal, "list")) {
#     "exact"
#   } else if (inherits(internal, "array")) {
#     "discrete"
#   } else {
#     stop(
#       "Internal structure of `",
#       deparse(substitute(object)),
#       "` cannot be determined."
#     )
#   }
# }

#' @rdname landscape
#' @export
pl_is_exact <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  pl$isExact()
}

#' @rdname landscape
#' @export
pl_str <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  # if (is.atomic(pl$getInternal())) "discrete" else "exact"
  if (pl$isExact()) "exact" else "discrete"
}

#' @rdname landscape
#' @export
pl_num_envelopes <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    pl_str(pl),
    exact = return(length(pl$getInternal())),
    discrete = return(dim(pl$getInternal())[[1L]])
  )
  # if (pl_str(pl) == "exact") {
  #   return(length(pl$getInternal()))
  # } else {
  #   return(dim(pl$getInternal())[[1L]])
  # }
}

#' @rdname landscape
#' @export
pl_limits <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    pl_str(pl),
    exact = {
      xvec <- Reduce(rbind, pl$getInternal())[, 1L, drop = TRUE]
      return(range(xvec))
    },
    discrete = {
      return(range(pl$getInternal()[1L, , 1L]))
    }
  )
}

#' @rdname landscape
#' @export
pl_support <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    pl_str(pl),
    exact = {
      xvec <- Reduce(rbind, pl$getInternal())[, 1L, drop = TRUE]
      return(range(xvec[! is.infinite(xvec)]))
    },
    discrete = {
      nz <- which(apply(pl$getInternal()[, , 2L], 2L, max) > 0)
      supp <- intersect(
        Reduce(union, list(nz - 1L, nz, nz + 1L)),
        seq(dim(pl$getInternal())[[2L]])
      )
      return(range(pl$getInternal()[1L, supp, 1L]))
    }
  )
}
