#' @title Persistence Landscapes
#' @description Compute persistence landscapes from persistence data.
#'
#' @name landscape
#' @include PersistenceLandscape.r
#' @param pd Persistence data (or diagram), stored as a 2-column matrix or as a
#'   persistence diagram object with element `pairs` being a list of 2-column
#'   matrices (birth-death pairs).
#' @param degree Non-negative integer; if input is a persistence diagram
#'   object, then the dimension for which to compute a landscape. (For degree
#'   $d$, the $(d+1)$th matrix in the list will be selected.)
#' @param exact Set to `TRUE` for exact computation, `FALSE` (default) for
#'   discrete.
#' @param dx Domain grid diameter for discrete PL (should be more intuitively
#'   named).
#' @param min_x,max_x Domain thresholds for discrete PL (should be more
#'   intuitively named).
#' @param threshold Numeric; the threshold used to compute the persistence
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
    exact = FALSE, dx = 0.1, min_x = 0, max_x = 10, threshold = NULL
) {
  
  # birth-death pairs matrix `diagram` and threshold `max_y`
  if (inherits(pd, "persistence")) {
    if (is.null(degree))
      stop("`landscape()` requires a homological degree (`degree = <int>`).")
    diagram <- pd$pairs[[degree + 1L]]
    if (! is.na(pd$threshold)) threshold <- pd$threshold
    max_y <- threshold %||% max(diagram[, 2L])
  } else if (is.atomic(pd)) {
    diagram <- pd
    stopifnot(ncol(diagram) >= 2L, is.numeric(diagram))
    max_y <- threshold %||% max(diagram[, 2L])
  }
  
  # content check
  if (is.null(diagram) || all(is.na(diagram))) {
    stop("`landscape()` requires non-empty, non-missing persistence data.")
  }
  
  # construct persistence landscape
  new(PersistenceLandscape, diagram, exact, min_x, max_x, dx, max_y)
}

# str_internal <- function(pl) {
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
str_internal <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  if (is.atomic(pl$getInternal())) "discrete" else "exact"
}

#' @rdname landscape
#' @export
num_levels <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    str_internal(pl),
    exact = return(length(pl$getInternal())),
    discrete = return(dim(pl$getInternal())[[1L]])
  )
  # if (str_internal(pl) == "exact") {
  #   return(length(pl$getInternal()))
  # } else {
  #   return(dim(pl$getInternal())[[1L]])
  # }
}

#' @rdname landscape
#' @export
str_range <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    str_internal(pl),
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
supp_range <- function(pl) {
  stopifnot(inherits(pl, "Rcpp_PersistenceLandscape"))
  switch (
    str_internal(pl),
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
