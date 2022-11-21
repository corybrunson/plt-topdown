#' @title Print Method for Persistence Landscapes.
#' @description A [methods::show()] S4 method for persistence landscape objects,
#'   used for [base::print()].
#'
#' @name plot
#' @aliases show,Rcpp_PersistenceLandscape-method
#' @include PersistenceLandscape.r
#' @param object A [PersistenceLandscape] object.
#' @example inst/examples/ex-PersistenceLandscape.r
#' @export
setMethod(
  "show",
  "Rcpp_PersistenceLandscape",
  function(object) {
    
    # internal structure
    fmt <- str_internal(object)
    # number of envelopes
    envn <- num_levels(object)
    # abscissal support
    xran <- fmt_ran(supp_range(object))
    
    # send summary line to console
    cat(sprintf(
      "Persistence landscape (%s format) of %i envelopes over (%s,%s)",
      fmt, envn, xran[[1L]], xran[[2L]]
    ))
    
  }
)

fmt_ran <- function(x) {
  stopifnot(is.atomic(x) && is.double(x) && length(x) == 2L)
  xsmall <- max(0L, ceiling(-log(abs(diff(x)), 10)))
  format(x, digits = max(1L, xsmall), nsmall = xsmall)
}