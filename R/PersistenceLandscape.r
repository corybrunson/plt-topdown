#' @title Exported C++ Class 'PersistenceLandscape'
#' @description Export, and create and manipulate objects of, the
#'   'PersistenceLandscape' C++ class.
#' @details The C++ class 'PersistenceLandscape' is exposed as the S4 class
#'   'Rcpp_PersistenceLandscape' via the `RCPP_MODULE()` macro provided by
#'   **[Rcpp][Rcpp::Rcpp-package]**. See
#'   <https://github.com/r-pkg-examples/rcpp-modules-student> for an
#'   introduction. New objects should be created from persistence data
#'   (diagrams) using [landscape()].
#'
#' @docType class
#' @name PersistenceLandscape
#' @aliases Rcpp_PersistenceLandscape
#' @aliases PersistenceLandscape-class
#' @aliases Rcpp_PersistenceLandscape-class
#' @include plt-package.r
#' @inheritParams base::as.vector
#' @inheritParams base::as.data.frame
#' @param exact Whether to export the exact or a discrete (default)
#'   representation; if `TRUE` but `x` has a discrete representation, then
#'   ignored with a warning.
#' @example inst/examples/ex-PersistenceLandscape.r
#' @export PersistenceLandscape
#' @export Rcpp_PersistenceLandscape
#' @exportClass PersistenceLandscape
#' @exportClass Rcpp_PersistenceLandscape

PersistenceLandscape <- setClass("PersistenceLandscape")
Rcpp_PersistenceLandscape <- setClass("Rcpp_PersistenceLandscape")

# register S4 class for S3 inheritance
setOldClass("Rcpp_PersistenceLandscape")

#' @rdname PersistenceLandscape
#' @export
as.vector.Rcpp_PersistenceLandscape <- function(x, mode = "any") {
  # get discrete representation for consistent abscissa values
  internal <- x$getDiscrete()
  
  # export concatenation of envelopes
  as.vector(t(internal[, , 2L]), mode = mode)
}

#' @rdname PersistenceLandscape
#' @export
as.data.frame.Rcpp_PersistenceLandscape <- function(
    x, row.names = NULL, optional = FALSE, exact = FALSE, ...
) {
  if (exact) {
    internal <- try(x$getExact(), silent = TRUE)
    if (inherits(internal, "try-error")) {
      warning("Cannot recover exact PL data from a discrete PL object.")
      # recurse
      return(as.data.frame.Rcpp_PersistenceLandscape(
        x, row.names = row.names, optional = optional, exact = FALSE
      ))
    }
    envelopes <- vapply(internal, nrow, 0L)
    envelope <- rep(seq_along(internal), envelopes)
    yfy <- Reduce(rbind, internal)
    df <- data.frame(envelope, yfy)
    names(df) <- c("envelope", "x", "fx")
    if (! is.null(row.names)) {
      rownames(df) <- row.names
    } else if (! optional) {
      id <- unlist(lapply(envelopes, seq))
      rownames(df) <- paste(envelope, id, sep = ".")
    }
    df
  } else {
    internal <- x$getDiscrete()
    envelope <- rep(seq(dim(internal)[[2L]]), each = dim(internal)[[1L]])
    y <- as.vector(t(internal[, , 1L]))
    fy <- as.vector(t(internal[, , 2L]))
    df <- data.frame(envelope = envelope, x = y, fx = fy)
    if (! is.null(row.names) && length(row.names) == nrow(df)) {
      # adapted from `as.data.frame.matrix()`
      .rowNamesDF(df, make.names = TRUE) <- row.names
    } else if (! optional) {
      # adapted from `as.data.frame.matrix()`
      attr(df, "row.names") <- paste(envelope, y, sep = ".")
    }
    df
  }
}
