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
#' @inheritParams base::summary
#' @inheritParams base::print
#' @param exact Whether to export the exact or a discrete (default)
#'   representation; if `TRUE` but `x` has a discrete representation, then
#'   ignored with a warning.
#' @param ... Additional arguments; currently ignored.
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
summary.Rcpp_PersistenceLandscape <- function(object, ...) {
  res <- list(
    str = pl_str(object),
    n.levels = pl_num_levels(object),
    limits = pl_limits(object),
    resolution = if (object$isExact()) NA_real_ else object$xBy(),
    support = pl_support(object),
    range = c(pl_min(object), pl_max(object)),
    magnitude = object %*% object,
    integral = pl_integral(object)
  )
  class(res) <- "summary.Rcpp_PersistenceLandscape"
  res
}

#' @rdname PersistenceLandscape
#' @export
print.summary.Rcpp_PersistenceLandscape <- function(
    x, digits = max(1L, getOption("digits") - 3L), ...
) {
  fmt <- function(s) {
    if (is.na(s) || is.infinite(s)) return(format(s))
    format(round(s, max(0, digits - log10(abs(s)))))
  }
  cat("Internal representation: ", x$str, "\n")
  cat("Number of levels: ", x$n.levels, "\n")
  cat(
    "Representation limits: (",
    fmt(x$limits[[1L]]), ",", fmt(x$limits[[2L]]), ")",
    if (x$str == "discrete")
      paste0(" at resolution ", format(round(x$resolution))),
    "\n"
  )
  cat("Landscape range: (",
      fmt(x$range[[1L]]), ",", fmt(x$range[[2L]]), ")", "\n")
  cat("Magnitude: ", fmt(x$magnitude), "\n")
  cat("Integral:  ", fmt(x$integral), "\n")
}

#' @rdname PersistenceLandscape
#' @export
as.vector.Rcpp_PersistenceLandscape <- function(x, mode = "any") {
  # get discrete representation for consistent abscissa values
  internal <- x$toDiscrete()
  
  # export concatenation of levels
  as.vector(t(internal[, , 2L]), mode = mode)
}

#' @rdname PersistenceLandscape
#' @export
as.data.frame.Rcpp_PersistenceLandscape <- function(
    x, row.names = NULL, optional = FALSE, exact = FALSE, ...
) {
  if (exact) {
    internal <- try(x$toExact(), silent = TRUE)
    if (inherits(internal, "try-error")) {
      warning("Cannot recover exact PL data from a discrete PL object.")
      # recurse
      return(as.data.frame.Rcpp_PersistenceLandscape(
        x, row.names = row.names, optional = optional, exact = FALSE
      ))
    }
    levs <- vapply(internal, nrow, 0L)
    lev <- rep(seq_along(internal), levs)
    yfy <- Reduce(rbind, internal)
    df <- data.frame(lev, yfy)
    names(df) <- c("level", "x", "fx")
    if (! is.null(row.names)) {
      rownames(df) <- row.names
    } else if (! optional) {
      id <- unlist(lapply(levs, seq))
      rownames(df) <- paste(lev, id, sep = ".")
    }
    df
  } else {
    internal <- x$toDiscrete()
    lev <- rep(seq(dim(internal)[[2L]]), each = dim(internal)[[1L]])
    y <- as.vector(t(internal[, , 1L]))
    fy <- as.vector(t(internal[, , 2L]))
    df <- data.frame(level = lev, x = y, fx = fy)
    if (! is.null(row.names) && length(row.names) == nrow(df)) {
      # adapted from `as.data.frame.matrix()`
      .rowNamesDF(df, make.names = TRUE) <- row.names
    } else if (! optional) {
      # adapted from `as.data.frame.matrix()`
      attr(df, "row.names") <- paste(lev, y, sep = ".")
    }
    df
  }
}
