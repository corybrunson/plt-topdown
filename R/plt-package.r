#' @title 'plt' package
#' @description An Rcpp interface to Dłotko's Persistence Landscapes Toolbox
#'
#' @docType package
#' @name plt-package
#' @aliases plt
#' @import methods Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib plt, .registration = TRUE
#' @details This package provides an R interface to the Persistence Landscapes
#'   Toolbox by Paweł Dłotko and Peter Bubenik. It is adapted from the
#'   'tdatools' package prepared by Jose Bouza.
"_PACKAGE"

# load {Rcpp} module(s) (can be done anywhere, will activate at load time)
Rcpp::loadModule("PersistenceLandscape", loadNow = TRUE)

Rcpp::exposeClass(
  "PersistenceLandscape",
  constructors = list(
    c("NumericMatrix", "bool", "double", "double", "double", "double")
  ),
  fields = character(0L),
  methods = c(
    "isExact", "xMin", "xMax", "xBy",
    "toExact", "toDiscrete", "getInternal",
    "add", "scale", "inner"
  ),
  header = '#include "PersistenceLandscapeInterface.h"',
  CppClass = "PersistenceLandscapeInterface"
)

# debugging notes
# error: static_assert failed due to requirement '!sizeof(PersistenceLandscapeInterface)' "cannot convert type to SEXP"
# https://github.com/RcppCore/Rcpp/blob/9633169d6a7fda656e7752ab32acbad51554358f/inst/include/Rcpp/internal/wrap.h#L507
