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

# load modules (can be done anywhere, will activate at load time)
Rcpp::loadModule("class_PersistenceLandscape", loadNow = TRUE)

Rcpp::exposeClass(
  "PersistenceLandscape",
  constructors = list(
    c("NumericMatrix", "bool", "double", "double", "double", "double", "double")
  ),
  fields = c(),
  methods = c(
    "isExact", "xMin", "xMax", "xBy",
    "getInternal", "toDiscrete", "toExact", "discretize",
    "expand", "contract"
    # "abs",
    # "minimum", "maximum", "moment", "computeNormOfLandscape",
    #"integral", "distance", "norm", "indicator", "indicator_form"
  ),
  header = '#include "PersistenceLandscape.h"',
  CppClass = "PersistenceLandscape"
)
