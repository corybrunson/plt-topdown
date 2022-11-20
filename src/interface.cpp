#include "pl.h"
#include <Rcpp.h>

using namespace Rcpp;
RCPP_EXPOSED_CLASS(PersistenceLandscapeInterface)
RCPP_MODULE(Landscape) {

  class_<PersistenceLandscapeInterface>("PersistenceLandscape")
      .constructor<NumericMatrix, bool, double, double, double, double>()

      .method("isExact",
              &PersistenceLandscapeInterface::isExact,
              "Queries whether the underlying PL representation is exact.")
      .method("getMin",
              &PersistenceLandscapeInterface::getMin,
              "Returns the infimum of the PL support.")
      .method("getMax",
              &PersistenceLandscapeInterface::getMax,
              "Returns the supremum of the PL support.")
      .method("getdx",
              &PersistenceLandscapeInterface::getdx,
              "Returns the resolution of the discrete representation.")
      .method("getExact",
              &PersistenceLandscapeInterface::getPersistenceLandscapeExact,
              "Returns the PL in exact representation.")
      .method("getDiscrete",
              &PersistenceLandscapeInterface::getPersistenceLandscapeDiscrete,
              "Returns the PL in discrete representation.")
      .method("getInternal", &PersistenceLandscapeInterface::getInternal,
              "Returns the internal tensor representation of the PL.")
      .method("add", &PersistenceLandscapeInterface::sum,
              "Adds this PL to another.")
      .method("scale", &PersistenceLandscapeInterface::scale,
              "Multiplies this PL by a scalar.")
      .method("inner", &PersistenceLandscapeInterface::inner,
              "Takes the inner product of this PL with another.");

  Rcpp::function("PLsum", &PLsum);
  Rcpp::function("PLscale", &PLscale);
  Rcpp::function("PLinner", &PLinner);
  Rcpp::function("PLaverage", &PLaverage);
}
