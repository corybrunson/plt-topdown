#include "PersistenceLandscapeInterface.h"
#include <Rcpp.h>

using namespace Rcpp;

RCPP_EXPOSED_CLASS(PersistenceLandscapeInterface)

RCPP_MODULE(Landscape) {

  // NB: Expose C++ class 'PersistenceLandscapeInterface' as R class
  // 'PersistenceLandscape'
  class_<PersistenceLandscapeInterface>("PersistenceLandscape")
  
  .constructor<NumericMatrix, bool, double, double, double, double>()
  
  // TODO: Implement validator. -JCB

  // TODO: Delete `isExact` method if `exact` field can be read-only exposed.
  // -JCB
  // .field_readonly("exact",
  //                 &PersistenceLandscapeInterface::exact,
  //                 "Representation of the underlying PL.")
  
  .method("isExact",
          &PersistenceLandscapeInterface::isExact,
          "Queries whether the underlying PL representation is exact.")
  .method("xMin",
          &PersistenceLandscapeInterface::xMin,
          "Returns the infimum (left endpoint) of the PL support.")
  .method("xMax",
          &PersistenceLandscapeInterface::xMax,
          "Returns the supremum (right endpoint) of the PL support.")
  .method("xBy",
          &PersistenceLandscapeInterface::xBy,
          "Returns the resolution of the discrete representation.")
  .method("getInternal",
          &PersistenceLandscapeInterface::getInternal,
          "Returns the internal tensor representation of the PL.")
  .method("toExact",
          &PersistenceLandscapeInterface::toExact,
          "Returns the PL in exact representation.")
  .method("toDiscrete",
          &PersistenceLandscapeInterface::toDiscrete,
          "Returns the PL in discrete representation.")
  .method("expand",
          &PersistenceLandscapeInterface::expand,
          "Expands the limits of this PL.")
  .method("add",
          &PersistenceLandscapeInterface::add,
          "Adds this PL to another.")
  .method("scale",
          &PersistenceLandscapeInterface::scale,
          "Multiplies this PL by a scalar.")
  .method("abs",
          &PersistenceLandscapeInterface::abs,
          "Takes the absolute value of this PL.")
  .method("inner",
          &PersistenceLandscapeInterface::inner,
          "Takes the inner product of this PL with another.")
  .method("minimum",
          &PersistenceLandscapeInterface::minimum,
          "Finds the minimum value of one level of this PL.")
  .method("maximum",
          &PersistenceLandscapeInterface::maximum,
          "Finds the maximum value of one level of this PL.")
  .method("moment",
          &PersistenceLandscapeInterface::moment,
          "Computes the n^th moment of one level of this PL.")
  .method("integral",
          &PersistenceLandscapeInterface::integral,
          "Computes the integral of this PL.")
  .method("distance",
          &PersistenceLandscapeInterface::distance,
          "Takes the p-distance between this PL and another.")
  .method("norm",
          &PersistenceLandscapeInterface::norm,
          "Computes the p-norm of this PL.")
  .method("indicator",
          &PersistenceLandscapeInterface::indicator,
          "Multiplies this PL by a level-indexed set of indicator functions.")
  .method("indicator_form",
          &PersistenceLandscapeInterface::indicator_form,
          "Computes the integral of the productof this PL with an indicator.")
  ;

  Rcpp::function("PLsum", &PLsum);
  Rcpp::function("PLdiff", &PLdiff);
  Rcpp::function("PLmean", &PLmean);
  Rcpp::function("PLvar", &PLvar);
  Rcpp::function("PLsd", &PLsd);
}
