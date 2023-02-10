// #include "PersistenceLandscape.h"
// #include <Rcpp.h>
// 
// using namespace Rcpp;
// 
// RCPP_EXPOSED_CLASS(PersistenceLandscape)
// 
// RCPP_MODULE(Landscape) {
//   
//   // NB: Expose C++ class 'PersistenceLandscape' as R class
//   // 'PersistenceLandscape'
//   class_<PersistenceLandscape>("PersistenceLandscape")
//   
//   // Match specifications in `PersistenceLandscape.h`.
//   // TODO: Implement a validator. -JCB
//   .constructor<NumericMatrix, bool, double, double, double, double, double>(""
//   "Creates a PL from a PD, in the form of a 2-column numeric matrix")
//   
//   // TODO: Delete `isExact` method if `exact` field can be read-only exposed.
//   // -JCB
//   // .field_readonly("exact",
//   // &PersistenceLandscape::exact,
//   // "Representation of the underlying PL.")
//   
//   .method("isExact",
//   &PersistenceLandscape::isExact,
//   "Queries whether the underlying PL representation is exact")
//   .method("xMin",
//   &PersistenceLandscape::xMin,
//   "Returns the infimum (left endpoint) of the PL support")
//   .method("xMax",
//   &PersistenceLandscape::xMax,
//   "Returns the supremum (right endpoint) of the PL support")
//   .method("xBy",
//   &PersistenceLandscape::xBy,
//   "Returns the resolution of the discrete representation")
//   .method("getInternal",
//   &PersistenceLandscape::getInternal,
//   "Returns the internal tensor representation of the PL")
//   .method("toExact",
//   &PersistenceLandscape::toExact,
//   "Returns the PL in exact representation")
//   .method("toDiscrete",
//   &PersistenceLandscape::toDiscrete,
//   "Returns the PL in discrete representation")
//   .method("discretize",
//   &PersistenceLandscape::discretize,
//   "Casts an exact PL to a discrete one")
//   .method("expand",
//   &PersistenceLandscape::expand,
//   "Expands the limits of this PL")
//   .method("contract",
//   &PersistenceLandscape::contract,
//   "Contracts the limits of this PL")
//   // .method("add",
//   // &PersistenceLandscape::add,
//   // "Adds this PL to another")
//   // .method("scale",
//   // &PersistenceLandscape::scale,
//   // "Multiplies this PL by a scalar")
//   .method("abs",
//   &PersistenceLandscape::abs,
//   "Takes the absolute value of this PL")
//   // .method("inner",
//   // &PersistenceLandscape::inner,
//   // "Takes the inner product of this PL with another")
//   .method("minimum",
//   &PersistenceLandscape::minimum,
//   "Finds the minimum value of one level of this PL")
//   .method("maximum",
//   &PersistenceLandscape::maximum,
//   "Finds the maximum value of one level of this PL")
//   .method("moment",
//   &PersistenceLandscape::moment,
//   "Computes the n^th moment of one level of this PL")
//   // .method("computeIntegralOfLandscape",
//   // &PersistenceLandscape::computeIntegralOfLandscape,
//   // "Computes the integral of this PL")
//   // .method("distance",
//   // &PersistenceLandscape::distance,
//   // "Takes the p-distance between this PL and another")
//   .method("computeNormOfLandscape",
//   &PersistenceLandscape::computeNormOfLandscape,
//   "Computes the p-norm of this PL")
//   // .method("multiplyByIndicatorFunction",
//   // &PersistenceLandscape::multiplyByIndicatorFunction,
//   // "Multiplies this PL by a level-indexed set of indicator functions")
//   // .method("computeIntegralOfLandscapeMultipliedByIndicatorFunction",
//   // &PersistenceLandscape::computeIntegralOfLandscapeMultipliedByIndicatorFunction,
//   // "Computes the integral of the productof this PL with an indicator")
//   ;
//   
//   Rcpp::function("PLsum", &PLsum);
//   // Rcpp::function("PLdiff", &PLdiff);
//   Rcpp::function("PLmean", &PLmean);
//   Rcpp::function("PLdist", &PLdist);
//   Rcpp::function("PLvar", &PLvar);
//   Rcpp::function("PLsd", &PLsd);
// }
