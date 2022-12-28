// #include "pl.h"
// #include <Rcpp.h>
// 
// using namespace Rcpp;
// 
// RCPP_EXPOSED_CLASS(PersistenceLandscapeInterface)
// 
// RCPP_MODULE(Landscape) {
// 
//   // NB: Expose C++ class 'PersistenceLandscapeInterface' as R class
//   // 'PersistenceLandscape'
//   class_<PersistenceLandscapeInterface>("PersistenceLandscape")
//   
//   .constructor<NumericMatrix, bool, double, double, double, double>()
//   
//   // TODO: Implement validator. -JCB
// 
//   // TODO: Delete `isExact` method if `exact` field can be read-only exposed.
//   // -JCB
//   // .field_readonly("exact",
//   //                 &PersistenceLandscapeInterface::exact,
//   //                 "Representation of the underlying PL.")
//   
//   .method("isExact",
//           &PersistenceLandscapeInterface::isExact,
//           "Queries whether the underlying PL representation is exact.")
//   .method("getMin",
//           &PersistenceLandscapeInterface::getMin,
//           "Returns the infimum of the PL support.")
//   .method("getMax",
//           &PersistenceLandscapeInterface::getMax,
//           "Returns the supremum of the PL support.")
//   .method("getdx",
//           &PersistenceLandscapeInterface::getdx,
//           "Returns the resolution of the discrete representation.")
//   .method("getExact",
//           &PersistenceLandscapeInterface::getExact,
//           "Returns the PL in exact representation.")
//   .method("getDiscrete",
//           &PersistenceLandscapeInterface::getDiscrete,
//           "Returns the PL in discrete representation.")
//   .method("getInternal",
//           &PersistenceLandscapeInterface::getInternal,
//           "Returns the internal tensor representation of the PL.")
//   .method("add",
//           &PersistenceLandscapeInterface::add,
//           "Adds this PL to another.")
//   .method("scale",
//           &PersistenceLandscapeInterface::scale,
//           "Multiplies this PL by a scalar.")
//   .method("inner",
//           &PersistenceLandscapeInterface::inner,
//           "Takes the inner product of this PL with another.")
//   .method("distance",
//           &PersistenceLandscapeInterface::distance,
//           "Takes the p-norm distance between this PL and another.")
//   ;
// 
//   // TODO: Decide whether to use these or R wrappers in
//   // 'landscape-operations.r'. -JCB
//   Rcpp::function("PLsum", &PLsum);
//   Rcpp::function("PLscale", &PLscale);
//   Rcpp::function("PLinner", &PLinner);
//   Rcpp::function("PLdistance", &PLdistance);
//   Rcpp::function("PLaverage", &PLaverage);
// }
