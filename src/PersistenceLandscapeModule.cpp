#include <Rcpp.h>
using namespace Rcpp ;
#include "PersistenceLandscape.h"

RCPP_MODULE(class_PersistenceLandscape) {


    class_<PersistenceLandscape>("PersistenceLandscape")

    .constructor<NumericMatrix,bool,double,double,double,double,double>()


    .method("isExact", &PersistenceLandscape::isExact)
    .method("xMin", &PersistenceLandscape::xMin)
    .method("xMax", &PersistenceLandscape::xMax)
    .method("xBy", &PersistenceLandscape::xBy)
    .method("getInternal", &PersistenceLandscape::getInternal)
    .method("toDiscrete", &PersistenceLandscape::toDiscrete)
    .method("toExact", &PersistenceLandscape::toExact)
    .method("discretize", &PersistenceLandscape::discretize)
    .method("expand", &PersistenceLandscape::expand)
    .method("contract", &PersistenceLandscape::contract)
    ;
}
