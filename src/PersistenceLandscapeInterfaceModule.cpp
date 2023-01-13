#include <Rcpp.h>
using namespace Rcpp ;
#include "PersistenceLandscapeInterface.h"

RCPP_MODULE(class_PersistenceLandscape) {


    class_<PersistenceLandscapeInterface>("PersistenceLandscape")

    .constructor<NumericMatrix,bool,double,double,double,double>()

    .field("exact", &PersistenceLandscapeInterface::exact)

    .method("isExact", &PersistenceLandscapeInterface::isExact)
    .method("xMin", &PersistenceLandscapeInterface::xMin)
    .method("xMax", &PersistenceLandscapeInterface::xMax)
    .method("xBy", &PersistenceLandscapeInterface::xBy)
    .method("toExact", &PersistenceLandscapeInterface::toExact)
    .method("toDiscrete", &PersistenceLandscapeInterface::toDiscrete)
    .method("getInternal", &PersistenceLandscapeInterface::getInternal)
    .method("add", &PersistenceLandscapeInterface::add)
    .method("scale", &PersistenceLandscapeInterface::scale)
    .method("inner", &PersistenceLandscapeInterface::inner)
    ;
}
