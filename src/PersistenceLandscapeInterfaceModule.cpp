#include <Rcpp.h>
using namespace Rcpp ;
#include "PersistenceLandscapeInterface.h"

RCPP_MODULE(class_PersistenceLandscape) {


    class_<PersistenceLandscapeInterface>("PersistenceLandscape")

    .constructor<NumericMatrix,bool,double,double,double,double>()


    .method("isExact", &PersistenceLandscapeInterface::isExact)
    .method("getMin", &PersistenceLandscapeInterface::getMin)
    .method("getMax", &PersistenceLandscapeInterface::getMax)
    .method("getdx", &PersistenceLandscapeInterface::getdx)
    .method("getExact", &PersistenceLandscapeInterface::getExact)
    .method("getDiscrete", &PersistenceLandscapeInterface::getDiscrete)
    .method("getInternal", &PersistenceLandscapeInterface::getInternal)
    .method("add", &PersistenceLandscapeInterface::add)
    .method("scale", &PersistenceLandscapeInterface::scale)
    .method("inner", &PersistenceLandscapeInterface::inner)
    ;
}
