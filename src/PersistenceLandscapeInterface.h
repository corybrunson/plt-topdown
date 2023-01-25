#include "PersistenceLandscape.h"
#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

// NOTE: Indexing of (level/envelope) arguments must at some point be shifted
// from starting at 1 (in R) to starting at 0 (in C++). Since the
// 'PersistenceLandscapeInterface' methods are exposed to R, this shift will be
// done in these methods, before they are passed to 'PersistenceLandscape'
// methods. -JCB

class PersistenceLandscapeInterface {
  
public:
  
  // Creates a PL from a PD, in the form of a 2-column numeric matrix.
  // Natural defaults are inferred in R through `landscape()`.
  PersistenceLandscapeInterface(
    NumericMatrix pd,
    bool exact = false,
    double min_x = 0, double max_x = 1,
    double dx = 0.01,
    double min_y = R_NegInf, double max_y = R_PosInf)
    : exact(exact), min_x(min_x), max_x(max_x), dx(dx) {
    
    // Initialize a PersistenceLandscape object.
    std::vector<std::pair<double, double>> bars;
    
    for (int i = 0; i < pd.nrow(); i++)
      bars.push_back(std::make_pair(pd(i,0),
                                    std::max(std::min(pd(i,1), max_y), min_y)));
    // auto pb = PersistenceBarcodes(bars);
    
    PersistenceLandscapeInterface::pl_raw =
      PersistenceLandscape(bars, exact, min_x, max_x, 2 * dx);
  }
  
  PersistenceLandscapeInterface(
    PersistenceLandscape pl,
    bool exact,
    double min_x, double max_x,
    double dx)
    : pl_raw(pl), exact(exact), min_x(min_x), max_x(max_x), dx(dx) {}
  
  bool isExact() const {
    return exact;
  }
  
  double xMin() const {
    return min_x;
  }
  
  double xMax() const {
    return max_x;
  }
  
  double xBy() const {
    return dx;
  }
  
  SEXP getInternal() {
    if (PersistenceLandscapeInterface::exact == false)
      return wrap(discretePersistenceLandscapeToR(pl_raw.land));
    else
      return wrap(exactPersistenceLandscapeToR(pl_raw.land));
  }
  
  NumericVector toDiscrete() {
    if (exact) {
      return discretePersistenceLandscapeToR(
        exactLandscapeToDiscrete(pl_raw.land, min_x, max_x, dx));
    } else {
      return discretePersistenceLandscapeToR(pl_raw.land);
    }
  }
  
  std::vector<NumericVector> toExact() {
    if (! exact) {
      stop("Can not convert a discrete PL to an exact PL.");
    } else {
      return exactPersistenceLandscapeToR(pl_raw.land);
    }
  }
  
  PersistenceLandscapeInterface discretize(
      double min_x, double max_x,
      double dx) {
    
    PersistenceLandscape out;
    
    if (exact) {
      out = exactLandscapeToDiscrete(pl_raw.land, min_x, max_x, dx);
    } else {
      warning("Can not yet re-discretize a discrete PL.");
      out = pl_raw.land;
    }
    
    return PersistenceLandscapeInterface(out, false, min_x, max_x, dx);
  }
  
  // REVIEW: Expands the limits of a PL -JCB
  PersistenceLandscapeInterface expand(
      double min_x, double max_x) {
    
    PersistenceLandscape out;
    
    if (exact)
      out = pl_raw;
    else
      out = PersistenceLandscape(expandDiscreteLandscape(
        pl_raw,
        min_x, max_x,
        dx));
    
    return PersistenceLandscapeInterface(out, exact, min_x, max_x, dx);
  }
  
  // REVIEW: Contracts the limits of a PL -JCB
  PersistenceLandscapeInterface contract(
      double min_x, double max_x) {
    
    PersistenceLandscape out;
    
    if (exact)
      out = pl_raw;
    else
      out = PersistenceLandscape(contractDiscreteLandscape(
        pl_raw,
        min_x, max_x,
        dx));
    
    return PersistenceLandscapeInterface(out, exact, min_x, max_x, dx);
  }
  
  // Adds this to another PL
  PersistenceLandscapeInterface add(
      const PersistenceLandscapeInterface &other) {
    
    PersistenceLandscape add_out;
    // appropriate settings for output PL
    bool out_exact = exact & other.isExact();
    double out_min = min(min_x, other.min_x);
    double out_max = max(max_x, other.max_x);
    
    // both landscapes are exact
    if (PersistenceLandscapeInterface::exact && other.isExact())
      add_out = pl_raw + other.pl_raw;
    
    // neither landscape is exact
    else if (! PersistenceLandscapeInterface::exact && ! other.isExact()) {
      if (! checkPairOfDiscreteLandscapes(*this, other)) {
        stop("Resolutions and limits are incompatible.");
      }
      // REVIEW: Allow more landscapes to be added. Take output limits to be
      // union of input limits. -JCB
      if (alignPairOfDiscreteLandscapes(*this, other)) {
        add_out =
          PersistenceLandscape(addDiscreteLandscapes(pl_raw, other.pl_raw));
      } else {
        add_out =
          PersistenceLandscape(addDiscreteLandscapes(
              expandDiscreteLandscape(pl_raw, out_min, out_max, dx),
              expandDiscreteLandscape(other.pl_raw, out_min, out_max, dx)));
      }
    }
    
    else{
      // Conversions:
      std::pair<bool, bool> conversions =
        operationOnPairOfLanscapesConversion(*this, other);
      
      // only first landscape is exact
      if (conversions.first == true) {
        auto conversion1 = exactLandscapeToDiscrete(
          this->pl_raw,
          other.min_x,
          other.max_x,
          other.dx);
        add_out = PersistenceLandscape(addDiscreteLandscapes(
          conversion1,
          other.pl_raw));
      }
      
      // only second landscape is exact
      else if (conversions.second == true) {
        auto conversion2 = exactLandscapeToDiscrete(
          other.pl_raw,
          PersistenceLandscapeInterface::min_x,
          PersistenceLandscapeInterface::max_x,
          PersistenceLandscapeInterface::dx);
        add_out = PersistenceLandscape(addDiscreteLandscapes(
          conversion2,
          PersistenceLandscapeInterface::pl_raw));
        // REVIEW: This is an alternative formulation that seems to still work
        // but not fix the bug.
        // add_out = PersistenceLandscape(addDiscreteLandscapes(
        //   this->pl_raw,
        //   conversion2));
      }
    }
    
    // REVIEW: This has been edited from {tdatools} by JCB.
    return PersistenceLandscapeInterface(add_out,
                                         out_exact,
                                         out_min, out_max,
                                         dx);
  }
  
  PersistenceLandscapeInterface scale(
      double scale) {
    
    PersistenceLandscape scale_out;
    
    if (exact)
      scale_out = scale * pl_raw;
    else
      scale_out = PersistenceLandscape(scaleDiscreteLandscape(pl_raw, scale));
    
    return PersistenceLandscapeInterface(scale_out, exact, min_x, max_x, dx);
  }
  
  PersistenceLandscapeInterface abs() {
    
    PersistenceLandscape abs_out = pl_raw.abs();
    
    return PersistenceLandscapeInterface(abs_out, exact, min_x, max_x, dx);
  }
  
  double inner(
      PersistenceLandscapeInterface &other) {
    
    double scalar_out;
    
    if (exact)
      scalar_out = computeInnerProduct(pl_raw, other.pl_raw);
    else
      scalar_out = innerProductDiscreteLandscapes(pl_raw, other.pl_raw, dx);
    
    return scalar_out;
  }
  
  double minimum(unsigned level) {
    if (level <= 0 || level > pl_raw.land.size())
      return NA_REAL;
    return pl_raw.findMin(level - 1);
  }
  
  double maximum(unsigned level) {
    if (level <= 0 || level > pl_raw.land.size())
      return NA_REAL;
    return pl_raw.findMax(level - 1);
  }
  
  double moment(
      unsigned n,
      double center,
      unsigned level) {
    
    if (level <= 0 || level > pl_raw.land.size())
      return NA_REAL;
    
    double moment_out = pl_raw.computeNthMoment(n, center, level - 1);
    
    return moment_out;
  }
  
  double integral(
      unsigned p) {
    
    double int_out;
    
    if (p == 1)
      int_out = pl_raw.computeIntegralOfLandscape();
    else
      int_out = pl_raw.computeIntegralOfLandscape(p);
    
    return int_out;
  }
  
  double distance(
      PersistenceLandscapeInterface &other,
      unsigned p) {
    
    double dist_out;
    
    if (p == 0)
      // `p = 0` encodes `p = Inf`
      dist_out = computeMaxNormDistanceBetweenLandscapes(pl_raw, other.pl_raw);
    else
      dist_out = computeDistanceBetweenLandscapes(pl_raw, other.pl_raw, p);
    
    return dist_out;
  }
  
  double norm(
      unsigned p) {
    
    double norm_out = pl_raw.computeNormOfLandscape(p);
    
    return norm_out;
  }
  
  PersistenceLandscapeInterface indicator(
      List indicator,
      unsigned r) {
    
    // Encode the list of vectors as a vector of pairs.
    std::vector<std::pair<double, double>> ind;
    for (size_t i = 0; i != indicator.length(); ++i) {
      std::vector<double> supp = indicator[i];
      ind.push_back(std::make_pair(supp[0], supp[1]));
    }
    
    PersistenceLandscape pl_out = pl_raw.multiplyByIndicatorFunction(ind, r);
    
    return PersistenceLandscapeInterface(pl_out, exact, min_x, max_x, dx);
  }
  
  double indicator_form(
      List indicator,
      unsigned r,
      unsigned p) {
    
    // Encode the list of vectors as a vector of pairs.
    std::vector<std::pair<double, double>> ind;
    for (size_t i = 0; i != indicator.length(); ++i) {
      std::vector<double> supp = indicator[i];
      ind.push_back(std::make_pair(supp[0], supp[1]));
    }
    
    double form_out;
    
    if (p == 1)
      form_out = pl_raw.computeIntegralOfLandscapeMultipliedByIndicatorFunction(
        ind, r);
    else
      form_out = pl_raw.computeIntegralOfLandscapeMultipliedByIndicatorFunction(
        ind, r, p);
    
    return form_out;
  }
  
  friend bool checkPairOfDiscreteLandscapes(
      PersistenceLandscapeInterface &l1,
      const PersistenceLandscapeInterface &l2);
  
  friend bool alignPairOfDiscreteLandscapes(
      PersistenceLandscapeInterface &l1,
      const PersistenceLandscapeInterface &l2);
  
  friend std::pair<bool, bool> operationOnPairOfLanscapesConversion(
      PersistenceLandscapeInterface &l1,
      const PersistenceLandscapeInterface &l2);
  
private:
  
  PersistenceLandscape pl_raw;
  bool exact;
  double min_x = 0;
  double max_x = 1;
  double dx = 0.01;
  
};

// REVIEW: Check only that resolutions are equal and that endpoints lie on the
// same grid of the shared resolution. -JCB
bool checkPairOfDiscreteLandscapes(
    PersistenceLandscapeInterface &l1,
    const PersistenceLandscapeInterface &l2) {
  // if (l1.min_x != l2.min_x || l1.max_x != l2.max_x || l1.dx != l2.dx)
  //   return false;
  if (l1.dx != l2.dx) return false;
  
  double mod_dx = fmod(l1.min_x - l2.min_x, l1.dx);
  // WARNING: Allow flexibility here since `addDiscreteLandscapes()` simply
  // takes the first input's abscissa. -JCB
  if (! almostEqual(mod_dx, 0.) && ! almostEqual(mod_dx, l1.dx)) return false;
  
  return true;
}

// REVIEW: Provided parameters are compatible, check whether bounds need to be
// aligned. -JCB
bool alignPairOfDiscreteLandscapes(
    PersistenceLandscapeInterface &l1,
    const PersistenceLandscapeInterface &l2) {
  if (l1.min_x != l2.min_x || l1.max_x != l2.max_x)
    return false;
  return true;
}

PersistenceLandscapeInterface PLsum(List pl_list) {
  
  PersistenceLandscapeInterface
  sum_out = as<PersistenceLandscapeInterface>(pl_list[0]);
  
  for (int i = 1; i < pl_list.size(); i++) {
    sum_out = sum_out.add(as<PersistenceLandscapeInterface>(pl_list[i]));
  }
  
  return sum_out;
}

List PLdiff(List pl_list) {
  
  List diff_out;
  
  for (int i = 1; i < pl_list.size(); i++) {
    PersistenceLandscapeInterface
    diff_i = as<PersistenceLandscapeInterface>(pl_list[i]);
    diff_i = diff_i.add(as<PersistenceLandscapeInterface>(
      pl_list[i - 1]).scale(-1));
    diff_out.push_back(diff_i);
  }
  
  return diff_out;
}

PersistenceLandscapeInterface PLmean(List pl_list) {
  
  PersistenceLandscapeInterface avg_out = PLsum(pl_list);
  
  return avg_out.scale(1.0 / pl_list.size());
}

NumericMatrix PLdist(List pl_list, unsigned p) {
  
  // empty matrix
  unsigned n = pl_list.size();
  NumericMatrix out(n, n);
  
  // not assuming symmetric distance calculation
  for (int i = 0; i != n; i++) {
    PersistenceLandscapeInterface
    pl_i = as<PersistenceLandscapeInterface>(pl_list[i]);
    for (int j = 0; j != n; j++) {
      if (j == i) {
        // same landscape
        out(i, j) = 0.;
      } else {
        // different landscape
        PersistenceLandscapeInterface
        pl_j = as<PersistenceLandscapeInterface>(pl_list[j]);
        out(i, j) = pl_i.distance(pl_j, p);
      }
    }
  }
  
  return out;
}

double PLvar(List pl_list, unsigned p) {
  
  // average landscape
  PersistenceLandscapeInterface avg = PLmean(pl_list);
  
  // sum-squared distance
  double ssd = 0;
  
  for (size_t i = 0; i != pl_list.size(); ++i) {
    
    PersistenceLandscapeInterface
    pl_i = as<PersistenceLandscapeInterface>(pl_list[i]);
    
    double d = avg.distance(pl_i, p);
    
    ssd += d * d;
  }
  
  // sample standard deviation
  double var_out = ssd / pl_list.size();
  // double var_out = ssd / (pl_list.size() - 1.0);
  return var_out;
}

double PLsd(List pl_list, unsigned p) {
  
  double sd_out = PLvar(pl_list, p);
  
  sd_out = sqrt(sd_out);
  return sd_out;
}

// For operations on two landscapes we need to know if the output will be
// discrete or exact. The rules are:
// Exact + Exact = Exact
// Discrete + Discrete = Discrete
// Exact + Discrete = Discrete
// This function is for the third case, where we will need to convert one of the
// exact PLs to a discrete PL. The two bools correspond to the two operands, and
// are 1 (true) if we need to convert the operand to discrete. In all other
// cases it is 0 (false).
std::pair<bool,bool> operationOnPairOfLanscapesConversion(
    PersistenceLandscapeInterface &pl1,
    const PersistenceLandscapeInterface &pl2) {
  
  bool p1 = 0;
  bool p2 = 0;
  
  if (pl1.exact != pl2.exact) {
    if (pl1.exact == true)
      p1 = 1;
    else
      // REVIEW: This has been edited from {tdatools} by JCB.
      p2 = 1;
  }
  
  return std::make_pair(p1,p2);
}
