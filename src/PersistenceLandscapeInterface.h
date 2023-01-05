#include "PersistenceLandscape.h"
#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

int TDIndex(int X, int Y, int Z, int x, int y, int z) {
  return x + X * (y + Y * z);
}

int TDIndex2(int X, int Y, int x, int y) {
  return x + X * y;
}

NumericVector discretePersistenceLandscapeToR(
    std::vector<std::vector<std::pair<double, double>>> input) {
  
  Dimension d(input.size(), input[0].size(), 2);
  NumericVector out(input.size() * input[0].size() * 2);
  NumericVector out_d(d);
  
  for (int j = 0; j < input.size(); j++) {
    for (int i = 0; i < input[0].size(); i++) {
      out[TDIndex(input.size(), input[0].size(), 2, j, i, 0)] =
        input[j][i].first;
      out[TDIndex(input.size(), input[0].size(), 2, j, i, 1)] =
        input[j][i].second;
      
      // REVIEW: Why is this needed for discrete landscapes? -JCB
      if (input[j][i].first == INT_MAX)
        out[TDIndex(input.size(), input[0].size(), 2, j, i, 0)] = R_PosInf;
      if (input[j][i].first == INT_MIN)
        out[TDIndex(input.size(), input[0].size(), 2, j, i, 0)] = R_NegInf;
    }
  }
  
  std::copy(out.begin(), out.end(), out_d.begin());
  return out_d;
}

std::vector<NumericVector> exactPersistenceLandscapeToR(
    std::vector<std::vector<std::pair<double, double>>> input) {
  
  std::vector<NumericVector> out_d;
  
  for (int j = 0; j < input.size(); j++) {
    Dimension d(input[j].size(), 2);
    NumericVector out(input[j].size() * 2);
    
    for (int i = 0; i < input[j].size(); i++) {
      out[TDIndex2(input[j].size(), 2, i, 0)] =
        input[j][i].first;
      out[TDIndex2(input[j].size(), 2, i, 1)] =
        input[j][i].second;
      
      if (input[j][i].first == INT_MAX)
        out[TDIndex2(input[j].size(), 2, i, 0)] = R_PosInf;
      
      if (input[j][i].first == INT_MIN)
        out[TDIndex2(input[j].size(), 2, i, 0)] = R_NegInf;
    }
    
    out.attr("dim") = d;
    out_d.push_back(out);
  }
  
  return out_d;
}

std::vector<std::vector<std::pair<double, double>>> addDiscreteLandscapes(
    const PersistenceLandscape &l1,
    const PersistenceLandscape &l2) {
  
  int min_level = std::min(l1.land.size(), l2.land.size());
  std::vector<std::vector<std::pair<double, double>>> out;
  for (int i = 0; i < min_level; i++) {
    int min_index = std::min(l1.land[i].size(), l2.land[i].size());
    std::vector<std::pair<double, double>> level_out;
    for (int j = 0; j < min_index; j++)
      level_out.push_back(std::make_pair(
          l1.land[i][j].first,
          l1.land[i][j].second + l2.land[i][j].second));
    
    for (; min_index < l1.land[i].size(); min_index++)
      level_out.push_back(l1.land[i][min_index]);
    
    for (; min_index < l2.land[i].size(); min_index++)
      level_out.push_back(l2.land[i][min_index]);
    
    out.push_back(level_out);
  }
  
  for (; min_level < l1.land.size(); min_level++)
    out.push_back(l1.land[min_level]);
  
  for (; min_level < l2.land.size(); min_level++)
    out.push_back(l2.land[min_level]);
  
  return out;
}

std::vector<std::vector<std::pair<double, double>>> scaleDiscreteLandscapes(
    double scale,
    PersistenceLandscape l) {
  
  std::vector<std::vector<std::pair<double, double>>> out;
  
  for (std::vector<std::pair<double, double>> level : l.land) {
    std::vector<std::pair<double, double>> level_out;
    for (std::pair<double, double> pair : level)
      level_out.push_back(std::make_pair(pair.first, scale * pair.second));
    // for (auto i : level_out) {
    //   // TODO: Delete this loop or figure out what belongs in it! -JCB
    // }
    out.push_back(level_out);
  }
  
  return out;
}

std::vector<std::vector<std::pair<double, double>>> exactLandscapeToDiscrete(
    PersistenceLandscape l,
    double min_x,
    double max_x,
    double dx) {
  
  std::vector<std::vector<std::pair<double, double>>> out;
  
  for (unsigned i = 0; i < l.land.size(); i++) {
    
    auto level = l.land[i];
    std::vector<std::pair<double, double>> level_out;
    
    double starting_y = l.computeValueAtAGivenPoint(i, min_x);
    std::pair<double, double> startingPoint = std::make_pair(min_x, starting_y);
    double x_buffer = startingPoint.first;
    double y_buffer = startingPoint.second;
    
    for (int i = 1; i < level.size(); i++) {
      std::pair<double, double> point = level[i];
      // Make sure slope is well defined:
      if (point.first != startingPoint.first) {
        double delta_x = level[i].first - startingPoint.first;
        double delta_y = level[i].second - startingPoint.second;
        double slope = delta_y / delta_x;
        
        while (x_buffer < point.first && x_buffer < max_x) {
          y_buffer += dx * slope;
          x_buffer += dx;
          level_out.push_back(std::make_pair(x_buffer, y_buffer));
        }
      }
      
      startingPoint = point;
    }
    
    out.push_back(level_out);
  }
  
  return out;
}

double innerProductDiscreteLandscapes(
    PersistenceLandscape l1,
    PersistenceLandscape l2,
    double dx) {
  
  int min_level = std::min(l1.land.size(), l2.land.size());
  double integral_buffer = 0;
  
  for (int i = 0; i < min_level; i++) {
    int min_index = std::min(l1.land[i].size(), l2.land[i].size());
    for (int j = 0; j < min_index; j++)
      integral_buffer += l1.land[i][j].second * l2.land[i][j].second;
  }
  
  return integral_buffer * dx;
}

std::vector<std::pair<double, double>> generateGrid(
    double start,
    double end,
    double dx) {
  
  std::vector<std::pair<double, double>> grid;
  
  for (double current = start; current < end; current += dx)
    grid.push_back(std::make_pair(current, 0));
  
  return grid;
}

class PersistenceLandscapeInterface {
public:
  // Creates PL from PD (in the form of a 2-column numeric matrix)
  // Natural defaults are inferred in R through `landscape()`.
  PersistenceLandscapeInterface(
    NumericMatrix pd,
    bool exact = false,
    double min_x = 0, double max_x = 1,
    double dx = 0.01,
    double max_y = 1000)
    : exact(exact), min_x(min_x), max_x(max_x), dx(dx) {
    
    // Initialize a PersistenceLandscape object.
    std::vector<std::pair<double, double>> bars;
    
    for (int i = 0; i < pd.nrow(); i++)
      bars.push_back(std::make_pair(pd(i, 0), std::min(pd(i, 1), max_y)));
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
  
  std::vector<NumericVector> toExact() {
    if (!exact) {
      stop("Error: Can not convert a discrete PL to an exact PL.");
    }
    else {
      return exactPersistenceLandscapeToR(pl_raw.land);
    }
  }
  
  NumericVector toDiscrete() {
    if (exact) {
      return discretePersistenceLandscapeToR(
        exactLandscapeToDiscrete(pl_raw.land, 0, max_x, dx));
    } else {
      // TODO: Allow this function to change the resolution.
      return discretePersistenceLandscapeToR(pl_raw.land);
    }
  }
  
  SEXP getInternal() {
    if (PersistenceLandscapeInterface::exact == false)
      return wrap(discretePersistenceLandscapeToR(pl_raw.land));
    else
      return wrap(exactPersistenceLandscapeToR(pl_raw.land));
  }
  
  bool isExact() const {
    return exact;
  }
  
  double getMin() const {
    return min_x;
  }
  
  double getMax() const {
    return max_x;
  }
  
  double getdx() const {
    return dx;
  }
  
  // Adds this to another PL
  PersistenceLandscapeInterface add(
      const PersistenceLandscapeInterface &other) {
    
    PersistenceLandscape add_out;
    
    // both landscapes are exact
    if (PersistenceLandscapeInterface::exact && other.isExact())
      add_out = pl_raw + other.pl_raw;
    
    // neither landscape is exact
    else if (!PersistenceLandscapeInterface::exact && !other.isExact()) {
      // FIXME: Drop these requirements once settings can be reconciled.
      if (!checkPairOfDiscreteLandscapes(*this, other)) {
        stop("Error: Persistence Landscape Properties Do Not Match.");
      }
      add_out =
        PersistenceLandscape(addDiscreteLandscapes(pl_raw, other.pl_raw));
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
    // appropriate settings for output PL
    bool out_exact = exact & other.isExact();
    double out_min = min(min_x, other.getMin());
    double out_max = max(max_x, other.getMax());
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
      scale_out = PersistenceLandscape(scaleDiscreteLandscapes(scale, pl_raw));
    
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
  
  double supremum(unsigned level) {
    return pl_raw.findMax(level);
  }
  
  double infimum(unsigned level) {
    return pl_raw.findMin(level);
  }
  
  double moment(
      unsigned n,
      double center,
      unsigned level) {
    
    double moment_out = pl_raw.computeNthMoment(n, center, level);
    
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

bool checkPairOfDiscreteLandscapes(
    PersistenceLandscapeInterface &l1,
    const PersistenceLandscapeInterface &l2) {
  if (l1.min_x != l2.min_x || l1.max_x != l2.max_x || l1.dx != l2.dx)
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
