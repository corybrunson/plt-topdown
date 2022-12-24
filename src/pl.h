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
  // TODO: Better defaults.
  // Natural defaults are inferred in R through `landscape()`.
  // Should C++ defaults be removed? -JCB
  PersistenceLandscapeInterface(
    NumericMatrix pd,
    bool exact = false,
    double min_pl = 0, double max_pl = 1,
    double dx = 0.01,
    double max_y = 1000)
    : exact(exact), min_pl(min_pl), max_pl(max_pl), dx(dx) {
    
    // Initalize a PersistenceLandscape object.
    std::vector<std::pair<double, double>> bars;
    
    for (int i = 0; i < pd.nrow(); i++)
      bars.push_back(std::make_pair(pd(i, 0), std::min(pd(i, 1), max_y)));
    // auto pb = PersistenceBarcodes(bars);
    
    PersistenceLandscapeInterface::pl_raw =
      PersistenceLandscape(bars, exact, min_pl, max_pl, 2 * dx);
  }
  
  PersistenceLandscapeInterface(
    PersistenceLandscape pl,
    bool exact,
    double min_pl, double max_pl,
    double dx)
    : pl_raw(pl), exact(exact), min_pl(min_pl), max_pl(max_pl), dx(dx) {}
  
  std::vector<NumericVector> getExact() {
    if (!exact) {
      stop("Error: Can not convert a discrete PL to an exact PL.");
    }
    else {
      return exactPersistenceLandscapeToR(pl_raw.land);
    }
  }
  
  NumericVector getDiscrete() {
    if (exact) {
      return discretePersistenceLandscapeToR(
        exactLandscapeToDiscrete(pl_raw.land, 0, max_pl, dx));
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
    return min_pl;
  }
  
  double getMax() const {
    return max_pl;
  }
  
  double getdx() const {
    return dx;
  }
  
  // Adds this to another PL
  PersistenceLandscapeInterface add(
      const PersistenceLandscapeInterface &other) {
    
    PersistenceLandscape pl_out;
    
    // both landscapes are exact
    if (PersistenceLandscapeInterface::exact && other.isExact())
      pl_out = pl_raw + other.pl_raw;
    
    // neither landscape is exact
    else if (!PersistenceLandscapeInterface::exact && !other.isExact()) {
      // FIXME: Drop these requirements once settings can be reconciled.
      if (!checkPairOfDiscreteLandscapes(*this, other)) {
        stop("Error: Persistence Landscape Properties Do Not Match.");
      }
      pl_out =
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
          other.min_pl,
          other.max_pl,
          other.dx);
        pl_out = PersistenceLandscape(addDiscreteLandscapes(
          conversion1,
          other.pl_raw));
      }
      
      // only second landscape is exact
      else if (conversions.second == true) {
        auto conversion2 = exactLandscapeToDiscrete(
          other.pl_raw,
          PersistenceLandscapeInterface::min_pl,
          PersistenceLandscapeInterface::max_pl,
          PersistenceLandscapeInterface::dx);
        pl_out = PersistenceLandscape(addDiscreteLandscapes(
          conversion2,
          PersistenceLandscapeInterface::pl_raw));
        // REVIEW: This is an alternative formulation that seems to still work
        // but not fix the bug.
        // pl_out = PersistenceLandscape(addDiscreteLandscapes(
        //   this->pl_raw,
        //   conversion2));
      }
    }
    
    // REVIEW: This has been edited from {tdatools} by JCB.
    // appropriate settings for output PL
    bool out_exact = exact & other.isExact();
    double out_min = min(min_pl, other.getMin());
    double out_max = max(max_pl, other.getMax());
    return PersistenceLandscapeInterface(pl_out,
                                         out_exact,
                                         out_min, out_max,
                                         dx);
  }
  
  double inner(
      PersistenceLandscapeInterface &other) {
    
    double scaler_out;
    
    if (exact)
      scaler_out = computeInnerProduct(pl_raw, other.pl_raw);
    else
      scaler_out = innerProductDiscreteLandscapes(pl_raw, other.pl_raw, dx);
    
    return scaler_out;
  }
  
  PersistenceLandscapeInterface scale(
      double scale) {
    
    PersistenceLandscape pl_out;
    
    if (exact)
      pl_out = scale * pl_raw;
    else
      pl_out = PersistenceLandscape(scaleDiscreteLandscapes(scale, pl_raw));
    
    return PersistenceLandscapeInterface(pl_out, exact, min_pl, max_pl, dx);
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
  double min_pl = 0;
  double max_pl = 1;
  double dx = 0.01;
  
};

PersistenceLandscapeInterface PLaverage(List p) {
  
  PersistenceLandscapeInterface out = as<PersistenceLandscapeInterface>(p[0]);
  
  for (int i = 1; i < p.size(); i++) {
    out = out.add(as<PersistenceLandscapeInterface>(p[i]));
  }
  
  return out.scale(1.0/p.size());
}

PersistenceLandscapeInterface PLsum(
    PersistenceLandscapeInterface p1,
    PersistenceLandscapeInterface p2) {
  return p1.add(p2);
}

PersistenceLandscapeInterface PLscale(
    double scale,
    PersistenceLandscapeInterface p) {
  return p.scale(scale);
}

double PLinner(
    PersistenceLandscapeInterface p1,
    PersistenceLandscapeInterface p2) {
  return p1.inner(p2);
}

bool checkPairOfDiscreteLandscapes(
    PersistenceLandscapeInterface &l1,
    const PersistenceLandscapeInterface &l2) {
  if (l1.min_pl != l2.min_pl || l1.max_pl != l2.max_pl || l1.dx != l2.dx)
    return false;
  return true;
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
