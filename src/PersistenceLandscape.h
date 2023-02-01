//    AOAOAOA
//    AOA
//    AOA
//    AOA
//    AOA
//    AOA
//    A
//    Copyright 2013-2014 University of Pennsylvania
//    Created by Pawel Dlotko
//
//    This file is part of Persistence Landscape Toolbox (PLT).
//
//    PLT is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published
//    by the Free Software Foundation, either version 3 of the License, or (at
//    your option) any later version.
//
//    PLT is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with PLT.  If not, see <http://www.gnu.org/licenses/>.

//    This file has been substantially modified by Jason Cory Brunson in
//    consultation with Dlotko.

#ifndef PERISTENCELANDSCAPE_H
#define PERISTENCELANDSCAPE_H

#include <Rcpp.h>
#include <limits>

using namespace std;
using namespace Rcpp;

// TODO: Make this a user-settable option through {Rcpp}.
double epsi = 0.000005;

inline bool almostEqual(double a, double b) {
  if (fabs(a - b) < epsi)
    return true;
  return false;
}

double birth(std::pair<double, double> a) { return a.first - a.second; }

double death(std::pair<double, double> a) { return a.first + a.second; }

// this function assumes birth-death coordinates
bool comparePoints2(
    std::pair<double, double> f,
    std::pair<double, double> s) {
  if (f.first < s.first) {
    return true;
  } else {
    // f.first >= s.first
    if (f.first > s.first) {
      return false;
    } else {
      // f.first == s.first
      if (f.second > s.second) {
        return true;
      } else {
        return false;
      }
    }
  }
}

// named functions for infix operators
inline double add(double x, double y) { return x + y; }
inline double sub(double x, double y) { return x - y; }

// get linear function value at a point
double functionValue(
    std::pair<double, double> p1,
    std::pair<double, double> p2,
    double x) {
  // we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points
  // between which we will put the line segment
  double a = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - a * p1.first;
  
  return (a * x + b);
}

std::vector<std::pair<double, double>> generateGrid(
    double start,
    double end,
    double dx) {
  
  std::vector<std::pair<double, double>> grid;
  
  for (double current = start; current < end; current += dx)
    grid.push_back(std::make_pair(current, 0.));
  
  return grid;
}

// Persistence landscape data transformations

// The following functions transform the data in which persistence landscapes
// are encoded, for export to R.

int TDIndex2(int X, int Y, int x, int y) {
  return x + X * y;
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

int TDIndex(int X, int Y, int Z, int x, int y, int z) {
  return x + X * (y + Y * z);
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

// Modified by Jose Bouza to accept exact/discrete constructor param.
class PersistenceLandscape {
  
public:
  
  // This is the public data member of the 'PersistenceLandscape' class. It is a
  // vector of levels, each of which is a vector of x,y-pairs, each of which is
  // a double. If the landscape is exact, then the levels may have different
  // sizes; if the landscape is discrete, then the levels have the same size and
  // the same x-values.
  std::vector<std::vector<std::pair<double, double>>> land;
  
  // This member function provides a shortcut for the size of a PL, interpreted
  // as the number of its levels.
  size_t size() const { return this->land.size(); }
  
  PersistenceLandscape(){}
  
  // Constructors (member functions that instantiate new class objects):
  
  // This constructor is specified here in the class and defined below. It takes
  // as its primary input a persistence diagram (data set) comprising
  // birth-death pairs.
  PersistenceLandscape(
    const std::vector<std::pair<double, double>> &diagram,
    bool exact = true,
    double min_x = 0, double max_x = 1,
    double dx = 0.01);
  
  // This constructor is both specified and defined here in the class.
  PersistenceLandscape(
    const NumericMatrix pd,
    bool exact = true,
    double min_x = 0, double max_x = 1,
    double dx = 0.01,
    double min_y = R_NegInf, double max_y = R_PosInf)
    : exact(exact), min_x(min_x), max_x(max_x), dx(dx) {
    
    // Convert numeric matrix in R to PD structure in C++.
    std::vector<std::pair<double, double>> diag;
    for (int i = 0; i < pd.nrow(); i++)
      diag.push_back(std::make_pair(pd(i,0),
                                    std::max(std::min(pd(i,1), max_y), min_y)));
    
    // PersistenceLandscape::diagram = diag;
    PersistenceLandscape(diag, exact, min_x, max_x, dx);
  };
  
  PersistenceLandscape(const PersistenceLandscape &original);
  
  // assignment operator overload
  PersistenceLandscape operator=(const PersistenceLandscape &original);
  
  PersistenceLandscape(
    std::vector<std::vector<std::pair<double, double>>>
    landscapePointsWithoutInfinities);
  
  double computeValueAtAGivenPoint(
      unsigned level,
      double x) const;
  
  PersistenceLandscape scalePersistenceLandscape(double x) const;
  
  unsigned removePairsOfLocalMaximumMinimumOfEpsPersistence(
      double errorTolerance);
  void reduceAllPairsOfLowPersistenceMaximaMinima(
      double epsilon);
  void reduceAlignedPoints(
      double tol = 0.000001);
  unsigned reducePoints(
      double tol,
      double (*penalty)(
          std::pair<double, double>,
          std::pair<double, double>,
          std::pair<double, double>));
  
  // Friendzone:
  
  // This is a general algorithm to perform linear operations on persistence
  // landscapes. It perform it by doing operations on landscape points.
  friend PersistenceLandscape operationOnTwoExactLandscapes(
      const PersistenceLandscape &land1,
      const PersistenceLandscape &land2,
      double (*oper)(double, double));
  
  friend PersistenceLandscape operationOnTwoLandscapes(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2,
      double (*oper)(double, double));
  
  friend PersistenceLandscape addTwoLandscapes(
      const PersistenceLandscape &land1,
      const PersistenceLandscape &land2) {
    return operationOnTwoLandscapes(land1, land2, add);
  }
  
  friend PersistenceLandscape subtractTwoLandscapes(
      const PersistenceLandscape &land1,
      const PersistenceLandscape &land2) {
    return operationOnTwoLandscapes(land1, land2, sub);
  }
  
  friend PersistenceLandscape operator+(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second) {
    return addTwoLandscapes(first, second);
  }
  
  friend PersistenceLandscape operator-(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second) {
    return subtractTwoLandscapes(first, second);
  }
  
  friend PersistenceLandscape operator*(
      const PersistenceLandscape &first,
      double con) {
    return first.scalePersistenceLandscape(con);
  }
  
  friend PersistenceLandscape operator*(
      double con,
      const PersistenceLandscape &first) {
    return first.scalePersistenceLandscape(con);
  }
  
  friend PersistenceLandscape operator/(
      const PersistenceLandscape &first,
      double con) {
    return first.scalePersistenceLandscape(1 / con);
  }
  
  PersistenceLandscape operator+=(const PersistenceLandscape &rhs) {
    *this = *this + rhs;
    return *this;
  }
  
  PersistenceLandscape operator-=(const PersistenceLandscape &rhs) {
    *this = *this - rhs;
    return *this;
  }
  
  PersistenceLandscape operator*=(double x) {
    *this = *this * x;
    return *this;
  }
  
  PersistenceLandscape operator/=(double x) {
    if (x == 0)
      throw("In operator /=, division by 0. Program terminated.");
    *this = *this * (1 / x);
    return *this;
  }
  
  // This operator is defined outside this declaration:
  // `PersistenceLandscape::operator==`
  bool operator==(const PersistenceLandscape &rhs) const;
  
  // integral of (the p^th power of) a landscape
  double computeIntegralOfLandscape() const;
  double computeIntegralOfLandscape(double p) const;
  
  PersistenceLandscape multiplyByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator,
      unsigned r) const;
  PersistenceLandscape multiplyByIndicatorFunction(
      List indicator,
      unsigned r) const;
  
  // integral of (the p^th power of) the product of a landscape with an
  // indicator function
  double computeIntegralOfLandscapeMultipliedByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator,
      unsigned r,
      double p) const;
  double computeIntegralOfLandscapeMultipliedByIndicatorFunction(
      List indicator,
      unsigned r,
      double p) const;
  
  double computeNormOfLandscape(int i) {
    PersistenceLandscape l;
    if (i != -1) {
      return computeFinNormDistanceBetweenLandscapes(*this, l, i);
    } else {
      return computeInfNormDistanceBetweenLandscapes(*this, l);
    }
  }
  
  double operator()(unsigned level, double x) const {
    return this->computeValueAtAGivenPoint(level, x);
  }
  
  friend bool checkPairOfDiscreteLandscapes(
      const PersistenceLandscape &l1,
      const PersistenceLandscape &l2);
  
  friend bool alignPairOfDiscreteLandscapes(
      const PersistenceLandscape &l1,
      const PersistenceLandscape &l2);
  
  friend std::pair<bool, bool> operationOnPairOfLanscapesConversion(
      const PersistenceLandscape &l1,
      const PersistenceLandscape &l2);
  
  friend double computeFinNormDistanceBetweenLandscapes(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second,
      unsigned p);
  
  friend double computeMaximalDistanceNonSymmetric(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2);
  
  friend double computeMaximalDistanceNonSymmetric(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2,
      unsigned &nrOfLand,
      double &x,
      double &y1, double &y2);
  // This function additionally returns integer n and double x, y1, y2 such that
  // the maximal distance is obtained between lambda_n's on a coordinate x such
  // that the value of the first landscape is y1, and the vale of the second
  // landscape is y2.
  
  friend double computeInfNormDistanceBetweenLandscapes(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second);
  
  friend double computeInfNormDistanceBetweenLandscapes(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second,
      unsigned &nrOfLand,
      double &x,
      double &y1, double &y2);
  
  friend double computeDistanceBetweenLandscapes(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second,
      unsigned p);
  
  PersistenceLandscape abs();
  
  double minimum(unsigned lambda) const;
  
  double maximum(unsigned lambda) const;
  
  friend double computeInnerProduct(
      const PersistenceLandscape &l1,
      const PersistenceLandscape &l2);
  
  double moment(
      unsigned p,
      double center,
      unsigned level) const;
  
  // Copied or adapted from class 'PersistenceLandscapeInterface':
  
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
    if (exact)
      return wrap(exactPersistenceLandscapeToR(this->land));
    else
      return wrap(discretePersistenceLandscapeToR(this->land));
  }
  
  NumericVector toDiscrete();
  
  std::vector<NumericVector> toExact();
  
  PersistenceLandscape discretize(
      double min_x, double max_x,
      double dx);
  
  // Expands the limits of a PL -JCB
  PersistenceLandscape expand(
      double min_x, double max_x);
  
  // Contracts the limits of a PL -JCB
  PersistenceLandscape contract(
      double min_x, double max_x);
  
private:
  
  // REVIEW: The PD should not be preserved in the PL object. -JCB
  std::vector<std::pair<double, double>> diagram;
  bool exact;
  double min_x;
  double max_x;
  double dx;
  
};

PersistenceLandscape::PersistenceLandscape(
  const std::vector<std::pair<double, double>> &diagram,
  bool exact,
  double min_x, double max_x,
  double grid_diameter) {
  
  // Adapted from constructor:
  // `PersistenceBarcodes(std::vector<std::pair<double, double>> bars)`
  std::vector<std::pair<double, double>> pd = diagram;
  unsigned nb = 0;
  for (size_t i = 0; i != pd.size(); ++i) {
    if (pd[i].second != R_PosInf & pd[i].second != R_NegInf) {
      ++nb;
    }
    // FIXME: Harmonize this step with extended persistence data. -JCB
    if (pd[i].second < pd[i].first) {
      double sec = pd[i].second;
      pd[i].second = pd[i].first;
      pd[i].first = sec;
    }
  }
  unsigned nr = 0;
  for (size_t i = 0; i != pd.size(); ++i) {
    if (pd[i].second != R_PosInf & pd[i].second != R_NegInf) {
      // This is a finite interval.
      pd[nr] = std::make_pair(pd[i].first, pd[i].second);
      ++nr;
    }
  }
  
  if (exact) {
    // This is a general algorithm to construct persistence landscapes.
    
    std::vector<std::pair<double, double>> pds;
    pds.insert(pds.begin(), pd.begin(), pd.end());
    std::sort(pds.begin(), pds.end(), comparePoints2);
    
    std::vector<std::pair<double, double>> characteristicPoints(pds.size());
    
    for (size_t i = 0; i != pds.size(); ++i) {
      characteristicPoints[i] =
        std::make_pair((pds[i].first + pds[i].second) / 2.0,
                       (pds[i].second - pds[i].first) / 2.0);
    }
    
    std::vector<std::vector<std::pair<double, double>>> persistenceLandscape;
    while (!characteristicPoints.empty()) {
      std::vector<std::pair<double, double>> lambda_n;
      lambda_n.push_back(std::make_pair(INT_MIN, 0.));
      lambda_n.push_back(std::make_pair(birth(characteristicPoints[0]), 0.));
      lambda_n.push_back(characteristicPoints[0]);
      
      int i = 1;
      std::vector<std::pair<double, double>> newCharacteristicPoints;
      while (i < characteristicPoints.size()) {
        size_t j = 1;
        if ((birth(characteristicPoints[i]) >=
            birth(lambda_n[lambda_n.size() - 1])) &&
            (death(characteristicPoints[i]) >
            death(lambda_n[lambda_n.size() - 1]))) {
          if (birth(characteristicPoints[i]) <
            death(lambda_n[lambda_n.size() - 1])) {
            std::pair<double, double> point =
              std::make_pair((birth(characteristicPoints[i]) +
              death(lambda_n[lambda_n.size() - 1])) /
                2.,
                (death(lambda_n[lambda_n.size() - 1]) -
                  birth(characteristicPoints[i])) /
                    2.);
            lambda_n.push_back(point);
            
            while ((i + j < characteristicPoints.size()) &&
                   (almostEqual(birth(point),
                                birth(characteristicPoints[i + j]))) &&
                                  (death(point) <=
                                  death(characteristicPoints[i + j]))) {
              newCharacteristicPoints.push_back(characteristicPoints[i + j]);
              
              ++j;
            }
            
            newCharacteristicPoints.push_back(point);
            
            while ((i + j < characteristicPoints.size()) &&
                   (birth(point) <= birth(characteristicPoints[i + j])) &&
                   (death(point) >= death(characteristicPoints[i + j]))) {
              newCharacteristicPoints.push_back(characteristicPoints[i + j]);
              ++j;
            }
            
          } else {
            lambda_n.push_back(
              std::make_pair(death(lambda_n[lambda_n.size() - 1]), 0.));
            lambda_n.push_back(
              std::make_pair(birth(characteristicPoints[i]), 0.));
          }
          lambda_n.push_back(characteristicPoints[i]);
        } else {
          newCharacteristicPoints.push_back(characteristicPoints[i]);
        }
        i = i + j;
      }
      lambda_n.push_back(
        std::make_pair(death(lambda_n[lambda_n.size() - 1]), 0.));
      lambda_n.push_back(std::make_pair(INT_MAX, 0.));
      
      // CHANGE
      characteristicPoints = newCharacteristicPoints;
      // characteristicPoints.swap(newCharacteristicPoints);
      
      lambda_n.erase(std::unique(lambda_n.begin(), lambda_n.end()),
                     lambda_n.end());
      
      this->land.push_back(lambda_n);
    }
    
  } else {
    
    // In this case we will build a landscape on a grid.
    double gridDiameter = grid_diameter;
    // REVIEW: Why create `minMax` rather than use `min_x` and `max_x`? -JCB
    std::pair<double, double> minMax = std::make_pair(min_x, max_x);
    size_t numberOfBins =
      2 * ((minMax.second - minMax.first) / gridDiameter) + 1;
    
    // The first element of a pair `std::pair< double, std::vector<double> >`
    // is an x-value. The second element is a vector of values of landscapes.
    std::vector<std::pair<double, std::vector<double>>>
      criticalValuesOnPointsOfGrid(numberOfBins);
    
    // Filling up the bins:
    
    // Now, the idea is to iterate on `this->land[lambda-1]` and use only points
    // over there. The problem is at the very beginning, when there is nothing
    // in `this->land`. That is why over here, we make a fake `this->land[0]`.
    // It will be later deleted before moving on.
    std::vector<std::pair<double, double>> aa;
    double x = minMax.first;
    for (size_t i = 0; i != numberOfBins; ++i) {
      std::vector<double> v;
      std::pair<double, std::vector<double>> p = std::make_pair(x, v);
      criticalValuesOnPointsOfGrid[i] = p;
      aa.push_back(std::make_pair(x, 0.));
      x += 0.5 * gridDiameter;
    }
    
    // For every persistent interval, sample on the grid.
    for (size_t intervalNo = 0; intervalNo != pd.size(); ++intervalNo) {
      // size_t beginn = (size_t)(2*( pd[intervalNo].first-minMax.first
      // )/( gridDiameter ))+1;
      size_t beginn = 0;
      
      while (beginn < criticalValuesOnPointsOfGrid.size()) {
        if (fabs(criticalValuesOnPointsOfGrid[beginn].first >
                   pd[intervalNo].first) &&
                   fabs(criticalValuesOnPointsOfGrid[beginn].first <
                     pd[intervalNo].second)) {
          criticalValuesOnPointsOfGrid[beginn].second.push_back(
              std::min(fabs(criticalValuesOnPointsOfGrid[beginn].first -
                pd[intervalNo].first),
                fabs(criticalValuesOnPointsOfGrid[beginn].first -
                  pd[intervalNo].second)));
        } else
          criticalValuesOnPointsOfGrid[beginn].second.push_back(0.0);
        
        ++beginn;
      }
    }
    
    // Now, the basic structure is created. We need to translate it to a
    // persistence landscape data structure. To do so, first we need to sort all
    // the vectors in `criticalValuesOnPointsOfGrid[i].second`.
    size_t maxNonzeroLambda = 0;
    for (size_t i = 0; i != criticalValuesOnPointsOfGrid.size(); ++i) {
      std::sort(criticalValuesOnPointsOfGrid[i].second.begin(),
                criticalValuesOnPointsOfGrid[i].second.end(),
                std::greater<double>());
      if (criticalValuesOnPointsOfGrid[i].second.size() > maxNonzeroLambda) {
        maxNonzeroLambda = criticalValuesOnPointsOfGrid[i].second.size();
      }
    }
    
    // Initialize to zero
    this->land.resize(maxNonzeroLambda, aa);
    
    // Add values
    for (unsigned int i = 0; i < criticalValuesOnPointsOfGrid.size(); i++) {
      for (size_t lambda = 0;
           lambda < criticalValuesOnPointsOfGrid[i].second.size(); ++lambda) {
        this->land[lambda][i] =
          std::make_pair(criticalValuesOnPointsOfGrid[i].first,
                         criticalValuesOnPointsOfGrid[i].second[lambda]);
      }
    }
  }
}

PersistenceLandscape::PersistenceLandscape(
  const PersistenceLandscape &original) {
  // Rcpp::Rcerr << "Running copy constructor \n";
  std::vector<std::vector<std::pair<double, double>>> land(
      original.land.size());
  for (size_t i = 0; i != original.land.size(); ++i) {
    land[i].insert(land[i].end(),
                   original.land[i].begin(),
                   original.land[i].end());
  }
  // CHANGE
  // this->land = land;
  this->land.swap(land);
}

PersistenceLandscape PersistenceLandscape::
  operator=(const PersistenceLandscape &original) {
    std::vector<std::vector<std::pair<double, double>>> land(
        original.land.size());
    for (size_t i = 0; i != original.land.size(); ++i) {
      land[i].insert(land[i].end(),
                     original.land[i].begin(),
                     original.land[i].end());
    }
    // CHANGE
    // this->land = land;
    this->land.swap(land);
    return *this;
  }

// REVIEW: What is this doing? -JCB
PersistenceLandscape::PersistenceLandscape(
  std::vector<std::vector<std::pair<double, double>>>
  landscapePointsWithoutInfinities) {
  for (size_t level = 0;
       level != landscapePointsWithoutInfinities.size();
       ++level) {
    std::vector<std::pair<double, double>> v;
    // v.push_back(std::make_pair(INT_MIN, 0.));
    v.insert(v.end(), landscapePointsWithoutInfinities[level].begin(),
             landscapePointsWithoutInfinities[level].end());
    // v.push_back(std::make_pair(INT_MAX, 0.));
    this->land.push_back(v);
  }
}

// Exact and discrete representations

std::vector<std::vector<std::pair<double, double>>> exactLandscapeToDiscrete(
    PersistenceLandscape l,
    double min_x, double max_x,
    double dx) {
  
  std::vector<std::vector<std::pair<double, double>>> out;
  
  std::pair<double, double> currentPoint;
  
  for (unsigned i = 0; i < l.land.size(); i++) {
    
    auto level = l.land[i];
    std::vector<std::pair<double, double>> level_out;
    
    // Start at the level value at `min_x`.
    double x_buff = min_x;
    double y_buff = l.computeValueAtAGivenPoint(i, min_x);
    currentPoint = std::make_pair(x_buff, y_buff);
    level_out.push_back(currentPoint);
    
    // Iterate over finite critical points...
    for (int j = 1; j < level.size() - 1; j++) {
      
      // Skip to the critical point rightward of `currentPoint` for which the
      // next critical point is just leftward of `currentPoint + dx`.
      if (level[j].first <= x_buff || level[j + 1].first < x_buff + dx)
        continue;
      
      // If the next critical point is at least `dx` rightward, then increment
      // linearly to it; else, compute and increment to the next level value.
      if (level[j].first < x_buff + dx) {
        
        x_buff += dx;
        y_buff = l.computeValueAtAGivenPoint(i, x_buff + dx);
        currentPoint = std::make_pair(x_buff, y_buff);
        level_out.push_back(currentPoint);
        
      } else {
        std::pair<double, double> nextPoint = level[j];
        
        // If change in x, increment linearly from current to just before next.
        if (nextPoint.first != currentPoint.first) {
          double delta_x = nextPoint.first - currentPoint.first;
          double delta_y = nextPoint.second - currentPoint.second;
          double slope = delta_y / delta_x;
          
          int n_incr = std::floor((std::min(nextPoint.first, max_x) - x_buff) /
                                  dx);
          for (int k = 0; k < n_incr; k++) {
            x_buff += dx;
            y_buff += dx * slope;
            level_out.push_back(std::make_pair(x_buff, y_buff));
          }
        }
        
        currentPoint = std::make_pair(x_buff, y_buff);
      }
    }
    
    // If `max_x` has not been reached, then increment along zero y values.
    if (x_buff + dx < max_x + epsi) {
      int n_incr = std::floor((max_x + epsi - x_buff) / dx);
      y_buff = 0;
      for (int k = 0; k < n_incr; k++) {
        x_buff += dx;
        level_out.push_back(std::make_pair(x_buff, y_buff));
      }
    }
    
    out.push_back(level_out);
  }
  
  return out;
}

NumericVector PersistenceLandscape::toDiscrete() {
  if (exact) {
    return discretePersistenceLandscapeToR(
      exactLandscapeToDiscrete(*this, min_x, max_x, dx));
  } else {
    return discretePersistenceLandscapeToR(this->land);
  }
}

std::vector<NumericVector> PersistenceLandscape::toExact() {
  if (! exact) {
    stop("Can not convert a discrete PL to an exact PL.");
  } else {
    return exactPersistenceLandscapeToR(this->land);
  }
}

PersistenceLandscape PersistenceLandscape::discretize(
    double min_x, double max_x,
    double dx) {
  
  PersistenceLandscape out;
  
  if (exact) {
    out.land = exactLandscapeToDiscrete(*this, min_x, max_x, dx);
  } else {
    warning("Can not yet re-discretize a discrete PL.");
    out.land = this->land;
  }
  out.exact = false;
  out.min_x = min_x;
  out.max_x = max_x;
  out.dx = dx;
  
  return out;
}

// REVIEW: Check only that resolutions are equal and that endpoints lie on the
// same grid of the shared resolution. -JCB
bool checkPairOfDiscreteLandscapes(
    const PersistenceLandscape &l1,
    const PersistenceLandscape &l2) {
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
    const PersistenceLandscape &l1,
    const PersistenceLandscape &l2) {
  if (l1.min_x != l2.min_x || l1.max_x != l2.max_x)
    return false;
  return true;
}

// NB: This operation does not verify that the input is correcctly encoded; in
// particular, it does not check that the resolution is fixed throughout. -JCB
std::vector<std::vector<std::pair<double, double>>> expandDiscreteLandscape(
    PersistenceLandscape l,
    double min_x, double max_x,
    double dx) {
  
  std::vector<std::vector<std::pair<double, double>>> out;
  
  // reality check
  // if (fabs(l.land[0][1].first - l.land[0][0].first - dx) > epsi)
  if (! almostEqual(l.land[0][1].first, l.land[0][0].first + dx))
    stop("Resolutions do not agree.");
  
  // numbers of additional grid points
  double diffLeft = l.land[0][0].first - min_x;
  double diffRight = max_x - l.land[0][l.land[0].size() - 1].first;
  int gridDiffLeft = std::ceil(std::max(0., (diffLeft) / dx));
  int gridDiffRight = std::ceil(std::max(0., diffRight / dx));
  
  for (std::vector<std::pair<double, double>> level : l.land) {
    std::vector<std::pair<double, double>> level_out;
    
    // prefix
    for (int i = gridDiffLeft; i > 0; i--)
      level_out.push_back(std::make_pair(
          level[0].first - i * dx,
          0));
    // existing
    for (std::pair<double, double> pair : level)
      level_out.push_back(std::make_pair(
          pair.first,
          pair.second));
    // suffix
    for (int i = 1; i <= gridDiffRight; i++)
      level_out.push_back(std::make_pair(
          level[level.size() - 1].first + i * dx,
          0));
    
    out.push_back(level_out);
  }
  
  return out;
}

std::vector<std::vector<std::pair<double, double>>> contractDiscreteLandscape(
    PersistenceLandscape l,
    double min_x, double max_x,
    double dx) {
  
  std::vector<std::vector<std::pair<double, double>>> out;
  
  // reality check
  // if (fabs(l.land[0][1].first - l.land[0][0].first - dx) > epsi)
  if (! almostEqual(l.land[0][1].first, l.land[0][0].first + dx))
    stop("Resolutions do not agree.");
  
  // TODO: If no grid points lie between `min_x` and `max_x`, then return an
  // empty landscape.
  // `PersistenceLandscape l;`
  
  // numbers of fewer grid points
  double diffLeft = min_x - l.land[0][0].first;
  double diffRight = l.land[0][l.land[0].size() - 1].first - max_x;
  int gridDiffLeft = std::ceil(std::max(0., (diffLeft) / dx));
  int gridDiffRight = std::ceil(std::max(0., diffRight / dx));
  
  for (std::vector<std::pair<double, double>> level : l.land) {
    std::vector<std::pair<double, double>> level_out;
    
    // truncated range
    for (int i = gridDiffLeft; i < l.land[0].size() - gridDiffRight; i++)
      level_out.push_back(std::make_pair(
          level[i].first,
          level[i].second));
    
    out.push_back(level_out);
  }
  
  return out;
}

// Expands the limits of a PL -JCB
PersistenceLandscape PersistenceLandscape::expand(
    double min_x, double max_x) {
  
  PersistenceLandscape out;
  
  if (exact)
    out.land = this->land;
  else
    out.land = expandDiscreteLandscape(
      this->land,
      min_x, max_x,
      dx);
  out.exact = exact;
  out.min_x = min_x;
  out.max_x = max_x;
  out.dx = dx;
  
  return out;
}

// Contracts the limits of a PL -JCB
PersistenceLandscape PersistenceLandscape::contract(
    double min_x, double max_x) {
  
  PersistenceLandscape out;
  
  if (exact)
    out.land = this->land;
  else
    out.land = contractDiscreteLandscape(
      this->land,
      min_x, max_x,
      dx);
  out.exact = exact;
  out.min_x = min_x;
  out.max_x = max_x;
  out.dx = dx;
  
  return out;
}

// this is O(log(n)) algorithm, where n is number of points in this->land.
double PersistenceLandscape::computeValueAtAGivenPoint(
    unsigned level,
    double x) const {
  // in such a case lambda_level = 0.
  if (level > this->land.size())
    return 0;
  
  // we know that the points in this->land[level] are ordered according to x
  // coordinate. Therefore, we can find the point by using bisection:
  unsigned coordBegin = 1;
  unsigned coordEnd = this->land[level].size() - 2;
  
  // in this case x is outside the support of the landscape, therefore the value
  // of the landscape is 0.
  if (x <= this->land[level][coordBegin].first)
    return 0;
  if (x >= this->land[level][coordEnd].first)
    return 0;
  
  while (coordBegin + 1 != coordEnd) {
    unsigned newCord = (unsigned)floor((coordEnd + coordBegin) / 2.0);
    
    if (this->land[level][newCord].first <= x) {
      coordBegin = newCord;
      if (this->land[level][newCord].first == x)
        return this->land[level][newCord].second;
    } else {
      coordEnd = newCord;
    }
  }
  
  return functionValue(this->land[level][coordBegin],
                       this->land[level][coordEnd], x);
}

bool PersistenceLandscape::operator==(const PersistenceLandscape &rhs) const {
  if (this->land.size() != rhs.land.size()) {
    return false;
  }
  for (size_t level = 0; level != this->land.size(); ++level) {
    if (this->land[level].size() != rhs.land[level].size()) {
      return false;
    }
    for (size_t i = 0; i != this->land[level].size(); ++i) {
      if (this->land[level][i] != rhs.land[level][i]) {
        return false;
      }
    }
  }
  return true;
}

// This function finds the minimum value at a level.
double PersistenceLandscape::minimum(
    unsigned level) const {
  // if (this->land.size() < level)
  //   return 0;
  level = level - 1;
  if (level < 0 || level >= this->land.size())
    return NA_REAL;
  
  double minimum = INT_MAX;
  for (size_t i = 0; i != this->land[level].size(); ++i) {
    if (this->land[level][i].second < minimum)
      minimum = this->land[level][i].second;
  }
  return minimum;
}

// This function finds the maximum value at a level.
double PersistenceLandscape::maximum(
    unsigned level) const {
  // if (this->land.size() < level)
  //   return 0;
  level = level - 1;
  if (level < 0 || level >= this->land.size())
    return NA_REAL;
  
  double maximum = INT_MIN;
  for (size_t i = 0; i != this->land[level].size(); ++i) {
    if (this->land[level][i].second > maximum)
      maximum = this->land[level][i].second;
  }
  return maximum;
}

// REVIEW: Assume only that the `dx` are equal and that they divide the
// difference between the `min_x`. -JCB
std::vector<std::vector<std::pair<double, double>>>
  operationOnDiscreteLandscapes(
    const PersistenceLandscape &l1,
    const PersistenceLandscape &l2,
    double (*oper)(double, double)) {
    
    int min_level = std::min(l1.land.size(), l2.land.size());
    std::vector<std::vector<std::pair<double, double>>> out;
    for (int i = 0; i < min_level; i++) {
      int min_index = std::min(l1.land[i].size(), l2.land[i].size());
      std::vector<std::pair<double, double>> level_out;
      for (int j = 0; j < min_index; j++)
        level_out.push_back(std::make_pair(
            l1.land[i][j].first,
            // l1.land[i][j].second + l2.land[i][j].second));
            oper(l1.land[i][j].second, l2.land[i][j].second)));
      
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

// REVIEW: Assume only that the `dx` are equal and that they divide the
// difference between the `min_x`. -JCB
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

std::vector<std::vector<std::pair<double, double>>> scaleDiscreteLandscape(
    PersistenceLandscape l,
    double scale) {
  
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

// This function computes the n^th moment of a level.
double PersistenceLandscape::moment(
    unsigned p,
    double center,
    unsigned level) const {
  if (p < 1)
    throw("Cannot compute p^th moment for p < 1.\n");
  level = level - 1;
  if (level < 0 || level >= this->land.size())
    return NA_REAL;
  
  double result = 0;
  if (this->land.size() > level) {
    for (size_t i = 2; i != this->land[level].size() - 1; ++i) {
      if (this->land[level][i].first - this->land[level][i - 1].first == 0)
        continue;
      // Between `this->land[level][i]` and `this->land[level][i-1]`, the
      // `lambda_level` is of the form a x + b. First we need to find a and b.
      double a =
        (this->land[level][i].second - this->land[level][i - 1].second) /
          (this->land[level][i].first - this->land[level][i - 1].first);
      double b =
        this->land[level][i - 1].second - a * this->land[level][i - 1].first;
      
      double x1 = this->land[level][i - 1].first;
      double x2 = this->land[level][i].first;
      
      // double first =
      // b*(pow((x2-center),(double)(p+1))/(p+1)-
      // pow((x1-center),(double)(p+1))/(p+1));
      // double second = a/(p+1)*((x2*pow((x2-center),(double)(p+1))) -
      // (x1*pow((x1-center),(double)(p+1))) )
      //              +
      //              a/(p+1)*( pow((x2-center),(double)(p+2))/(p+2) -
      //              pow((x1-center),(double)(p+2))/(p+2) );
      // result += first;
      // result += second;
      
      double first = a / (p + 2) *
        (pow((x2 - center), (double)(p + 2)) -
        pow((x1 - center), (double)(p + 2)));
      double second = center / (p + 1) *
        (pow((x2 - center), (double)(p + 1)) -
        pow((x1 - center), (double)(p + 1)));
      double third = b / (p + 1) *
        (pow((x2 - center), (double)(p + 1)) -
        pow((x1 - center), (double)(p + 1)));
      
      result += first + second + third;
    }
  }
  return result;
}

double PersistenceLandscape::computeIntegralOfLandscape() const {
  double result = 0;
  
  for (size_t i = 0; i != this->land.size(); ++i) {
    // REVIEW: Handle exact and discrete cases differently. -JCB
    int infs = this->land[i][0].first == INT_MIN;
    // It suffices to compute every planar integral and then sum them up for
    // each `lambda_n`.
    // for (size_t nr = 2; nr != this->land[i].size() - 1; ++nr) {
    for (size_t nr = 1 + infs; nr != this->land[i].size() - infs; ++nr) {
      result += 0.5 * (this->land[i][nr].first - this->land[i][nr - 1].first) *
        (this->land[i][nr].second + this->land[i][nr - 1].second);
    }
  }
  
  return result;
}

std::pair<double, double> computeParametersOfALine(
    std::pair<double, double> p1,
    std::pair<double, double> p2) {
  // p1.second = a * p1.first + b => b = p1.second - a * p1.first
  // p2.second = a * p2.first + b = a * p2.first + p1.second - a * p1.first =
  // p1.second + a * ( p2.first - p1.first )
  // =>
  // ( p2.second - p1.second ) / ( p2.first - p1.first ) = a
  // b = p1.second - a * p1.first
  double a = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - a * p1.first;
  return std::make_pair(a, b);
}

double PersistenceLandscape::computeIntegralOfLandscape(
    double p) const {
  double result = 0;
  
  for (size_t i = 0; i != this->land.size(); ++i) {
    // REVIEW: Handle exact and discrete cases differently. -JCB
    int infs = this->land[i][0].first == INT_MIN;
    // for (size_t nr = 2; nr != this->land[i].size() - 1; ++nr) {
    for (size_t nr = 1 + infs; nr != this->land[i].size() - infs; ++nr) {
      // In this interval, the landscape has a form f(x) = ax + b. We want to
      // compute integral of (ax + b)^p = 1 / a * (ax + b)^{p + 1} / (p + 1)
      std::pair<double, double> coef =
        computeParametersOfALine(this->land[i][nr], this->land[i][nr - 1]);
      double a = coef.first;
      // double b = coef.second;
      
      if (this->land[i][nr].first == this->land[i][nr - 1].first)
        continue;
      // REVIEW: Debug discrepancy with R implementation. -JCB
      // if (a != 0) {
      // if (fabs(a) > epsi) {
      if (! almostEqual(a, 0.)) {
        // REVIEW: Simplify this formula:
        // result += 1 / (a * (p + 1)) *
        //   (pow((a * this->land[i][nr].first + b), p + 1) -
        //   pow((a * this->land[i][nr - 1].first + b), p + 1));
        result += 1 / (a * (p + 1)) *
          (pow(this->land[i][nr].second, p + 1) -
          pow(this->land[i][nr - 1].second, p + 1));
      } else {
        result += (this->land[i][nr].first - this->land[i][nr - 1].first) *
          (pow(this->land[i][nr].second, p));
      }
    }
  }
  
  return result;
}

// The `indicator` function is a vector of pairs. Its length is the number of
// levels on which it may be nonzero. See Section 3.6 of Bubenik (2015).
PersistenceLandscape PersistenceLandscape::multiplyByIndicatorFunction(
    std::vector<std::pair<double, double>> indicator,
    unsigned r) const {
  
  PersistenceLandscape result;
  
  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    
    double lev_c = pow(pow(lev + 1, -1), r);
    std::vector<std::pair<double, double>> lambda_n;
    
    // left (lower) limit
    if (exact) {
      lambda_n.push_back(std::make_pair(INT_MIN, 0.));
    }
    
    // if the indicator has at least `lev` levels...
    if (indicator.size() > lev) {
      
      if (exact) {
        // original method, for exact landscapes
        
        // loop over the critical points...
        for (size_t nr = 0; nr != this->land[lev].size(); ++nr) {
          
          // critical point lies before left endpoint; exclude
          if (this->land[lev][nr].first < indicator[lev].first) {
            continue;
          }
          
          // critical point lies just after right endpoint; interpolate
          if (this->land[lev][nr].first > indicator[lev].second) {
            lambda_n.push_back(std::make_pair(
                indicator[lev].second,
                functionValue(this->land[lev][nr - 1], this->land[lev][nr],
                              indicator[lev].second) * lev_c));
            lambda_n.push_back(std::make_pair(indicator[lev].second, 0.));
            break;
          }
          
          // critical point lies just after left endpoint; interpolate
          if ((this->land[lev][nr].first >= indicator[lev].first) &&
              (this->land[lev][nr - 1].first <= indicator[lev].first)) {
            lambda_n.push_back(std::make_pair(indicator[lev].first, 0.));
            lambda_n.push_back(std::make_pair(
                indicator[lev].first,
                functionValue(this->land[lev][nr - 1], this->land[lev][nr],
                              indicator[lev].first) * lev_c));
          }
          
          // critical point lies between left and right endpoints; include
          // lambda_n.push_back(this->land[lev][nr]);
          lambda_n.push_back(std::make_pair(
              this->land[lev][nr].first,
              this->land[lev][nr].second * lev_c));
        }
        
      } else {
        // method for discrete landscapes
        
        // loop over grid...
        for (size_t nr = 0; nr != this->land[lev].size(); ++nr) {
          
          if (this->land[lev][nr].first >= indicator[lev].first &&
              this->land[lev][nr].first <= indicator[lev].second) {
            // critical point lies inside endpoints; include
            
            // lambda_n.push_back(this->land[lev][nr]);
            lambda_n.push_back(std::make_pair(
                this->land[lev][nr].first,
                this->land[lev][nr].second * lev_c));
          } else {
            // critical point lies outside endpoints; exclude
            
            lambda_n.push_back(std::make_pair(this->land[lev][nr].first, 0.));
          }
        }
        
      }
      
    }
    
    // right (upper) limit
    if (exact) {
      lambda_n.push_back(std::make_pair(INT_MAX, 0.));
    }
    
    // cases with no critical points
    if (lambda_n.size() > 2) {
      result.land.push_back(lambda_n);
    }
  }
  return result;
}

// REVIEW: May want to make this only a single function.
PersistenceLandscape PersistenceLandscape::multiplyByIndicatorFunction(
    List indicator,
    unsigned r) const {
  
  PersistenceLandscape result;
  
  // Encode the list of vectors as a vector of pairs.
  std::vector<std::pair<double, double>> ind;
  for (size_t i = 0; i != indicator.length(); ++i) {
    std::vector<double> supp = indicator[i];
    ind.push_back(std::make_pair(supp[0], supp[1]));
  }
  
  return this->multiplyByIndicatorFunction(ind, r);
}

double
  PersistenceLandscape::computeIntegralOfLandscapeMultipliedByIndicatorFunction(
    std::vector<std::pair<double, double>> indicator,
    unsigned r,
    // This function computes the integral of the p^th power of a landscape.
    double p) const {
    
    PersistenceLandscape l = this->multiplyByIndicatorFunction(indicator, r);
    
    double result;
    if (p == 1)
      result = l.computeIntegralOfLandscape();
    else
      result = l.computeIntegralOfLandscape(p);
    
    return result;
  }

double
  PersistenceLandscape::computeIntegralOfLandscapeMultipliedByIndicatorFunction(
    List indicator,
    unsigned r,
    double p) const {
  
  // Encode the list of vectors as a vector of pairs.
  std::vector<std::pair<double, double>> ind;
  for (size_t i = 0; i != indicator.length(); ++i) {
    std::vector<double> supp = indicator[i];
    ind.push_back(std::make_pair(supp[0], supp[1]));
  }
  
  double result;
  result =
    this->computeIntegralOfLandscapeMultipliedByIndicatorFunction(ind, r, p);
  
  return result;
}

// This is a standard function which pairs maxima and minima which are not more
// than epsilon apart. This algorithm does not reduce all of them, just makes
// one pass through data. In order to reduce all of them, use the function
// `reduceAllPairsOfLowPersistenceMaximaMinima( double epsilon )`. WARNING! THIS
// PROCEDURE MODIFIES THE LANDSCAPE!!!
unsigned PersistenceLandscape::removePairsOfLocalMaximumMinimumOfEpsPersistence(
    double epsilon) {
  unsigned numberOfReducedPairs = 0;
  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    // Make sure that the loop in below is not infinite.
    if (2 > this->land[lev].size() - 3)
      continue;
    for (size_t nr = 2; nr != this->land[lev].size() - 3; ++nr) {
      if ((fabs(this->land[lev][nr].second - this->land[lev][nr + 1].second) <
        epsilon) &&
        (this->land[lev][nr].second != this->land[lev][nr + 1].second)) {
        // Right now we modify only the values of a points. That means that
        // slopes of lines in the landscape change a bit. This is the easiest
        // computational way of doing this. But I am not sure if this is the
        // best way of doing such a reduction of nonessential critical points.
        // Think about this!
        if (this->land[lev][nr].second < this->land[lev][nr + 1].second) {
          this->land[lev][nr].second = this->land[lev][nr + 1].second;
        } else {
          this->land[lev][nr + 1].second = this->land[lev][nr].second;
        }
        ++numberOfReducedPairs;
      }
    }
  }
  return numberOfReducedPairs;
}

// This procedure reduces all critical points of low persistence.
void PersistenceLandscape::reduceAllPairsOfLowPersistenceMaximaMinima(
    double epsilon) {
  unsigned numberOfReducedPoints = 1;
  while (numberOfReducedPoints) {
    numberOfReducedPoints =
      this->removePairsOfLocalMaximumMinimumOfEpsPersistence(epsilon);
  }
}

// It may happened that some landscape points obtained as the result of an
// algorithm lie on a line. In this case, the following procedure allows to
// remove unnecessary points.
void PersistenceLandscape::reduceAlignedPoints(
    // This parameter says how much the intercept and slope may be different to
    // consider points aligned.
    double tol) {
  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    size_t nr = 1;
    std::vector<std::pair<double, double>> lambda_n;
    lambda_n.push_back(this->land[lev][0]);
    while (nr != this->land[lev].size() - 2) {
      // First, compute the intercept and slope of a line crossing
      // `this->land[lev][nr]` and `this->land[lev][nr+1]`.
      std::pair<double, double> res = computeParametersOfALine(
        this->land[lev][nr], this->land[lev][nr + 1]);
      lambda_n.push_back(this->land[lev][nr]);
      
      double a = res.first;
      double b = res.second;
      int i = 1;
      while (nr + i != this->land[lev].size() - 2) {
        std::pair<double, double> res1 = computeParametersOfALine(
          this->land[lev][nr], this->land[lev][nr + i + 1]);
        if ((fabs(res1.first - a) < tol) &&
            (fabs(res1.second - b) < tol)) {
          ++i;
        } else {
          break;
        }
      }
      nr += i;
    }
    lambda_n.push_back(this->land[lev][this->land[lev].size() - 2]);
    lambda_n.push_back(this->land[lev][this->land[lev].size() - 1]);
    
    // If something was reduced, then replace `this->land[lev]` with the new
    // `lambda_n`.
    if (lambda_n.size() < this->land[lev].size()) {
      if (lambda_n.size() > 4) {
        this->land[lev].swap(lambda_n);
      }
      /*else
       {
       this->land[lev].clear();
       }*/
    }
  }
}

// Yet another function to smooth up the data. The idea of this one is as
// follows. Let us take a landscape point A which is not (+infty,0), (-infty,0)
// of (a,0), (b,0), where a and b denotes the points which support of the
// function begins and ends. Let B and C will be the landscape points after A.
// Suppose B and C are also no one as above. The question we are asking here is
// -- can we remove the point B and draw a line from A to C such that the
// difference in a landscape will be not greater than epsilon? To measure the
// penalty of removing B, the funcion penalty. In below, the simplese example is
// given:

double penalty(
    std::pair<double, double> A,
    std::pair<double, double> B,
    std::pair<double, double> C) {
  return fabs(functionValue(A, C, B.first) - B.second);
}

unsigned PersistenceLandscape::reducePoints(
    double tol,
    double (*penalty)(std::pair<double, double>, std::pair<double, double>,
            std::pair<double, double>)) {
  unsigned numberOfPointsReduced = 0;
  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    size_t nr = 1;
    std::vector<std::pair<double, double>> lambda_n;
    lambda_n.push_back(this->land[lev][0]);
    while (nr <= this->land[lev].size() - 2) {
      lambda_n.push_back(this->land[lev][nr]);
      if (penalty(this->land[lev][nr], this->land[lev][nr + 1],
                  this->land[lev][nr + 2]) < tol) {
        ++nr;
        ++numberOfPointsReduced;
      }
      ++nr;
    }
    lambda_n.push_back(this->land[lev][this->land[lev].size() - 2]);
    lambda_n.push_back(this->land[lev][this->land[lev].size() - 1]);
    
    // if something was reduced, then replace this->land[lev] with the new
    // lambda_n.
    if (lambda_n.size() < this->land[lev].size()) {
      if (lambda_n.size() > 4) {
        // CHANGE
        // this->land[lev] = lambda_n;
        this->land[lev].swap(lambda_n);
      } else {
        this->land[lev].clear();
      }
    }
  }
  return numberOfPointsReduced;
}

double findZeroOfALineSegmentBetweenThoseTwoPoints(
    std::pair<double, double> p1,
    std::pair<double, double> p2) {
  if (p1.first == p2.first)
    return p1.first;
  if (p1.second * p2.second > 0) {
    std::ostringstream errMessage;
    errMessage << "In function findZeroOfALineSegmentBetweenThoseTwoPoints the "
    "agguments are: ("
    << p1.first << "," << p1.second << ") and (" << p2.first << ","
    << p2.second
    << "). There is no zero in line between those two points. "
    "Program terminated.";
    std::string errMessageStr = errMessage.str();
    const char *err = errMessageStr.c_str();
    throw(err);
  }
  // we assume here, that x \in [ p1.first, p2.first ] and p1 and p2 are points
  // between which we will put the line segment
  double a = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - a * p1.first;
  // cerr << "Line crossing points : (" << p1.first << "," << p1.second << ")
  // oraz (" << p2.first << "," << p2.second << ") : \n"; cerr << "a : " << a <<
  // " , b : " << b << " , x : " << x << endl;
  return -b / a;
}

PersistenceLandscape PersistenceLandscape::abs() {
  PersistenceLandscape result;
  for (size_t level = 0; level != this->land.size(); ++level) {
    std::vector<std::pair<double, double>> lambda_n;
    // REVIEW: Try to prevent operations from infinitizing endpoints. -JCB
    // if (this->land[level][0].first == INT_MIN) {
    //   lambda_n.push_back(std::make_pair(INT_MIN, 0.));
    // } else {
    //   lambda_n.push_back(std::make_pair(this->land[level][0].first,
    //                                     fabs(this->land[level][0].second)));
    // }
    lambda_n.push_back(std::make_pair(this->land[level][0].first,
                                      fabs(this->land[level][0].second)));
    for (size_t i = 1; i != this->land[level].size(); ++i) {
      // if a line segment between this->land[level][i-1] and
      // this->land[level][i] crosses the x-axis, then we have to add one
      // landscape point to result
      if ((this->land[level][i - 1].second) *
          (this->land[level][i].second) < 0) {
        double zero = findZeroOfALineSegmentBetweenThoseTwoPoints(
          this->land[level][i - 1], this->land[level][i]);
        
        lambda_n.push_back(std::make_pair(zero, 0.));
        lambda_n.push_back(std::make_pair(this->land[level][i].first,
                                          fabs(this->land[level][i].second)));
      } else {
        lambda_n.push_back(std::make_pair(this->land[level][i].first,
                                          fabs(this->land[level][i].second)));
      }
    }
    result.land.push_back(lambda_n);
  }
  return result;
}

PersistenceLandscape scaleExactLandscape(
    const PersistenceLandscape &pl,
    double x) {
  std::vector<std::vector<std::pair<double, double>>> result(pl.land.size());
  for (size_t lev = 0; lev != pl.land.size(); ++lev) {
    std::vector<std::pair<double, double>> lambda_lev(pl.land[lev].size());
    for (size_t i = 0; i != pl.land[lev].size(); ++i) {
      lambda_lev[i] = std::make_pair(pl.land[lev][i].first,
                                     x * pl.land[lev][i].second);
    }
    result[lev] = lambda_lev;
  }
  PersistenceLandscape res;
  // CHANGE
  // res.land = result;
  res.land.swap(result);
  return res;
}

PersistenceLandscape PersistenceLandscape::scalePersistenceLandscape(
    double scale) const {
  
  PersistenceLandscape result;
  
  if (exact)
    result = scaleExactLandscape(*this, scale);
  else
    result = PersistenceLandscape(scaleDiscreteLandscape(*this, scale));
  
  // REVIEW: Is this necessary?
  result.exact = this->exact;
  result.min_x = this->min_x;
  result.max_x = this->max_x;
  result.dx = this->dx;
  
  return result;
}

// Add or subtract exact landscapes.
PersistenceLandscape operationOnTwoExactLandscapes(
    const PersistenceLandscape &land1,
    const PersistenceLandscape &land2,
    double (*oper)(double, double)) {
  
  PersistenceLandscape result;
  std::vector<std::vector<std::pair<double, double>>> land(
      std::max(land1.land.size(), land2.land.size()));
  result.land = land;
  
  for (size_t i = 0; i != std::min(land1.land.size(), land2.land.size()); ++i) {
    std::vector<std::pair<double, double>> lambda_n;
    int p = 0;
    int q = 0;
    while ((p + 1 < land1.land[i].size()) && (q + 1 < land2.land[i].size())) {
      
      if (land1.land[i][p].first < land2.land[i][q].first) {
        lambda_n.push_back(std::make_pair(
            land1.land[i][p].first,
            oper(land1.land[i][p].second,
                 functionValue(land2.land[i][q - 1], land2.land[i][q],
                               land1.land[i][p].first))));
        ++p;
        continue;
      }
      if (land1.land[i][p].first > land2.land[i][q].first) {
        lambda_n.push_back(std::make_pair(
            land2.land[i][q].first,
            oper(functionValue(land1.land[i][p], land1.land[i][p - 1],
                               land2.land[i][q].first),
                               land2.land[i][q].second)));
        ++q;
        continue;
      }
      if (land1.land[i][p].first == land2.land[i][q].first) {
        lambda_n.push_back(std::make_pair(
            land2.land[i][q].first,
            oper(land1.land[i][p].second, land2.land[i][q].second)));
        ++p;
        ++q;
      }
    }
    while ((p + 1 < land1.land[i].size()) && (q + 1 >= land2.land[i].size())) {
      lambda_n.push_back(std::make_pair(land1.land[i][p].first,
                                        oper(land1.land[i][p].second, 0)));
      ++p;
    }
    while ((p + 1 >= land1.land[i].size()) && (q + 1 < land2.land[i].size())) {
      lambda_n.push_back(std::make_pair(land2.land[i][q].first,
                                        oper(0, land2.land[i][q].second)));
      ++q;
    }
    // REVIEW: Try to prevent operations from infinitizing endpoints. -JCB
    if (land1.land[i][p].first == land2.land[i][q].first) {
      if (land2.land[i][q].first == INT_MAX) {
        lambda_n.push_back(std::make_pair(INT_MAX, 0.));
      } else {
        lambda_n.push_back(std::make_pair(land1.land[i][p].first,
                                          oper(land1.land[i][p].second,
                                               land2.land[i][q].second)));
      }
    } else {
      // REVIEW: Is there a better option in this case? -JCB
      lambda_n.push_back(std::make_pair(INT_MAX, 0.));
    }
    // CHANGE
    // result.land[i] = lambda_n;
    result.land[i].swap(lambda_n);
  }
  if (land1.land.size() > std::min(land1.land.size(), land2.land.size())) {
    for (size_t i = std::min(land1.land.size(), land2.land.size());
         i != std::max(land1.land.size(), land2.land.size()); ++i) {
      std::vector<std::pair<double, double>> lambda_n(land1.land[i]);
      for (size_t nr = 0; nr != land1.land[i].size(); ++nr) {
        lambda_n[nr] = std::make_pair(land1.land[i][nr].first,
                                      oper(land1.land[i][nr].second, 0));
      }
      // CHANGE
      // result.land[i] = lambda_n;
      result.land[i].swap(lambda_n);
    }
  }
  if (land2.land.size() > std::min(land1.land.size(), land2.land.size())) {
    for (size_t i = std::min(land1.land.size(), land2.land.size());
         i != std::max(land1.land.size(), land2.land.size()); ++i) {
      std::vector<std::pair<double, double>> lambda_n(land2.land[i]);
      for (size_t nr = 0; nr != land2.land[i].size(); ++nr) {
        lambda_n[nr] = std::make_pair(land2.land[i][nr].first,
                                      oper(0, land2.land[i][nr].second));
      }
      // CHANGE
      // result.land[i] = lambda_n;
      result.land[i].swap(lambda_n);
    }
  }
  return result;
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
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2) {
  
  bool p1 = 0;
  bool p2 = 0;
  
  if (pl1.exact != pl2.exact) {
    if (pl1.exact == true)
      p1 = 1;
    else
      // REVIEW: This has been edited from {tdatools} by JCB.
      p2 = 1;
  }
  
  return std::make_pair(p1, p2);
}

// Add or subtract two PLs.
PersistenceLandscape operationOnTwoLandscapes(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2,
    double (*oper)(double, double)) {
  
  PersistenceLandscape result;
  // appropriate settings for output PL
  bool res_exact = pl1.exact & pl2.exact;
  double res_min = min(pl1.min_x, pl2.min_x);
  double res_max = max(pl1.max_x, pl2.max_x);
  
  // both landscapes are exact
  if (pl1.exact && pl2.exact)
    // result = pl1 + pl2;
    operationOnDiscreteLandscapes(pl1, pl2, oper);
  
  // neither landscape is exact
  else if (! pl1.exact && ! pl2.exact) {
    if (! checkPairOfDiscreteLandscapes(pl1, pl2)) {
      stop("Resolutions and limits are incompatible.");
    }
    // REVIEW: Allow more landscapes to be added. Take output limits to be
    // union of input limits. -JCB
    if (alignPairOfDiscreteLandscapes(pl1, pl2)) {
      result =
        PersistenceLandscape(addDiscreteLandscapes(pl1, pl2));
    } else {
      result =
        PersistenceLandscape(addDiscreteLandscapes(
            expandDiscreteLandscape(pl2, res_min, res_max, pl1.dx),
            expandDiscreteLandscape(pl2, res_min, res_max, pl1.dx)));
    }
  }
  
  else {
    // Conversions:
    std::pair<bool, bool> conversions =
      operationOnPairOfLanscapesConversion(pl1, pl2);
    
    // only first landscape is exact
    if (conversions.first == true) {
      auto conversion1 = exactLandscapeToDiscrete(
        pl1,
        pl2.min_x,
        pl2.max_x,
        pl2.dx);
      result = PersistenceLandscape(addDiscreteLandscapes(
        conversion1,
        pl2));
    }
    
    // only second landscape is exact
    else if (conversions.second == true) {
      auto conversion2 = exactLandscapeToDiscrete(
        pl2,
        pl1.min_x,
        pl1.max_x,
        pl1.dx);
      result = PersistenceLandscape(addDiscreteLandscapes(
        conversion2,
        pl1));
      // REVIEW: This is an alternative formulation that seems to still work
      // but not fix the bug.
      // result = PersistenceLandscape(addDiscreteLandscapes(
      //   pl1,
      //   conversion2));
    }
  }
  
  result.exact = res_exact;
  result.min_x = res_min;
  result.max_x = res_max;
  result.dx = pl1.dx;
  
  return result;
}

double computeMaximalDistanceNonSymmetric(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2,
    unsigned &nrOfLand, double &x,
    double &y1, double &y2) {
  // this distance is not symmetric. It compute ONLY distance between inflection
  // points of pl1 and pl2.
  double maxDist = 0;
  int minimalNumberOfLevels = std::min(pl1.land.size(), pl2.land.size());
  for (int level = 0; level != minimalNumberOfLevels; ++level) {
    int p2Count = 0;
    for (int i = 1; i != pl1.land[level].size() - 1;
    ++i) // w tym przypadku nie rozwarzam punktow w nieskocznosci
    {
      while (true) {
        if ((pl1.land[level][i].first >= pl2.land[level][p2Count].first) &&
            (pl1.land[level][i].first <= pl2.land[level][p2Count + 1].first))
          break;
        p2Count++;
      }
      double val = functionValue(pl2.land[level][p2Count],
                                 pl2.land[level][p2Count + 1],
                                                pl1.land[level][i].first) -
                                                  pl1.land[level][i].second;
      val = fabs(val);
      
      // Rcpp::Rcerr << "functionValue( pl2.land[level][p2Count] ,
      // pl2.land[level][p2Count+1] , pl1.land[level][i].first ) : " <<
      // functionValue( pl2.land[level][p2Count] , pl2.land[level][p2Count+1] ,
      // pl1.land[level][i].first ) << "\n"; Rcpp::Rcerr <<
      // "pl1.land[level][i].second : " << pl1.land[level][i].second << "\n";
      // Rcpp::Rcerr << "pl1.land[level][i].first :" << pl1.land[level][i].first
      // << "\n"; std::cin.ignore();
      
      if (maxDist <= val) {
        maxDist = val;
        nrOfLand = level;
        x = pl1.land[level][i].first;
        y1 = pl1.land[level][i].second;
        y2 = functionValue(pl2.land[level][p2Count],
                           pl2.land[level][p2Count + 1],
                                          pl1.land[level][i].first);
      }
    }
  }
  
  if (minimalNumberOfLevels < pl1.land.size()) {
    for (int level = minimalNumberOfLevels; level != pl1.land.size(); ++level) {
      for (int i = 0; i != pl1.land[level].size(); ++i) {
        if (maxDist < pl1.land[level][i].second) {
          maxDist = pl1.land[level][i].second;
          nrOfLand = level;
          x = pl1.land[level][i].first;
          y1 = pl1.land[level][i].second;
          y2 = 0;
        }
      }
    }
  }
  return maxDist;
}

double computeInfNormDistanceBetweenLandscapes(
    const PersistenceLandscape &first,
    const PersistenceLandscape &second,
    unsigned &nrOfLand, double &x,
    double &y1, double &y2) {
  unsigned nrOfLandFirst;
  double xFirst, y1First, y2First;
  double dFirst = computeMaximalDistanceNonSymmetric(
    first, second, nrOfLandFirst, xFirst, y1First, y2First);
  
  unsigned nrOfLandSecond;
  double xSecond, y1Second, y2Second;
  double dSecond = computeMaximalDistanceNonSymmetric(
    second, first, nrOfLandSecond, xSecond, y1Second, y2Second);
  
  if (dFirst > dSecond) {
    nrOfLand = nrOfLandFirst;
    x = xFirst;
    y1 = y1First;
    y2 = y2First;
  } else {
    nrOfLand = nrOfLandSecond;
    x = xSecond;
    // this twist in below is neccesary!
    y2 = y1Second;
    y1 = y2Second;
    // y1 = y1Second;
    // y2 = y2Second;
  }
  return std::max(dFirst, dSecond);
}

double computeMaximalDistanceNonSymmetric(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2) {
  // this distance is not symmetric. It compute ONLY distance between inflection
  // points of pl1 and pl2.
  double maxDist = 0;
  int minimalNumberOfLevels = std::min(pl1.land.size(), pl2.land.size());
  for (int level = 0; level != minimalNumberOfLevels; ++level) {
    int p2Count = 0;
    for (int i = 1; i != pl1.land[level].size() - 1;
    ++i) // w tym przypadku nie rozwarzam punktow w nieskocznosci
    {
      while (true) {
        if ((pl1.land[level][i].first >= pl2.land[level][p2Count].first) &&
            (pl1.land[level][i].first <= pl2.land[level][p2Count + 1].first))
          break;
        p2Count++;
      }
      double val = functionValue(pl2.land[level][p2Count],
                                 pl2.land[level][p2Count + 1],
                                                pl1.land[level][i].first) -
                                                  pl1.land[level][i].second;
      val = fabs(val);
      if (maxDist <= val)
        maxDist = val;
    }
  }
  
  if (minimalNumberOfLevels < pl1.land.size()) {
    for (int level = minimalNumberOfLevels; level != pl1.land.size(); ++level) {
      for (int i = 0; i != pl1.land[level].size(); ++i) {
        if (maxDist < pl1.land[level][i].second)
          maxDist = pl1.land[level][i].second;
      }
    }
  }
  return maxDist;
}

double computeFinNormDistanceBetweenLandscapes(
    const PersistenceLandscape &first,
    const PersistenceLandscape &second,
    unsigned p) {
  // This is what we want to compute:
  // ( \int_{- \infty}^{+\infty} | first - second |^p )^(1/p)
  // We will do it one step at a time:
  
  // first - second
  PersistenceLandscape diff = first - second;
  
  // | first - second |
  diff = diff.abs();
  
  // \int_{- \infty}^{+\infty} | first-second |^p
  double result;
  if (p == 1) {
    result = diff.computeIntegralOfLandscape();
  } else {
    result = diff.computeIntegralOfLandscape(p);
  }
  
  // ( \int_{- \infty}^{+\infty} | first - second |^p )^(1/p)
  return pow(result, 1 / (double)p);
}

double computeInfNormDistanceBetweenLandscapes(
    const PersistenceLandscape &first,
    const PersistenceLandscape &second) {
  return std::max(computeMaximalDistanceNonSymmetric(first, second),
                  computeMaximalDistanceNonSymmetric(second, first));
}

double computeDistanceBetweenLandscapes(
    const PersistenceLandscape &first,
    const PersistenceLandscape &second,
    unsigned p) {
  
  double dist_out;
  
  if (p == R_PosInf)
    dist_out = computeInfNormDistanceBetweenLandscapes(first, second);
  else
    dist_out = computeFinNormDistanceBetweenLandscapes(first, second, p);
  
  return dist_out;
}

bool comparePairsForMerging(
    std::pair<double, unsigned> first,
    std::pair<double, unsigned> second) {
  return (first.first < second.first);
}

double innerProductExactLandscapes(
    const PersistenceLandscape &l1,
    const PersistenceLandscape &l2) {
  double result = 0;
  
  for (size_t level = 0; level != std::min(l1.size(), l2.size()); ++level) {
    if (l1.land[level].size() * l2.land[level].size() == 0)
      continue;
    
    // endpoints of the interval on which we will compute the inner product of
    // two locally linear functions:
    double x1 = INT_MIN;
    double x2;
    if (l1.land[level][1].first < l2.land[level][1].first) {
      x2 = l1.land[level][1].first;
    } else {
      x2 = l2.land[level][1].first;
    }
    
    // iterators for the landscapes l1 and l2
    size_t l1It = 0;
    size_t l2It = 0;
    
    while ((l1It < l1.land[level].size() - 1) &&
           (l2It < l2.land[level].size() - 1)) {
      // compute the value of a inner product on a interval [x1,x2]
      
      double a, b, c, d;
      
      a = (l1.land[level][l1It + 1].second - l1.land[level][l1It].second) /
        (l1.land[level][l1It + 1].first - l1.land[level][l1It].first);
      b = l1.land[level][l1It].second - a * l1.land[level][l1It].first;
      c = (l2.land[level][l2It + 1].second - l2.land[level][l2It].second) /
        (l2.land[level][l2It + 1].first - l2.land[level][l2It].first);
      d = l2.land[level][l2It].second - c * l2.land[level][l2It].first;
      
      double contributionFromThisPart =
        (a * c * x2 * x2 * x2 / 3 + (a * d + b * c) * x2 * x2 / 2 +
        b * d * x2) -
        (a * c * x1 * x1 * x1 / 3 + (a * d + b * c) * x1 * x1 / 2 +
        b * d * x1);
      
      result += contributionFromThisPart;
      
      // we have two intervals in which functions are constant:
      //[l1.land[level][l1It].first , l1.land[level][l1It+1].first]
      // and
      //[l2.land[level][l2It].first , l2.land[level][l2It+1].first]
      // We also have an interval [x1,x2]. Since the intervals in the landscapes
      // cover the whole R, then it is clear that x2 is either
      // l1.land[level][l1It+1].first of l2.land[level][l2It+1].first or both.
      // Lets test it.
      if (x2 == l1.land[level][l1It + 1].first) {
        if (x2 == l2.land[level][l2It + 1].first) {
          // in this case, we increment both:
          ++l2It;
        }
        ++l1It;
      } else {
        // in this case we increment l2It
        ++l2It;
      }
      // Now, we shift x1 and x2:
      x1 = x2;
      if (l1.land[level][l1It + 1].first < l2.land[level][l2It + 1].first) {
        x2 = l1.land[level][l1It + 1].first;
      } else {
        x2 = l2.land[level][l2It + 1].first;
      }
    }
  }
  return result;
}

double computeInnerProduct(
    const PersistenceLandscape &l1,
    const PersistenceLandscape &l2) {
  
  double result;
  
  if (l1.exact)
    result = innerProductExactLandscapes(l1.land, l2.land);
  else
    result = innerProductDiscreteLandscapes(l1.land, l2.land, l1.dx);
  
  return result;
}

// Functions to be exposed to R:

PersistenceLandscape PLsum(List pl_list) {
  
  PersistenceLandscape sum_pl = as<PersistenceLandscape>(pl_list[0]);
  
  for (int i = 1; i < pl_list.size(); i++) {
    sum_pl = sum_pl + as<PersistenceLandscape>(pl_list[i]);
  }
  
  return sum_pl;
}

// List PLdiff(List pl_list) {
//   
//   List diff_pls;
//   
//   for (int i = 1; i < pl_list.size(); i++) {
//     PersistenceLandscape diff_i = as<PersistenceLandscape>(pl_list[i]) -
//       as<PersistenceLandscape>(pl_list[i - 1]);
//     diff_pls.push_back(diff_i);
//   }
//   
//   return diff_pls;
// }

PersistenceLandscape PLmean(List pl_list) {
  
  double num = pl_list.size();
  PersistenceLandscape avg_pl = PLsum(pl_list);
  
  return avg_pl/num;
}

NumericMatrix PLdist(List pl_list, unsigned p) {
  
  // empty matrix
  int n = pl_list.size();
  NumericMatrix dist_mat(n, n);
  
  // not assuming symmetric distance calculation
  for (int i = 0; i != n; i++) {
    PersistenceLandscape pl_i = as<PersistenceLandscape>(pl_list[i]);
    for (int j = 0; j != n; j++) {
      if (j == i) {
        // same landscape
        dist_mat(i, j) = 0.;
      } else {
        // different landscape
        PersistenceLandscape pl_j = as<PersistenceLandscape>(pl_list[j]);
        // dist_mat(i, j) = pl_i.distance(pl_j, p);
        dist_mat(i, j) = computeDistanceBetweenLandscapes(pl_i, pl_j, p);
      }
    }
  }
  
  return dist_mat;
}

double PLvar(List pl_list, unsigned p) {
  
  // average landscape
  PersistenceLandscape avg = PLmean(pl_list);
  
  // sum-squared distance
  double ssd = 0;
  
  for (size_t i = 0; i != pl_list.size(); ++i) {
    
    PersistenceLandscape pl_i = as<PersistenceLandscape>(pl_list[i]);
    
    double d = computeDistanceBetweenLandscapes(avg, pl_i, p);
    ssd += d * d;
  }
  
  // sample standard deviation
  double var_out = ssd / pl_list.size();
  return var_out;
}

double PLsd(List pl_list, unsigned p) {
  
  double sd_out = PLvar(pl_list, p);
  
  sd_out = sqrt(sd_out);
  return sd_out;
}

#endif
