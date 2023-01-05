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
using namespace std;

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

// Modified by Jose Bouza to accept exact/discrete constructor param.
class PersistenceLandscape {
  
public:
  
  PersistenceLandscape(){}
  
  PersistenceLandscape(
    const std::vector<std::pair<double, double>> &diagram,
    bool exact = true,
    double min_x = 0,
    double max_x = 1,
    double dx = 0.01);
  
  std::vector<std::vector<std::pair<double, double>>> land;
  
  size_t size() const { return this->land.size(); }
  
  PersistenceLandscape(const PersistenceLandscape &original);
  
  PersistenceLandscape operator=(const PersistenceLandscape &original);
  
  PersistenceLandscape(
    std::vector<std::vector<std::pair<double, double>>>
    landscapePointsWithoutInfinities);
  
  double computeValueAtAGivenPoint(
      unsigned level,
      double x) const;
  
  PersistenceLandscape multiplyLanscapeByRealNumber(double x) const;
  
  // integral of (the p^th power of) a landscape
  double computeIntegralOfLandscape() const;
  double computeIntegralOfLandscape(double p) const;
  
  PersistenceLandscape multiplyByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator,
      unsigned r) const;
  
  // integral of (the p^th power of) the product of a landscape with an
  // indicator function
  double computeIntegralOfLandscapeMultipliedByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator,
      unsigned r) const;
  double computeIntegralOfLandscapeMultipliedByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator,
      unsigned r,
      double p) const;
  
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
  friend PersistenceLandscape operationOnPairOfLandscapes(
      const PersistenceLandscape &land1,
      const PersistenceLandscape &land2,
      double (*oper)(double, double));
  
  friend PersistenceLandscape addTwoLandscapes(
      const PersistenceLandscape &land1,
      const PersistenceLandscape &land2) {
    return operationOnPairOfLandscapes(land1, land2, add);
  }
  
  friend PersistenceLandscape subtractTwoLandscapes(
      const PersistenceLandscape &land1,
      const PersistenceLandscape &land2) {
    return operationOnPairOfLandscapes(land1, land2, sub);
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
    return first.multiplyLanscapeByRealNumber(con);
  }
  
  friend PersistenceLandscape operator*(
      double con,
      const PersistenceLandscape &first) {
    return first.multiplyLanscapeByRealNumber(con);
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

  bool operator==(const PersistenceLandscape &rhs) const;

  double computeMaximum() const {
    double maxValue = 0;
    if (this->land.size()) {
      maxValue = -INT_MAX;
      for (size_t i = 0; i != this->land[0].size(); ++i) {
        if (this->land[0][i].second > maxValue)
          maxValue = this->land[0][i].second;
      }
    }
    return maxValue;
  }

  double computeNormOfLandscape(int i) {
    PersistenceLandscape l;
    if (i != -1) {
      return computeDistanceBetweenLandscapes(*this, l, i);
    } else {
      return computeMaxNormDistanceBetweenLandscapes(*this, l);
    }
  }

  double operator()(unsigned level, double x) const {
    return this->computeValueAtAGivenPoint(level, x);
  }
  
  friend double computeMaxNormDistanceBetweenLandscapes(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second);
  
  friend double computeMaxNormDistanceBetweenLandscapes(
      const PersistenceLandscape &first,
      const PersistenceLandscape &second,
      unsigned &nrOfLand,
      double &x,
      double &y1, double &y2);
  
  friend double computeDistanceBetweenLandscapes(
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

  double minimalNonzeroPoint(unsigned l) const {
    if (this->land.size() < l)
      return INT_MAX;
    return this->land[l][1].first;
  }

  double maximalNonzeroPoint(unsigned l) const {
    if (this->land.size() < l)
      return INT_MIN;
    return this->land[l][this->land[l].size() - 2].first;
  }

  PersistenceLandscape abs();

  double findMin(unsigned lambda) const;
  
  double findMax(unsigned lambda) const;
  
  friend double computeInnerProduct(
      const PersistenceLandscape &l1,
      const PersistenceLandscape &l2);
  
  double computeNthMoment(
      unsigned p,
      double center,
      unsigned level) const;
  
  // These are two functions to generate histograms of Betti numbers across the
  // filtration values.
  std::vector<std::pair<double, unsigned>>
  generateBettiNumbersHistogram() const;
  
private:
  
  bool exact;
  
};

// REVIEW: What is this doing? -JCB
PersistenceLandscape::PersistenceLandscape(
  std::vector<std::vector<std::pair<double, double>>>
  landscapePointsWithoutInfinities) {
  for (size_t level = 0; level != landscapePointsWithoutInfinities.size();
       ++level) {
    std::vector<std::pair<double, double>> v;
    // v.push_back(std::make_pair(INT_MIN,0));
    v.insert(v.end(), landscapePointsWithoutInfinities[level].begin(),
             landscapePointsWithoutInfinities[level].end());
    // v.push_back(std::make_pair(INT_MAX,0));
    this->land.push_back(v);
  }
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

// This function finds the maximum value at the level `lambda`.
double PersistenceLandscape::findMax(
    unsigned lambda) const {
  if (this->land.size() < lambda)
    return 0;
  double maximum = INT_MIN;
  for (size_t i = 0; i != this->land[lambda].size(); ++i) {
    if (this->land[lambda][i].second > maximum)
      maximum = this->land[lambda][i].second;
  }
  return maximum;
}

// This function finds the minimum value at the level `lambda`.
double PersistenceLandscape::findMin(
    unsigned lambda) const {
  if (this->land.size() < lambda)
    return 0;
  double minimum = INT_MAX;
  for (size_t i = 0; i != this->land[lambda].size(); ++i) {
    if (this->land[lambda][i].second < minimum)
      minimum = this->land[lambda][i].second;
  }
  return minimum;
}

// This function computes the n^th moment of the level `lambda`.
double PersistenceLandscape::computeNthMoment(
    unsigned p,
    double center,
    unsigned level) const {
  if (p < 1) {
    Rcpp::Rcerr << "Cannot compute p^th moment for  p = " << p
                << ". The program will now terminate \n";
    throw("Cannot compute p^th moment. The program will now terminate \n");
  }
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

// The `indicator` function is a vector of pairs. Its length is the number of
// envelopes on which it may be nonzero. See Section 3.6 of Bubenik (2015).
PersistenceLandscape PersistenceLandscape::multiplyByIndicatorFunction(
    std::vector<std::pair<double, double>> indicator,
    unsigned r) const {
  
  PersistenceLandscape result;
  
  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    
    double lev_c = pow(pow(lev + 1, -1), r);
    std::vector<std::pair<double, double>> lambda_n;
    
    // left (lower) limit
    if (exact) {
      lambda_n.push_back(std::make_pair(INT_MIN, 0));
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
            lambda_n.push_back(std::make_pair(indicator[lev].second, 0));
            break;
          }
          
          // critical point lies just after left endpoint; interpolate
          if ((this->land[lev][nr].first >= indicator[lev].first) &&
              (this->land[lev][nr - 1].first <= indicator[lev].first)) {
            lambda_n.push_back(std::make_pair(indicator[lev].first, 0));
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
            
            lambda_n.push_back(std::make_pair(this->land[lev][nr].first, 0));
          }
        }
        
      }
      
    }
    
    // right (upper) limit
    if (exact) {
      lambda_n.push_back(std::make_pair(INT_MAX, 0));
    }
    
    // cases with no critical points
    if (lambda_n.size() > 2) {
      result.land.push_back(lambda_n);
    }
  }
  return result;
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
    // TODO: Harmonize this step with extended persistence data. -JCB
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
      lambda_n.push_back(std::make_pair(INT_MIN, 0));
      lambda_n.push_back(std::make_pair(birth(characteristicPoints[0]), 0));
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
                                   2,
                               (death(lambda_n[lambda_n.size() - 1]) -
                                birth(characteristicPoints[i])) /
                                   2);
            lambda_n.push_back(point);

            while ((i + j < characteristicPoints.size()) &&
                   (almostEqual(birth(point),
                                birth(characteristicPoints[i + j]))) &&
                   (death(point) <= death(characteristicPoints[i + j]))) {
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
                std::make_pair(death(lambda_n[lambda_n.size() - 1]), 0));
            lambda_n.push_back(
                std::make_pair(birth(characteristicPoints[i]), 0));
          }
          lambda_n.push_back(characteristicPoints[i]);
        } else {
          newCharacteristicPoints.push_back(characteristicPoints[i]);
        }
        i = i + j;
      }
      lambda_n.push_back(
          std::make_pair(death(lambda_n[lambda_n.size() - 1]), 0));
      lambda_n.push_back(std::make_pair(INT_MAX, 0));

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
    // in `this->land`. That is why over here, we make a fate `this->land[0]`.
    // It will be later deleted before moving on.
    std::vector<std::pair<double, double>> aa;
    double x = minMax.first;
    for (size_t i = 0; i != numberOfBins; ++i) {
      std::vector<double> v;
      std::pair<double, std::vector<double>> p = std::make_pair(x, v);
      criticalValuesOnPointsOfGrid[i] = p;
      aa.push_back(std::make_pair(x, 0));
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
      if (fabs(a) > epsi) {
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

double
PersistenceLandscape::computeIntegralOfLandscapeMultipliedByIndicatorFunction(
    std::vector<std::pair<double, double>> indicator,
    unsigned r) const {
  PersistenceLandscape l = this->multiplyByIndicatorFunction(indicator, r);
  return l.computeIntegralOfLandscape();
}

double
PersistenceLandscape::computeIntegralOfLandscapeMultipliedByIndicatorFunction(
    std::vector<std::pair<double, double>> indicator,
    unsigned r,
    // This function computes the integral of the p^th power of a landscape.
    double p) const {
  PersistenceLandscape l = this->multiplyByIndicatorFunction(indicator, r);
  return l.computeIntegralOfLandscape(p);
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

PersistenceLandscape PersistenceLandscape::abs() {
  PersistenceLandscape result;
  for (size_t level = 0; level != this->land.size(); ++level) {
    std::vector<std::pair<double, double>> lambda_n;
    // REVIEW: Try to prevent operations from infinitizing endpoints. -JCB
    // if (this->land[level][0].first == INT_MIN) {
    //   lambda_n.push_back(std::make_pair(INT_MIN, 0));
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
        
        lambda_n.push_back(std::make_pair(zero, 0));
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

PersistenceLandscape PersistenceLandscape::multiplyLanscapeByRealNumber(
    double x) const {
  std::vector<std::vector<std::pair<double, double>>> result(this->land.size());
  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    std::vector<std::pair<double, double>> lambda_lev(this->land[lev].size());
    for (size_t i = 0; i != this->land[lev].size(); ++i) {
      lambda_lev[i] = std::make_pair(this->land[lev][i].first,
                                     x * this->land[lev][i].second);
    }
    result[lev] = lambda_lev;
  }
  PersistenceLandscape res;
  // CHANGE
  // res.land = result;
  res.land.swap(result);
  return res;
}

PersistenceLandscape operationOnPairOfLandscapes(
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
        lambda_n.push_back(std::make_pair(INT_MAX, 0));
      } else {
        lambda_n.push_back(std::make_pair(land1.land[i][p].first,
                                          oper(land1.land[i][p].second,
                                               land2.land[i][q].second)));
      }
    } else {
      // REVIEW: Is there a better option in this case? -JCB
      lambda_n.push_back(std::make_pair(INT_MAX, 0));
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
      double val = fabs(functionValue(pl2.land[level][p2Count],
                                      pl2.land[level][p2Count + 1],
                                      pl1.land[level][i].first) -
                        pl1.land[level][i].second);

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

double computeMaxNormDistanceBetweenLandscapes(
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
      double val = fabs(functionValue(pl2.land[level][p2Count],
                                      pl2.land[level][p2Count + 1],
                                      pl1.land[level][i].first) -
                        pl1.land[level][i].second);
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

double computeDistanceBetweenLandscapes(
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

double computeMaxNormDistanceBetweenLandscapes(
    const PersistenceLandscape &first,
    const PersistenceLandscape &second) {
  return std::max(computeMaximalDistanceNonSymmetric(first, second),
                  computeMaximalDistanceNonSymmetric(second, first));
}

bool comparePairsForMerging(
    std::pair<double, unsigned> first,
    std::pair<double, unsigned> second) {
  return (first.first < second.first);
}

std::vector<std::pair<double, unsigned>>
PersistenceLandscape::generateBettiNumbersHistogram() const {
  
  std::vector<std::pair<double, unsigned>> resultRaw;

  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    std::vector<std::pair<double, unsigned>> rangeOfLandscapeInThisDimension;
    if (lev > 0) {
      for (size_t i = 1; i != this->land[lev].size() - 1; ++i) {
        if (this->land[lev][i].second == 0) {
          rangeOfLandscapeInThisDimension.push_back(
              std::make_pair(this->land[lev][i].first, lev + 1));
        }
      }
    } else {
      // lev == 0.
      bool first = true;
      for (size_t i = 1; i != this->land[lev].size() - 1; ++i) {
        if (this->land[lev][i].second == 0) {
          if (first) {
            rangeOfLandscapeInThisDimension.push_back(
                std::make_pair(this->land[lev][i].first, 0));
          }
          rangeOfLandscapeInThisDimension.push_back(
              std::make_pair(this->land[lev][i].first, lev + 1));
          if (!first) {
            rangeOfLandscapeInThisDimension.push_back(
                std::make_pair(this->land[lev][i].first, 0));
          }
          first = !first;
        }
      }
    }
    std::vector<std::pair<double, unsigned>> resultRawNew(
        resultRaw.size() + rangeOfLandscapeInThisDimension.size());
    std::merge(resultRaw.begin(), resultRaw.end(),
               rangeOfLandscapeInThisDimension.begin(),
               rangeOfLandscapeInThisDimension.end(), resultRawNew.begin(),
               comparePairsForMerging);
    resultRaw.swap(resultRawNew);
  }

  // now we should make it into a step function by adding a points in the jumps:
  std::vector<std::pair<double, unsigned>> result;
  if (resultRaw.size() == 0)
    return result;
  for (size_t i = 1; i != resultRaw.size(); ++i) {
    result.push_back(resultRaw[i - 1]);
    if (resultRaw[i - 1].second <= resultRaw[i].second) {
      result.push_back(
          std::make_pair(resultRaw[i].first, resultRaw[i - 1].second));
    } else {
      result.push_back(
          std::make_pair(resultRaw[i - 1].first, resultRaw[i].second));
    }
  }
  result.erase(unique(result.begin(), result.end()), result.end());

  /*
      //cleaning for Cathy
      std::vector< std::pair< double , unsigned > > resultNew;
      size_t i = 0;
      while ( i != result.size() )
      {
          int j = 1;
          resultNew.push_back( std::make_pair(result[i].first , maxBetti) );
          unsigned maxBetti = result[i].second;
          while ( (i+j<=result.size() ) && (result[i].first ==
     result[i+j].first) )
          {
              if ( maxBetti < result[i+j].second ){maxBetti =
     result[i+j].second;}
              ++j;
          }
          //i += std::max(j,1);
          resultNew.push_back( std::make_pair(result[i].first , maxBetti) );
          i += j;
      }
      result.swap(resultNew);
  */
  std::vector<std::pair<double, unsigned>> resultNew;
  size_t i = 0;
  while (i != result.size()) {
    double x = result[i].first;
    double maxBetti = result[i].second;
    double minBetti = result[i].second;
    while ((i != result.size()) && (fabs(result[i].first - x) < 0.000001)) {
      if (maxBetti < result[i].second)
        maxBetti = result[i].second;
      if (minBetti > result[i].second)
        minBetti = result[i].second;
      ++i;
    }
    if (minBetti != maxBetti) {
      if ((resultNew.size() == 0) ||
          (resultNew[resultNew.size() - 1].second <= minBetti)) {
        // going up
        resultNew.push_back(std::make_pair(x, minBetti));
        resultNew.push_back(std::make_pair(x, maxBetti));
      } else {
        // going down
        resultNew.push_back(std::make_pair(x, maxBetti));
        resultNew.push_back(std::make_pair(x, minBetti));
      }
    } else {
      resultNew.push_back(std::make_pair(x, minBetti));
    }
  }
  result.swap(resultNew);

  return result;
}

double computeInnerProduct(
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

#endif
