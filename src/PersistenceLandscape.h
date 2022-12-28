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

#pragma once

#ifndef PERISTENCELANDSCAPE_H
#define PERISTENCELANDSCAPE_H

// #include "Configure.h"
// #include "PersistenceBarcode.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <limits>
#include <list>
#include <unistd.h>
#include <vector>
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

// functions used in PersistenceLandscape( const PersistenceBarcodes& pb )
// constructor:
bool comparePoints(std::pair<double, double> f, std::pair<double, double> s) {
  double differenceBirth = birth(f) - birth(s);
  if (differenceBirth < 0)
    differenceBirth *= -1;
  double differenceDeath = death(f) - death(s);
  if (differenceDeath < 0)
    differenceDeath *= -1;

  if ((differenceBirth < epsi) && (differenceDeath < epsi)) {
    return false;
  }
  if ((differenceBirth < epsi)) {
    // consider birth points the same. If we are here, we know that death points
    // are NOT the same
    if (death(f) < death(s)) {
      return true;
    }
    return false;
  }
  if (differenceDeath < epsi) {
    // we consider death points the same and since we are here, the birth points
    // are not the same!
    if (birth(f) < birth(s)) {
      return false;
    }
    return true;
  }

  if (birth(f) > birth(s)) {
    return false;
  }
  if (birth(f) < birth(s)) {
    return true;
  }
  // if this is true, we assume that death(f)<=death(s) -- othervise I have had
  // a lot of roundoff problems here!
  if ((death(f) <= death(s))) {
    return false;
  }
  return true;
}

// this function assumes birth-death coords
bool comparePoints2(std::pair<double, double> f, std::pair<double, double> s) {
  if (f.first < s.first) {
    return true;
  } else { // f.first >= s.first
    if (f.first > s.first) {
      return false;
    } else { // f.first == s.first
      if (f.second > s.second) {
        return true;
      } else {
        return false;
      }
    }
  }
}

class vectorSpaceOfPersistenceLandscapes;

// functions used to add and subtract landscapes
inline double add(double x, double y) { return x + y; }
inline double sub(double x, double y) { return x - y; }

// function used in computeValueAtAGivenPoint
double functionValue(std::pair<double, double> p1, std::pair<double, double> p2,
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
  PersistenceLandscape(const std::vector<std::pair<double, double>> &diagram,
                       bool exact = true,
                       double min_x = 0, double max_x = 10, double dx = 0.01);
  PersistenceLandscape operator=(const PersistenceLandscape &org);
  PersistenceLandscape(const PersistenceLandscape &);
  
  PersistenceLandscape(std::vector<std::vector<std::pair<double, double>>>
                           landscapePointsWithoutInfinities);
  std::vector<std::vector<std::pair<double, double>>>
  gimmeProperLandscapePoints();

  double computeIntegralOfLandscape() const;
  double computeIntegralOfLandscape(double p)
      const; // this function compute integral of p-th power of landscape.
  double computeIntegralOfLandscapeMultipliedByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator) const;
  double computeIntegralOfLandscapeMultipliedByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator,
      double p) const; // this function compute integral of p-th power of
                       // landscape multiplied by the indicator function.
  PersistenceLandscape multiplyByIndicatorFunction(
      std::vector<std::pair<double, double>> indicator) const;

  unsigned
  removePairsOfLocalMaximumMinimumOfEpsPersistence(double errorTolerance);
  void reduceAllPairsOfLowPersistenceMaximaMinima(double epsilon);
  void reduceAlignedPoints(double tollerance = 0.000001);
  unsigned reducePoints(double tollerance,
                        double (*penalty)(std::pair<double, double>,
                                          std::pair<double, double>,
                                          std::pair<double, double>));
  double computeValueAtAGivenPoint(unsigned level, double x) const;
  // friend std::ostream &operator<<(std::ostream &out,
  //                                 PersistenceLandscape &land);

  typedef std::vector<std::pair<double, double>>::iterator lDimIterator;
  lDimIterator lDimBegin(unsigned dim) {
    if (dim > this->land.size())
      throw("Calling lDimIterator in a dimension higher that dimension of "
            "landscape");
    return this->land[dim].begin();
  }
  lDimIterator lDimEnd(unsigned dim) {
    if (dim > this->land.size())
      throw("Calling lDimIterator in a dimension higher that dimension of "
            "landscape");
    return this->land[dim].end();
  }

  PersistenceLandscape multiplyLanscapeByRealNumberNotOverwrite(double x) const;
  void multiplyLanscapeByRealNumberOverwrite(double x);

  void plot(char *filename, size_t from = -1, size_t to = -1,
            double xRangeBegin = -1, double xRangeEnd = -1,
            double yRangeBegin = -1, double yRangeEnd = -1);

  // Friendzone:

  // this is a general algorithm to perform linear operations on persisntece
  // lapscapes. It perform it by doing operations on landscape points.
  friend PersistenceLandscape
  operationOnPairOfLandscapes(const PersistenceLandscape &land1,
                              const PersistenceLandscape &land2,
                              double (*oper)(double, double));

  friend PersistenceLandscape
  addTwoLandscapes(const PersistenceLandscape &land1,
                   const PersistenceLandscape &land2) {
    return operationOnPairOfLandscapes(land1, land2, add);
  }
  friend PersistenceLandscape
  subtractTwoLandscapes(const PersistenceLandscape &land1,
                        const PersistenceLandscape &land2) {
    return operationOnPairOfLandscapes(land1, land2, sub);
  }

  friend PersistenceLandscape operator+(const PersistenceLandscape &first,
                                        const PersistenceLandscape &second) {
    return addTwoLandscapes(first, second);
  }

  friend PersistenceLandscape operator-(const PersistenceLandscape &first,
                                        const PersistenceLandscape &second) {
    return subtractTwoLandscapes(first, second);
  }

  friend PersistenceLandscape operator*(const PersistenceLandscape &first,
                                        double con) {
    return first.multiplyLanscapeByRealNumberNotOverwrite(con);
  }

  friend PersistenceLandscape operator*(double con,
                                        const PersistenceLandscape &first) {
    return first.multiplyLanscapeByRealNumberNotOverwrite(con);
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
      return computeDistanceOfLandscapes(*this, l, i);
    } else {
      return computeMaxNormDiscanceOfLandscapes(*this, l);
    }
  }

  double operator()(unsigned level, double x) const {
    return this->computeValueAtAGivenPoint(level, x);
  }

  friend double
  computeMaxNormDiscanceOfLandscapes(const PersistenceLandscape &first,
                                     const PersistenceLandscape &second);
  friend double computeMaxNormDiscanceOfLandscapes(
      const PersistenceLandscape &first, const PersistenceLandscape &second,
      unsigned &nrOfLand, double &x, double &y1, double &y2);

  friend double computeDistanceOfLandscapes(const PersistenceLandscape &first,
                                            const PersistenceLandscape &second,
                                            unsigned p);

  friend double
  computeMaximalDistanceNonSymmetric(const PersistenceLandscape &pl1,
                                     const PersistenceLandscape &pl2);

  friend double computeMaximalDistanceNonSymmetric(
      const PersistenceLandscape &pl1, const PersistenceLandscape &pl2,
      unsigned &nrOfLand, double &x, double &y1, double &y2);
  // this function additionally returns integer n and double x, y1, y2 such that
  // the maximal distance is obtained betwenn lambda_n's on a coordinate x such
  // that the value of the first landscape is y1, and the vale of the second
  // landscape is y2.

  friend class vectorSpaceOfPersistenceLandscapes;

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

  size_t size() const { return this->land.size(); }

  double findMax(unsigned lambda) const;

  friend double computeInnerProduct(const PersistenceLandscape &l1,
                                    const PersistenceLandscape &l2);

  // this function compute n-th moment of lambda_level
  double computeNthMoment(unsigned n, double center, unsigned level) const;

  // those are two new functions to generate histograms of Betti numbers across
  // the filtration values.
  std::vector<std::pair<double, unsigned>>
  generateBettiNumbersHistogram() const;
  std::vector<std::vector<std::pair<double, double>>> land;

private:
  bool exact;
};

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

std::vector<std::vector<std::pair<double, double>>>
PersistenceLandscape::gimmeProperLandscapePoints() {
  std::vector<std::vector<std::pair<double, double>>> result;
  for (size_t level = 0; level != this->land.size(); ++level) {
    std::vector<std::pair<double, double>> v(this->land[level].begin() + 1,
                                             this->land[level].end() - 1);
    result.push_back(v);
  }
  return result;
}

inline bool check_if_file_exist(const char *name) {
  return (access(name, F_OK) != -1);
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

// this function find maximum of lambda_n
double PersistenceLandscape::findMax(unsigned lambda) const {
  if (this->land.size() < lambda)
    return 0;
  double maximum = INT_MIN;
  for (size_t i = 0; i != this->land[lambda].size(); ++i) {
    if (this->land[lambda][i].second > maximum)
      maximum = this->land[lambda][i].second;
  }
  return maximum;
}

// this function compute n-th moment of lambda_level
double PersistenceLandscape::computeNthMoment(unsigned n, double center,
                                              unsigned level) const {
  if (n < 1) {
    cerr << "Cannot compute n-th moment for  n = " << n
         << ". The program will now terminate \n";
    throw("Cannot compute n-th moment. The program will now terminate \n");
  }
  double result = 0;
  if (this->land.size() > level) {
    for (size_t i = 2; i != this->land[level].size() - 1; ++i) {
      if (this->land[level][i].first - this->land[level][i - 1].first == 0)
        continue;
      // between this->land[level][i] and this->land[level][i-1] the
      // lambda_level is of the form ax+b. First we need to find a and b.
      double a =
          (this->land[level][i].second - this->land[level][i - 1].second) /
          (this->land[level][i].first - this->land[level][i - 1].first);
      double b =
          this->land[level][i - 1].second - a * this->land[level][i - 1].first;

      double x1 = this->land[level][i - 1].first;
      double x2 = this->land[level][i].first;

      // double first =
      // b*(pow((x2-center),(double)(n+1))/(n+1)-
      // pow((x1-center),(double)(n+1))/(n+1));
      // double second = a/(n+1)*((x2*pow((x2-center),(double)(n+1))) -
      // (x1*pow((x1-center),(double)(n+1))) )
      //              +
      //              a/(n+1)*( pow((x2-center),(double)(n+2))/(n+2) -
      //              pow((x1-center),(double)(n+2))/(n+2) );
      // result += first;
      // result += second;

      double first = a / (n + 2) *
                     (pow((x2 - center), (double)(n + 2)) -
                      pow((x1 - center), (double)(n + 2)));
      double second = center / (n + 1) *
                      (pow((x2 - center), (double)(n + 1)) -
                       pow((x1 - center), (double)(n + 1)));
      double third = b / (n + 1) *
                     (pow((x2 - center), (double)(n + 1)) -
                      pow((x1 - center), (double)(n + 1)));

      result += first + second + third;
    }
  }
  return result;
} // computeNthMoment

bool multiplyByIndicatorFunctionBDG = false;
PersistenceLandscape PersistenceLandscape::multiplyByIndicatorFunction(
    std::vector<std::pair<double, double>> indicator) const {
  PersistenceLandscape result;
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    if (multiplyByIndicatorFunctionBDG) {
      Rcpp::Rcout << "dim : " << dim << "\n";
    }
    std::vector<std::pair<double, double>> lambda_n;
    lambda_n.push_back(std::make_pair(0, INT_MIN));
    if (indicator.size() > dim) {
      if (multiplyByIndicatorFunctionBDG) {
        Rcpp::Rcout << "There is nonzero indicator in this dimension\n";
        Rcpp::Rcout << "[ " << indicator[dim].first << " , "
                  << indicator[dim].second << "] \n";
      }
      for (size_t nr = 0; nr != this->land[dim].size(); ++nr) {
        if (multiplyByIndicatorFunctionBDG) {
          Rcpp::Rcout << "this->land[dim][nr] : " << this->land[dim][nr].first
                    << " , " << this->land[dim][nr].second << "\n";
        }
        if (this->land[dim][nr].first < indicator[dim].first) {
          if (multiplyByIndicatorFunctionBDG) {
            Rcpp::Rcout << "Below treshold\n";
          }
          continue;
        }
        if (this->land[dim][nr].first > indicator[dim].second) {
          if (multiplyByIndicatorFunctionBDG) {
            Rcpp::Rcout << "Just pass above treshold \n";
          }
          lambda_n.push_back(std::make_pair(
              indicator[dim].second,
              functionValue(this->land[dim][nr - 1], this->land[dim][nr],
                            indicator[dim].second)));
          lambda_n.push_back(std::make_pair(indicator[dim].second, 0));
          break;
        }
        if ((this->land[dim][nr].first >= indicator[dim].first) &&
            (this->land[dim][nr - 1].first <= indicator[dim].first)) {
          if (multiplyByIndicatorFunctionBDG) {
            Rcpp::Rcout << "Entering the indicator \n";
          }
          lambda_n.push_back(std::make_pair(indicator[dim].first, 0));
          lambda_n.push_back(std::make_pair(
              indicator[dim].first,
              functionValue(this->land[dim][nr - 1], this->land[dim][nr],
                            indicator[dim].first)));
        }

        if (multiplyByIndicatorFunctionBDG) {
          Rcpp::Rcout << "We are here\n";
        }
        lambda_n.push_back(std::make_pair(this->land[dim][nr].first,
                                          this->land[dim][nr].second));
      }
    }
    lambda_n.push_back(std::make_pair(0, INT_MIN));
    if (lambda_n.size() > 2) {
      result.land.push_back(lambda_n);
    }
  }
  return result;
}

PersistenceLandscape::PersistenceLandscape(
    const PersistenceLandscape &oryginal) {
  // Rcpp::Rcerr << "Running copy constructor \n";
  std::vector<std::vector<std::pair<double, double>>> land(
      oryginal.land.size());
  for (size_t i = 0; i != oryginal.land.size(); ++i) {
    land[i].insert(land[i].end(), oryginal.land[i].begin(),
                   oryginal.land[i].end());
  }
  // CHANGE
  // this->land = land;
  this->land.swap(land);
}

PersistenceLandscape PersistenceLandscape::
operator=(const PersistenceLandscape &oryginal) {
  std::vector<std::vector<std::pair<double, double>>> land(
      oryginal.land.size());
  for (size_t i = 0; i != oryginal.land.size(); ++i) {
    land[i].insert(land[i].end(), oryginal.land[i].begin(),
                   oryginal.land[i].end());
  }
  // CHANGE
  // this->land = land;
  this->land.swap(land);
  return *this;
}

// TODO -- removewhen the problem is respved
bool check(unsigned i, std::vector<std::pair<double, double>> v) {
  if ((i >= v.size())) {
    Rcpp::Rcout << "you want to get to index : " << i
              << " while there are only  : " << v.size() << " indices \n";
    std::cin.ignore();
    return true;
  }
  return false;
}
// if ( check( , ) ){Rcpp::Rcerr << "OUT OF MEMORY \n";}

// // TEST: Attempt to have `PersistenceLandscape()` handle matrices directly.
// // This function is adapted from the constructor:
// // `PersistenceBarcodes(std::vector<std::pair<double, double>> bars)`
// std::vector<std::pair<double, double>> PDSort(
//     std::vector<std::pair<double, double>> pd) {
//   extern double infty;
//   unsigned nb = 0;
//   for (size_t i = 0; i != pd.size(); ++i) {
//     if (pd[i].second != infty) {
//       ++nb;
//     }
//     if (pd[i].second < pd[i].first) {
//       double sec = pd[i].second;
//       pd[i].second = pd[i].first;
//       pd[i].first = sec;
//     }
//   }
//   std::vector<std::pair<double, double>> pds(nb);
//   unsigned nr = 0;
//   for (size_t i = 0; i != pd.size(); ++i) {
//     if (pd[i].second != infty) {
//       // this is a finite interval
//       pds[nr] = std::make_pair(pd[i].first, pd[i].second);
//       ++nr;
//     }
//   }
//   return pds;
// }

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
    // TODO: Harmonize this step with extended persistence data.
    if (pd[i].second < pd[i].first) {
      double sec = pd[i].second;
      pd[i].second = pd[i].first;
      pd[i].first = sec;
    }
  }
  unsigned nr = 0;
  for (size_t i = 0; i != pd.size(); ++i) {
    if (pd[i].second != R_PosInf & pd[i].second != R_NegInf) {
      // this is a finite interval
      pd[nr] = std::make_pair(pd[i].first, pd[i].second);
      ++nr;
    }
  }
  
  if (exact) {
    // this is a general algorithm to construct persistence landscapes.
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
    // in this case useGridInComputations is true, therefore we will build a
    // landscape on a grid.
    double gridDiameter = grid_diameter;
    // REVIEW: Why create `minMax` rather than use `min_x` and `max_x`? -JCB
    std::pair<double, double> minMax = std::make_pair(min_x, max_x);
    size_t numberOfBins =
        2 * ((minMax.second - minMax.first) / gridDiameter) + 1;

    // first element of a pair std::pair< double , std::vector<double> > is a
    // x-value. Second element is a vector of values of landscapes.
    std::vector<std::pair<double, std::vector<double>>>
        criticalValuesOnPointsOfGrid(numberOfBins);
    // filling up the bins:

    // Now, the idea is to iterate on this->land[lambda-1] and use only points
    // over there. The problem is at the very beginning, when there is nothing
    // in this->land. That is why over here, we make a fate this->land[0]. It
    // will be later deteted before moving on.
    std::vector<std::pair<double, double>> aa;
    double x = minMax.first;
    for (size_t i = 0; i != numberOfBins; ++i) {
      std::vector<double> v;
      std::pair<double, std::vector<double>> p = std::make_pair(x, v);
      criticalValuesOnPointsOfGrid[i] = p;
      aa.push_back(std::make_pair(x, 0));
      x += 0.5 * gridDiameter;
    }

    // for every peristent interval, sample on grid.
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

    // now, the basic structure is created. We need to translate it to a
    // persistence landscape data structure. To do so, first we need to sort all
    // the vectors in criticalValuesOnPointsOfGrid[i].second
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
    // it suffices to compute every planar integral and then sum them ap for
    // each lambda_n
    // for (size_t nr = 2; nr != this->land[i].size() - 1; ++nr) {
    for (size_t nr = 1 + infs; nr != this->land[i].size() - infs; ++nr) {
      result += 0.5 * (this->land[i][nr].first - this->land[i][nr - 1].first) *
        (this->land[i][nr].second + this->land[i][nr - 1].second);
    }
  }
  return result;
}

std::pair<double, double>
computeParametersOfALine(std::pair<double, double> p1,
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

double PersistenceLandscape::computeIntegralOfLandscape(double p) const {
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
    std::vector<std::pair<double, double>> indicator) const {
  PersistenceLandscape l = this->multiplyByIndicatorFunction(indicator);
  return l.computeIntegralOfLandscape();
}

double
PersistenceLandscape::computeIntegralOfLandscapeMultipliedByIndicatorFunction(
    std::vector<std::pair<double, double>> indicator, double p)
    const // this function compute integral of p-th power of landscape.
{
  PersistenceLandscape l = this->multiplyByIndicatorFunction(indicator);
  return l.computeIntegralOfLandscape(p);
}

// This is a standard function which pairs maxima and minima which are not more
// than epsilon apart. This algorithm do not reduce all of them, just make one
// passage through data. In order to reduce all of them use the function
// reduceAllPairsOfLowPersistenceMaximaMinima( double epsilon ) WARNING! THIS
// PROCEDURE MODIFIES THE LANDSCAPE!!!
unsigned PersistenceLandscape::removePairsOfLocalMaximumMinimumOfEpsPersistence(
    double epsilon) {
  unsigned numberOfReducedPairs = 0;
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    if (2 > this->land[dim].size() - 3)
      continue; // to make sure that the loop in below is not infinite.
    for (size_t nr = 2; nr != this->land[dim].size() - 3; ++nr) {
      if ((fabs(this->land[dim][nr].second - this->land[dim][nr + 1].second) <
           epsilon) &&
          (this->land[dim][nr].second != this->land[dim][nr + 1].second)) {
        // right now we modify only the lalues of a points. That means that
        // angles of lines in the landscape changes a bit. This is the easiest
        // computational way of doing this. But I am not sure if this is the
        // best way of doing such a reduction of nonessential critical points.
        // Think about this!
        if (this->land[dim][nr].second < this->land[dim][nr + 1].second) {
          this->land[dim][nr].second = this->land[dim][nr + 1].second;
        } else {
          this->land[dim][nr + 1].second = this->land[dim][nr].second;
        }
        ++numberOfReducedPairs;
      }
    }
  }
  return numberOfReducedPairs;
}

// this procedure redue all critical points of low persistence.
void PersistenceLandscape::reduceAllPairsOfLowPersistenceMaximaMinima(
    double epsilon) {
  unsigned numberOfReducedPoints = 1;
  while (numberOfReducedPoints) {
    numberOfReducedPoints =
        this->removePairsOfLocalMaximumMinimumOfEpsPersistence(epsilon);
  }
}

// It may happened that some landscape points obtained as a aresult of an
// algorithm lies in a line. In this case, the following procedure allows to
// remove unnecesary points.
bool reduceAlignedPointsBDG = false;
void PersistenceLandscape::reduceAlignedPoints(
    double tollerance) // this parapeter says how much the coeficients a and b
                       // in a formula y=ax+b may be different to consider
                       // points aligned.
{
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    size_t nr = 1;
    std::vector<std::pair<double, double>> lambda_n;
    lambda_n.push_back(this->land[dim][0]);
    while (nr != this->land[dim].size() - 2) {
      // first, compute a and b in formula y=ax+b of a line crossing
      // this->land[dim][nr] and this->land[dim][nr+1].
      std::pair<double, double> res = computeParametersOfALine(
          this->land[dim][nr], this->land[dim][nr + 1]);
      // if (reduceAlignedPointsBDG) {
      //   Rcpp::Rcout << "Considering points : "
      //             << this->land[dim][nr] << " and "
      //             << this->land[dim][nr + 1] << std::endl;
      //   Rcpp::Rcout << "Adding : " << this->land[dim][nr] << " to lambda_n."
      //             << std::endl;
      // }
      lambda_n.push_back(this->land[dim][nr]);

      double a = res.first;
      double b = res.second;
      int i = 1;
      while (nr + i != this->land[dim].size() - 2) {
        // if (reduceAlignedPointsBDG) {
        //   Rcpp::Rcout << "Checking if : " << this->land[dim][nr + i + 1]
        //             << " is aligned with them " << std::endl;
        // }
        std::pair<double, double> res1 = computeParametersOfALine(
            this->land[dim][nr], this->land[dim][nr + i + 1]);
        if ((fabs(res1.first - a) < tollerance) &&
            (fabs(res1.second - b) < tollerance)) {
          if (reduceAlignedPointsBDG) {
            Rcpp::Rcout << "It is aligned " << std::endl;
          }
          ++i;
        } else {
          if (reduceAlignedPointsBDG) {
            Rcpp::Rcout << "It is NOT aligned " << std::endl;
          }
          break;
        }
      }
      if (reduceAlignedPointsBDG) {
        Rcpp::Rcout << "We are out of the while loop. The number of aligned "
                     "points is : "
                  << i << std::endl; // std::cin.ignore();
      }
      nr += i;
    }
    // if (reduceAlignedPointsBDG) {
    //   Rcpp::Rcout << "Out  of main while loop, done with this dimension "
    //             << std::endl;
    //   Rcpp::Rcout << "Adding : "
    //             << this->land[dim][this->land[dim].size() - 2]
    //             << " to lamnda_n " << std::endl;
    //   Rcpp::Rcout << "Adding : "
    //             << this->land[dim][this->land[dim].size() - 1]
    //             << " to lamnda_n " << std::endl;
    //   std::cin.ignore();
    // }
    lambda_n.push_back(this->land[dim][this->land[dim].size() - 2]);
    lambda_n.push_back(this->land[dim][this->land[dim].size() - 1]);

    // if something was reduced, then replace this->land[dim] with the new
    // lambda_n.
    if (lambda_n.size() < this->land[dim].size()) {
      if (lambda_n.size() > 4) {
        this->land[dim].swap(lambda_n);
      }
      /*else
      {
          this->land[dim].clear();
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

double penalty(std::pair<double, double> A, std::pair<double, double> B,
               std::pair<double, double> C) {
  return fabs(functionValue(A, C, B.first) - B.second);
} // penalty

unsigned PersistenceLandscape::reducePoints(
    double tollerance,
    double (*penalty)(std::pair<double, double>, std::pair<double, double>,
                      std::pair<double, double>)) {
  unsigned numberOfPointsReduced = 0;
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    size_t nr = 1;
    std::vector<std::pair<double, double>> lambda_n;
    lambda_n.push_back(this->land[dim][0]);
    while (nr <= this->land[dim].size() - 2) {
      lambda_n.push_back(this->land[dim][nr]);
      if (penalty(this->land[dim][nr], this->land[dim][nr + 1],
                  this->land[dim][nr + 2]) < tollerance) {
        ++nr;
        ++numberOfPointsReduced;
      }
      ++nr;
    }
    lambda_n.push_back(this->land[dim][this->land[dim].size() - 2]);
    lambda_n.push_back(this->land[dim][this->land[dim].size() - 1]);

    // if something was reduced, then replace this->land[dim] with the new
    // lambda_n.
    if (lambda_n.size() < this->land[dim].size()) {
      if (lambda_n.size() > 4) {
        // CHANGE
        // this->land[dim] = lambda_n;
        this->land[dim].swap(lambda_n);
      } else {
        this->land[dim].clear();
      }
    }
  }
  return numberOfPointsReduced;
}

double
findZeroOfALineSegmentBetweenThoseTwoPoints(std::pair<double, double> p1,
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
double PersistenceLandscape::computeValueAtAGivenPoint(unsigned level,
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

void PersistenceLandscape::multiplyLanscapeByRealNumberOverwrite(double x) {
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    for (size_t i = 0; i != this->land[dim].size(); ++i) {
      this->land[dim][i].second *= x;
    }
  }
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

PersistenceLandscape
PersistenceLandscape::multiplyLanscapeByRealNumberNotOverwrite(double x) const {
  std::vector<std::vector<std::pair<double, double>>> result(this->land.size());
  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    std::vector<std::pair<double, double>> lambda_dim(this->land[dim].size());
    for (size_t i = 0; i != this->land[dim].size(); ++i) {
      lambda_dim[i] = std::make_pair(this->land[dim][i].first,
                                     x * this->land[dim][i].second);
    }
    result[dim] = lambda_dim;
  }
  PersistenceLandscape res;
  // CHANGE
  // res.land = result;
  res.land.swap(result);
  return res;
} // multiplyLanscapeByRealNumberOverwrite

PersistenceLandscape
operationOnPairOfLandscapes(const PersistenceLandscape &land1,
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
} // operationOnPairOfLandscapes

double computeMaximalDistanceNonSymmetric(const PersistenceLandscape &pl1,
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

double computeMaxNormDiscanceOfLandscapes(const PersistenceLandscape &first,
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

double computeMaximalDistanceNonSymmetric(const PersistenceLandscape &pl1,
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

double computeDistanceOfLandscapes(const PersistenceLandscape &first,
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
  if (p != 1) {
    result = diff.computeIntegralOfLandscape(p);
  } else {
    result = diff.computeIntegralOfLandscape();
  }

  // ( \int_{- \infty}^{+\infty} | first - second |^p )^(1/p)
  return pow(result, 1 / (double)p);
}

double computeMaxNormDiscanceOfLandscapes(const PersistenceLandscape &first,
                                          const PersistenceLandscape &second) {
  return std::max(computeMaximalDistanceNonSymmetric(first, second),
                  computeMaximalDistanceNonSymmetric(second, first));
}

bool comparePairsForMerging(std::pair<double, unsigned> first,
                            std::pair<double, unsigned> second) {
  return (first.first < second.first);
}

std::vector<std::pair<double, unsigned>>
PersistenceLandscape::generateBettiNumbersHistogram() const {
  
  std::vector<std::pair<double, unsigned>> resultRaw;

  for (size_t dim = 0; dim != this->land.size(); ++dim) {
    std::vector<std::pair<double, unsigned>> rangeOfLandscapeInThisDimension;
    if (dim > 0) {
      for (size_t i = 1; i != this->land[dim].size() - 1; ++i) {
        if (this->land[dim][i].second == 0) {
          rangeOfLandscapeInThisDimension.push_back(
              std::make_pair(this->land[dim][i].first, dim + 1));
        }
      }
    } else {
      // dim == 0.
      bool first = true;
      for (size_t i = 1; i != this->land[dim].size() - 1; ++i) {
        if (this->land[dim][i].second == 0) {
          if (first) {
            rangeOfLandscapeInThisDimension.push_back(
                std::make_pair(this->land[dim][i].first, 0));
          }
          rangeOfLandscapeInThisDimension.push_back(
              std::make_pair(this->land[dim][i].first, dim + 1));
          if (!first) {
            rangeOfLandscapeInThisDimension.push_back(
                std::make_pair(this->land[dim][i].first, 0));
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
} // generateBettiNumbersHistogram

double computeInnerProduct(const PersistenceLandscape &l1,
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
