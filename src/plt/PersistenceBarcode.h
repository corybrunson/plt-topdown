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

#ifndef PERISTENCEBARCODES_H
#define PERISTENCEBARCODES_H

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <Rcpp.h>

#include "Configure.h"

// code taken from http://ranger.uta.edu/~weems/NOTES5311/hungarian.c
//#include "HungarianC.h"

double computeDistanceOfPointsInPlane(std::pair<double, double> p1,
                                      std::pair<double, double> p2) {
  // Rcpp::Rcerr << "Computing distance of points : (" << p1.first << "," <<
  // p1.second << ") and (" << p2.first << "," << p2.second << ")\n"; Rcpp::Rcerr
  // << "Distance : " << sqrt( (p1.first-p2.first)*(p1.first-p2.first) +
  // (p1.second-p2.second)*(p1.second-p2.second) ) << "\n";
  return sqrt((p1.first - p2.first) * (p1.first - p2.first) +
              (p1.second - p2.second) * (p1.second - p2.second));
} // computeDistanceOfPointsInPlane

std::pair<double, double> projectionToDiagonal(std::pair<double, double> p) {
  return std::make_pair(0.5 * (p.first + p.second), 0.5 * (p.first + p.second));
}

// to write individual bar
std::ostream &operator<<(std::ostream &out, std::pair<double, double> p) {
  out << p.first << " " << p.second;
  return out;
}

// this class represent a barcodes in a given dimension.
class PersistenceBarcodes {
public:
  PersistenceBarcodes(){};
  PersistenceBarcodes(const PersistenceBarcodes &orgyginal);
  PersistenceBarcodes(const char *filename);
  PersistenceBarcodes(const char *filename, double begin, double step);
  PersistenceBarcodes(std::vector<std::pair<double, double>>);
  PersistenceBarcodes(const char *filename, unsigned dimensionOfBarcode);
  PersistenceBarcodes(std::vector<std::pair<double, double>>,
                      unsigned dimensionOfBarcode);
  PersistenceBarcodes operator=(const PersistenceBarcodes &rhs);

  void addPair(double b, double d) {
    this->barcodes.push_back(std::make_pair(b, d));
  }

  double computeLandscapeIntegralFromBarcodes();

  // friend std::pair< double , std::vector< std::pair<
  // std::pair<double,double>,std::pair<double,double> > > >
  // computeBottleneckDistance( const PersistenceBarcodes& first, const
  // PersistenceBarcodes& second , unsigned p );

  friend std::ostream &operator<<(std::ostream &out, PersistenceBarcodes &bar) {
    for (size_t i = 0; i != bar.barcodes.size(); ++i) {
      out << bar.barcodes[i].first << " " << bar.barcodes[i].second
          << std::endl;
    }
    return out;
  }
  size_t size() const { return this->barcodes.size(); }
  bool empty() const { return (this->barcodes.size() == 0); }
  unsigned dim() const { return this->dimensionOfBarcode; }

  // iterators
  typedef std::vector<std::pair<double, double>>::iterator bIterator;
  bIterator bBegin() { return this->barcodes.begin(); }
  bIterator bEnd() { return this->barcodes.end(); }

  friend class PersistenceLandscape;

  std::pair<double, double> minMax() const;

  void sort();
  bool compare(PersistenceBarcodes &b);

  // TODO uncomment
  // private:
  std::vector<std::pair<double, double>> barcodes;
  unsigned dimensionOfBarcode;
};

bool comparePairs(const std::pair<double, double> &f,
                  const std::pair<double, double> &s) {
  if (f.first < s.first)
    return true;
  if (f.first > s.first)
    return false;
  if (f.second < s.second)
    return true;
  return false;
}
void PersistenceBarcodes::sort() {
  std::sort(this->barcodes.begin(), this->barcodes.end(), comparePairs);
} // sort
bool PersistenceBarcodes::compare(PersistenceBarcodes &b) {
  bool dbg = true;
  if (dbg) {
    Rcpp::Rcerr << "this->barcodes.size() : " << this->barcodes.size()
              << std::endl;
    Rcpp::Rcerr << "b.barcodes.size() : " << b.barcodes.size() << std::endl;
  }
  if (this->barcodes.size() != b.barcodes.size())
    return false;
  this->sort();
  b.sort();
  for (size_t i = 0; i != this->barcodes.size(); ++i) {
    if (this->barcodes[i] != b.barcodes[i]) {
      Rcpp::Rcerr << "this->barcodes[" << i << "] = " << this->barcodes[i]
                << std::endl;
      Rcpp::Rcerr << "b.barcodes[" << i << "] = " << b.barcodes[i] << std::endl;
      getchar();
      return false;
    }
  }
  return true;
}

size_t minn(size_t f, size_t s) {
  if (f < s)
    return f;
  return s;
}

PersistenceBarcodes::PersistenceBarcodes(const PersistenceBarcodes &orgyginal) {
  this->dimensionOfBarcode = orgyginal.dimensionOfBarcode;
  this->barcodes.insert(this->barcodes.end(), orgyginal.barcodes.begin(),
                        orgyginal.barcodes.end());
}

PersistenceBarcodes PersistenceBarcodes::
operator=(const PersistenceBarcodes &rhs) {
  // Rcpp::Rcout << "Before : " << this->barcodes.size() << "\n";

  this->dimensionOfBarcode = rhs.dimensionOfBarcode;
  this->barcodes.clear();
  this->barcodes.insert(this->barcodes.begin(), rhs.barcodes.begin(),
                        rhs.barcodes.end());

  // Rcpp::Rcout << "after : " << this->barcodes.size() << "\n";
  // std::cin.ignore();

  return *this;
}

std::pair<double, double> PersistenceBarcodes::minMax() const {
  double bmin = INT_MAX;
  double bmax = INT_MIN;
  for (size_t i = 0; i != this->barcodes.size(); ++i) {
    if (this->barcodes[i].first < bmin)
      bmin = this->barcodes[i].first;
    if (this->barcodes[i].second > bmax)
      bmax = this->barcodes[i].second;
  }
  return std::make_pair(bmin, bmax);
} // minMax

double PersistenceBarcodes::computeLandscapeIntegralFromBarcodes() {
  double result = 0;
  for (size_t i = 0; i != this->barcodes.size(); ++i) {
    result += (this->barcodes[i].second - this->barcodes[i].first) *
              (this->barcodes[i].second - this->barcodes[i].first);
  }
  result *= 0.25;
  return result;
}

// TODO2 -- consider adding some instructions to remove anything that is not
// numeric from the input stream.
PersistenceBarcodes::PersistenceBarcodes(const char *filename) {
  bool dbg = false;
  // cerr << "PersistenceBarcodes::PersistenceBarcodes(const char* filename)
  // \n";
  this->dimensionOfBarcode = 0;
  std::ifstream read;
  read.open(filename);
  if (!read.good()) {
    std::ostringstream errMessage;
    cerr
        << "In constructor PersistenceBarcodes(const char* filename). Filename "
        << filename << " do not exist \n";
    errMessage
        << "In constructor PersistenceBarcodes(const char* filename). Filename "
        << filename << " do not exist \n";
    std::string errMessageStr = errMessage.str();
    const char *err = errMessageStr.c_str();
    throw(err);
  } else {
    if (dbg) {
      Rcpp::Rcerr << "Reading file : " << filename << std::endl;
      Rcpp::Rcerr << "areThereInfiniteIntervals : " << areThereInfiniteIntervals
                << std::endl;
    }

    while (read.good()) {
      double begin, end;
      std::string line;
      std::getline(read, line);

      if (line.size()) {
        stringstream s;
        s << line;
        s >> begin;
        s >> end;
      } else {
        break;
      }

      if (end < begin) {
        if (!areThereInfiniteIntervals) {
          double z = end;
          end = begin;
          begin = z;
        } else {
          // in this case there are infinite intervals, so we need to check if
          // end != infty
          if (end != infty) {
            double z = end;
            end = begin;
            begin = z;
          }
        }
      }

      if (dbg) {
        Rcpp::Rcerr << "Reading interval : " << begin << " " << end << std::endl;
      }

      // if ( (!areThereInfiniteIntervals) && (end != infty) )
      if (!areThereInfiniteIntervals) {
        if (begin != end) {
          // TODO -- if you want to cut the data, comment the line in below.
          this->barcodes.push_back(std::make_pair(begin, end));
          // TODO - This is a correction that allows cutting of all the bars
          // which born before whereToCut value.
          /*
          double whereToCut = 2;
          if ( end > whereToCut )
          {
              if ( begin > whereToCut )
              {
                  this->barcodes.push_back( std::make_pair( begin,end ) );
              }
              else
              {
                  //in this case begin <= whereToCut and end > whereToCut
                  this->barcodes.push_back( std::make_pair( whereToCut,end ) );
              }
          }
          */
        }
      } else {
        if (dbg) {
          Rcpp::Rcerr
              << "There are infinite intervals. The vaue of infinity is : "
              << infty << std::endl;
        }
        if (end != infty) {
          this->barcodes.push_back(std::make_pair(begin, end));
        } else {
          if (dbg) {
            Rcpp::Rcerr << "The endpoint is infinity \n";
          }
          // we have here infinite interval.
          if (shallInfiniteBarcodesBeIgnored) {
            this->barcodes.push_back(std::make_pair(begin, end));
          } else {
            this->barcodes.push_back(std::make_pair(begin, valueOfInfinity));
          }
        }
      }
    }
    read.close();
  }
}

PersistenceBarcodes::PersistenceBarcodes(const char *filename, double bbegin,
                                         double step) {
  extern double infty;
  this->dimensionOfBarcode = 0;
  std::ifstream read;
  read.open(filename);
  if (!read.good()) {
    std::ostringstream errMessage;
    errMessage
        << "In constructor PersistenceBarcodes(const char* filename). Filename "
        << filename << " do not exist \n";
    std::string errMessageStr = errMessage.str();
    const char *err = errMessageStr.c_str();
    throw(err);
  }
  while (read.good()) {
    double begin, end;
    read >> begin;
    read >> end;

    if (end != infty) {
      if (end < begin) {
        double z = end;
        end = begin;
        begin = z;
      }
      if (!read.good())
        break;
      if (begin != end) {
        this->barcodes.push_back(
            std::make_pair(bbegin + begin * step, bbegin + end * step));
      }
    }
    /*
   else
   {
       //we have here infinite interval:
       this->barcodes.push_back( std::make_pair( begin,INT_MAX ) );
   }*/
  }
  read.close();
}

PersistenceBarcodes::PersistenceBarcodes(
    std::vector<std::pair<double, double>> bars) {
  extern double infty;
  this->dimensionOfBarcode = 0;
  unsigned sizeOfBarcode = 0;
  for (size_t i = 0; i != bars.size(); ++i) {
    if (bars[i].second != infty) {
      ++sizeOfBarcode;
    }
    if (bars[i].second < bars[i].first) {
      double sec = bars[i].second;
      bars[i].second = bars[i].first;
      bars[i].first = sec;
    }
  }
  std::vector<std::pair<double, double>> barcodes(sizeOfBarcode);
  unsigned nr = 0;
  for (size_t i = 0; i != bars.size(); ++i) {
    if (bars[i].second != infty) {
      // this is a finite interval
      barcodes[nr] = std::make_pair(bars[i].first, bars[i].second);
      ++nr;
    }
    // to keep it all compact for now I am removing infinite intervals from
    // consideration.
    /*else
    {
        //this is infinite interval:
        barcodes[i] =  std::make_pair( bars[i].first , INT_MAX );
    }*/
  }
  // CHANGE
  // this->barcodes = barcodes;
  this->barcodes.swap(barcodes);
}

PersistenceBarcodes::PersistenceBarcodes(const char *filename,
                                         unsigned dimensionOfBarcode) {
  *this = PersistenceBarcodes(filename);
  this->dimensionOfBarcode = dimensionOfBarcode;
}

PersistenceBarcodes::PersistenceBarcodes(
    std::vector<std::pair<double, double>> vect, unsigned dimensionOfBarcode) {
  *this = PersistenceBarcodes(vect);
  this->dimensionOfBarcode = dimensionOfBarcode;
}

#endif
