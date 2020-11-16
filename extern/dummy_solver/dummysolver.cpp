//
// dummysolver.cpp
//
// (c)2020 Applied Scientific Research, Inc.
//         Mark J. Stock <markjstock@gmail.com>

#include "dummysolver.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdint>

namespace DummySolver {


//
// receive node locations and element indices for internal cells, walls, and open boundaries
//
void
Solver::init_d_(std::vector<double> _pts,
                std::vector<uint32_t> _cidx,
                std::vector<uint32_t> _widx,
                std::vector<uint32_t> _oidx) {
  return;
}


//
// return all solution nodes
//
std::vector<double>
Solver::getsolnpts_d_() {
  return std::vector<double>();
}


//
// return solution nodes closest to open boundaries
//
std::vector<double>
Solver::getopenpts_d_() {
  return std::vector<double>();
}


//
// set velocity BCs at open boundaries
//
void
Solver::setopenvels_d_(std::vector<double> _vels) {
  return;
}


//
// solve system to the given time
//
void
Solver::solveto_d_(const double _endtime) {
  return;
}


//
// return scalar vorticity at all solution nodes
//
std::vector<double>
Solver::getallvorts_d_() {
  return std::vector<double>();
}


} // end namespace
