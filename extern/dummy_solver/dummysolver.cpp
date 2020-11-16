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
// setters for HO solver runtime parameters
//
void
Solver::set_re_d_(const double _re) {
  std::cout << "DummySolver set re= " << _re << std::endl;
  reynolds = _re;
}

void
Solver::set_elemorder_d_(const uint8_t _eo) {
  std::cout << "DummySolver set element order= " << _eo << std::endl;
  elem_order = _eo;
}

void
Solver::set_timeorder_d_(const uint8_t _to) {
  std::cout << "DummySolver set time integration order= " << _to << std::endl;
  time_order = _to;
}

void
Solver::set_numsteps_d_(const uint32_t _ns) {
  std::cout << "DummySolver set num substeps= " << _ns << std::endl;
  num_substeps = _ns;
}


//
// receive node locations and element indices for internal cells, walls, and open boundaries
//
void
Solver::init_d_(std::vector<double> _pts,
                std::vector<uint32_t> _cidx,
                std::vector<uint32_t> _widx,
                std::vector<uint32_t> _oidx) {

  std::cout << "DummySolver initializing with " << _pts.size()/2 << " nodes" << std::endl;
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
  std::cout << "DummySolver solving to t= " << _endtime << " with " << num_substeps << " substeps" << std::endl;
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
