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
  std::cout << "  DummySolver set re= " << _re << std::endl;
  reynolds = _re;
}

void
Solver::set_elemorder_d_(const uint8_t _eo) {
  std::cout << "  DummySolver set element order= " << _eo << std::endl;
  elem_order = _eo;
}

void
Solver::set_timeorder_d_(const uint8_t _to) {
  std::cout << "  DummySolver set time integration order= " << _to << std::endl;
  time_order = _to;
}

void
Solver::set_numsteps_d_(const uint32_t _ns) {
  std::cout << "  DummySolver set num substeps= " << _ns << std::endl;
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

  //
  // the solver receives geometry nodes and elements and boundaries - save them
  //

  // first the node locations - used for all elements
  assert(_pts.size() % 2 == 0 && "WARN (Solver::init_d_) input _pts is not even");
  nodes = _pts;
  N_nodes = nodes.size()/2;
  std::cout << "  received " << N_nodes << " nodes" << std::endl;

  // 2d volume elements, assume 4 nodes per element
  //const size_t nper = _cidx.size() / 4;
  assert(_cidx.size() % 4 == 0 && "WARN (Solver::init_d_) input _cidx not a multiple of 4");
  elems = _cidx;
  N_elements = elems.size()/4;
  std::cout << "  received " << N_elements << " 2d elements" << std::endl;

  // the boundaries
  assert(_widx.size() % 2 == 0 && "WARN (Solver::init_d_) input _widx not a multiple of 2");
  wbdry = _widx;
  std::cout << "  received " << wbdry.size()/2 << " 1d wall boundary elements" << std::endl;

  assert(_oidx.size() % 2 == 0 && "WARN (Solver::init_d_) input _oidx not a multiple of 2");
  obdry = _oidx;
  std::cout << "  received " << obdry.size()/2 << " 1d open boundary elements" << std::endl;

  //
  // Use that information to generate *solution* nodes and elements
  // for this dummy case, assume all are 0th order elements: one node at the center
  //

  // first the 2d volume elements
  // take each element and make a single solution node from it
  for (size_t i=0; i<N_elements; ++i) {
    snodes.push_back(0.25*(nodes[2*elems[4*i]]+nodes[2*elems[4*i+1]]+
                           nodes[2*elems[4*i+2]]+nodes[2*elems[4*i+3]]));
    snodes.push_back(0.25*(nodes[2*elems[4*i]+1]+nodes[2*elems[4*i+1]+1]+
                           nodes[2*elems[4*i+2]+1]+nodes[2*elems[4*i+3]+1]));
    selems.push_back((uint32_t)i);
  }
  N_snodes = snodes.size()/2;
  N_selements = selems.size();

  // identify which of the solution nodes are the ones nearest the open boundary
  // HACK - just pick some
  for (size_t i=0; i<N_elements; i+=10) {
    sopts.push_back((uint32_t)i);
  }

  return;
}


//
// return all solution nodes
//
std::vector<double>
Solver::getsolnpts_d_() {
  return snodes;
}


//
// return solution nodes closest to open boundaries
//
std::vector<double>
Solver::getopenpts_d_() {
  // must assemble vector of just those nodes
  std::vector<double> opts;
  for (size_t i=0; i<sopts.size(); i+=10) {
    opts.push_back(snodes[2*sopts[i]]);
    opts.push_back(snodes[2*sopts[i]+1]);
  }
  return opts;
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
