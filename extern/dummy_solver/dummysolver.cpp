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
#include <cmath>

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
  std::cout << "  DummySolver set element order= " << (int)_eo << std::endl;
  elem_order = _eo;
}

void
Solver::set_timeorder_d_(const uint8_t _to) {
  std::cout << "  DummySolver set time integration order= " << (int)_to << std::endl;
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

  std::cout << "  DummySolver initializing" << std::endl;

  //
  // the solver receives geometry nodes and elements and boundaries - save them
  //

  // first the node locations - used for all elements
  assert(_pts.size() % 2 == 0 && "WARN (Solver::init_d_) input _pts is not even");
  nodes = _pts;
  N_nodes = nodes.size()/2;
  std::cout << "    received " << N_nodes << " nodes" << std::endl;

  // 2d volume elements, assume 4 nodes per element
  //const size_t nper = _cidx.size() / 4;
  assert(_cidx.size() % 4 == 0 && "WARN (Solver::init_d_) input _cidx not a multiple of 4");
  elems = _cidx;
  N_elements = elems.size()/4;
  std::cout << "    received " << N_elements << " 2d elements" << std::endl;

  // the boundaries
  assert(_widx.size() % 2 == 0 && "WARN (Solver::init_d_) input _widx not a multiple of 2");
  wbdry = _widx;
  std::cout << "    received " << wbdry.size()/2 << " 1d wall boundary elements" << std::endl;

  assert(_oidx.size() % 2 == 0 && "WARN (Solver::init_d_) input _oidx not a multiple of 2");
  obdry = _oidx;
  std::cout << "    received " << obdry.size()/2 << " 1d open boundary elements" << std::endl;

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
  std::cout << "    generated " << N_snodes << " solution nodes" << std::endl;

  // identify which of the solution nodes are the ones nearest the open boundary
  // HACK - assume it's the nodes with rad > 0.95
  for (size_t i=0; i<selems.size(); ++i) {
    const double xp = snodes[2*selems[i]];
    const double yp = snodes[2*selems[i]+1];
    if (xp*xp+yp*yp > 0.9) sopts.push_back((uint32_t)i);
  }
  assert(obdry.size()/2 == sopts.size() && "ERROR (Solver::init_d_) hack does not appear to be working");
  std::cout << "    of which " << sopts.size() << " are on the open boundary" << std::endl;

  curr_time = 0.0;

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
  for (size_t i=0; i<sopts.size(); ++i) {
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
  std::cout << "  DummySolver set velocities at open soln nodes: " << _vels.size()/2 << std::endl;
  assert(_vels.size() == sopts.size()*2 && "ERROR (Solver::setopenvels_d_) bad incoming velocity vector length");

  if (false) {
    std::cout << "    x  y  ux  uy" << std::endl;
    for (size_t i=0; i<_vels.size()/2; ++i) {
      std::cout << "    " << snodes[2*sopts[i]] << " " << snodes[2*sopts[i]+1] << "  " << _vels[2*i] << " " << _vels[2*i+1] << std::endl;
    }
  }

  return;
}


//
// set initial vorticity at all solution nodes
//
void
Solver::setsolnvort_d_(std::vector<double> _vort) {
  std::cout << "  DummySolver set vorticity at all soln nodes: " << _vort.size() << std::endl;
  assert(_vort.size() == N_snodes && "ERROR (Solver::setsolnvort_d_) bad incoming vorticity vector length");

  if (false) {
    std::cout << "    x  y  vort" << std::endl;
    for (size_t i=0; i<_vort.size(); ++i) {
      std::cout << "    " << snodes[2*i] << " " << snodes[2*i+1] << "  " << _vort[i] << std::endl;
    }
  }

  return;
}


//
// solve system to the given time
//
void
Solver::solveto_d_(const double _endtime) {
  std::cout << "DummySolver solving to t= " << _endtime << " with " << num_substeps << " substeps" << std::endl;

  double this_dt = (_endtime - curr_time) / (double)num_substeps;
  for (size_t step=0; step<num_substeps; ++step) {
    std::cout << "  substep " << step << " at t= " << curr_time << std::endl;
    curr_time += this_dt;
  }
  std::cout << "  solver time is now " << curr_time << std::endl;

  return;
}


//
// return scalar vorticity at all solution nodes
//
std::vector<double>
Solver::getallvorts_d_() {
  std::cout << "  DummySolver returning vorticity at all soln nodes" << std::endl;

  // HACK - return a diffuse-like vorticity
  std::vector<double> vort;
  vort.resize(N_snodes);

  for (size_t i=0; i<vort.size(); ++i) {
    const double xp = snodes[2*i];
    const double yp = snodes[2*i+1];
    const double r = std::sqrt(xp*xp + yp*yp);
    // this is a model of an exponential-like bump from r=0.5 to ~0.75
    const double factor = 5.0 * ((r<0.74 && r>0.501) ? 0.5 - 0.5*std::cos(10.0*std::pow(r-0.5,0.333333)) : 0.0);
    // and this makes it negative vort on top and positive beneath
    vort[i] = factor * (-yp) / r;
  }

  return vort;
}


} // end namespace
