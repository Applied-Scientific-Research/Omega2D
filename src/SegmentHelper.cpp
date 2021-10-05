/*
 * SegmentHelper.cpp - Useful routines for segments (group of edges or Surface elems)
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "Omega2D.h"

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>

//
// formula for a parabola between 0..1 peaking in the middle at 1.5 (like viscous flow through an inlet or pipe)
// this may someday go into a "functions" class to aid in defining many inlet velocity distributions
//
float get_viscous_parabola(const float x) {
  return 6.*x*(1.-x);
}

//
// given a set of surface elements (edges), put them in order end-to-end
//   and find the flow through each element based on its parametric position
//
std::vector<float> find_parabolic_vels(const std::vector<float>& _x,
                                       const std::vector<Int>& _idx,
                                       const size_t _n) {

  const size_t nper = _idx.size() / _n;
  assert(nper == 2 && "idx array in find_parabolic_vels is not 2, calculation loses accuracy");

  // find the total arc length, assuming straight edges
  float totlen = 0.0;
  for (size_t i=0; i<_n; ++i) {
    const Int n0 = _idx[2*i];
    const Int n1 = _idx[2*i+1];
    const float dx = _x[2*n0] - _x[2*n1];
    const float dy = _x[2*n0+1] - _x[2*n1+1];
    totlen += std::sqrt(dx*dx + dy*dy);
  }
  std::cout << "    this segment has length " << totlen << std::endl;

  // find the start node by finding which edge has a 0-node that is no other edges' 1-node
  // brute force search is fine here, as _n should always be less than 1000 or so
  // but there are smarter ways to do it
  size_t first_elem = 9999999;
  for (size_t i=0; i<_n; ++i) {
    const Int eleminode0 = _idx[2*i];
    bool itsthere = false;
    for (size_t j=0; j<_n; ++j) {
      const Int elemjnode1 = _idx[2*j+1];
      if (eleminode0 == elemjnode1) {
        itsthere = true;
        break;
      }
    }
    if (not itsthere) {
      first_elem = i;
      break;
    }
  }
  std::cout << "    first elem " << first_elem << " starts at " << _x[2*_idx[2*first_elem]] << " " << _x[2*_idx[2*first_elem]+1] << std::endl;

  // now, march through all elements, putting them in order
  std::vector<float> vels(_n, 0.);
  const float gl2 = 1./std::sqrt(3.);
  float arclen = 0.0;
  size_t ie = first_elem;
  for (size_t i=0; i<_n; ++i) {

    // how long is this element
    const Int n0 = _idx[2*ie];
    const Int n1 = _idx[2*ie+1];
    const float dx = _x[2*n0] - _x[2*n1];
    const float dy = _x[2*n0+1] - _x[2*n1+1];
    const float thislen = std::sqrt(dx*dx + dy*dy);

    // find parametric start and end points
    const float s0 = arclen / totlen;
    const float s1 = (arclen + thislen) / totlen;

    // use Gauss-Legendre quadrature to find mean flow vel through this element
    const float sc = 0.5*(s0+s1);
    const float sd = 0.5*(s1-s0);
    vels[ie] = 0.5 * (get_viscous_parabola(sc+gl2*sd) + get_viscous_parabola(sc-gl2*sd));
    //std::cout << "      vel on " << ie << " is " << vels[ie] << std::endl;

    // set up to find next one
    for (size_t j=0; j<_n; ++j) {
      if (_idx[2*j] == n1) {
        ie = j;
        break;
      }
    }

    // update and continue
    arclen += thislen;
  }

  std::cout << "    vels set on " << _n << " edges on this segment" << std::endl;

  return vels;
}
