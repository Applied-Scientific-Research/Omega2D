/*
 * RHS.h - Non-class velocity-to-right-hand-side calculations
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "Points.h"
#include "Surfaces.h"

#include <iostream>
#include <vector>
#include <cassert>


// No need to return anything because points currently cannot be reactive
template <class S>
std::vector<S> vels_to_rhs_points (Points<S> const& targ) {
  std::cout << "    NOT converting vels to RHS vector for " << targ.to_string() << std::endl;

  //auto start = std::chrono::steady_clock::now();
  //float flops = 0.0;

  // size the return vector
  //size_t ntarg  = targ.get_n();
  std::vector<S> rhs;
  rhs.resize(0);

  //auto end = std::chrono::steady_clock::now();
  //std::chrono::duration<double> elapsed_seconds = end-start;
  //const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  //printf("    points_on_points_coeff: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return rhs;
}

template <class S>
std::vector<S> vels_to_rhs_panels (Surfaces<S> const& targ) {
  std::cout << "    convert vels to RHS vector for" << targ.to_string() << std::endl;

  // pull references to the element arrays
  const std::array<Vector<S>,Dimensions>& tu = targ.get_vel();
  const std::array<Vector<S>,Dimensions>& tt = targ.get_tang();		// tangential vectors
  const std::array<Vector<S>,Dimensions>& tn = targ.get_norm();		// normal vectors
  //const Vector<S>&                      ttbc = targ.get_tang_bcs();
  const Vector<S>&                      tnbc = targ.get_norm_bcs();

  assert(tu[0].size() == tt[0].size() && "Input array sizes do not match");
  assert(tu[0].size() == tn[0].size() && "Input array sizes do not match");

  // find array sizes
  const size_t ntarg = targ.get_npanels();
  const size_t nunk  = targ.num_unknowns_per_panel();

  // prepare the rhs vector
  std::vector<S> rhs;
  rhs.resize(ntarg*nunk);

  // convert velocity and boundary condition to RHS values

  if (nunk == 1) {
    // we have vortex strengths only - use the tangential boundary condition
    for (size_t i=0; i<ntarg; i++) {
      // dot product of normalized tangent with local velocity
      rhs[i] = -(tu[0][i]*tt[0][i] + tu[1][i]*tt[1][i]);

      // DO NOT include the influence of the boundary condition here, do it before shedding
      //rhs[i] -= ttbc[i];

      //std::cout << "  elem " << i << " vel is " << tu[0][i] << " " << tu[1][i] << std::endl;
      //std::cout << "       " << " tan vec is " << tt[0][i] << " " << tt[1][i] << std::endl;
      //std::cout << "       " << " rhs is " << rhs[i] << std::endl;
    }

  } else {
    // we have unknown vortex and source strengths - use tangential and normal boundary conditions
    assert(tu[0].size() == tnbc.size() && "Input array sizes do not match");

    for (size_t i=0; i<ntarg; i++) {
      // dot product of normalized tangent with local velocity
      rhs[2*i]   = -(tu[0][i]*tt[0][i] + tu[1][i]*tt[1][i]);
      // like above, do not add the tangential BC here, because it will be accounted for during shedding ?!?

      // dot prod of normal vector with local vel
      rhs[2*i+1] = -(tu[0][i]*tn[0][i] + tu[1][i]*tn[1][i]);
      // but here we need to add the normal velocity boundary condition
      rhs[2*i+1] += tnbc[i];

      //std::cout << "  elem " << i << " vel is " << tu[0][i] << " " << tu[1][i] << std::endl;
      //std::cout << "       " << " tan vec is " << tt[0][i] << " " << tt[1][i] << std::endl;
      //std::cout << "       " << " nrm vec is " << tn[0][i] << " " << tn[1][i] << std::endl;
      //std::cout << "       " << " tnbc is " << tnbc[i] << std::endl;
      //std::cout << "       " << " rhs is " << rhs[2*i] << " " << rhs[2*i+1] << std::endl;
    }
  }

  return rhs;
}

// No need to return anything because brick elems currently cannot be reactive
template <class S>
std::vector<S> vels_to_rhs_elems (Volumes<S> const& targ) {
  std::cout << "    NOT converting vels to RHS vector for " << targ.to_string() << std::endl;

  // size the return vector
  //size_t ntarg  = targ.get_nelems();
  std::vector<S> rhs;
  rhs.resize(0);

  return rhs;
}


// helper struct for dispatching through a variant
template <class S>
struct RHSVisitor {
  // source collection, target collection
  std::vector<S> operator()(Points<S> const& targ)   { return vels_to_rhs_points<S>(targ); } 
  std::vector<S> operator()(Surfaces<S> const& targ) { return vels_to_rhs_panels<S>(targ); } 
  std::vector<S> operator()(Volumes<S> const& targ)  { return vels_to_rhs_elems<S>(targ); } 
};

