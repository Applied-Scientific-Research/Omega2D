/*
 * RHS.h - Non-class velocity-to-right-hand-side calculations
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

  //auto start = std::chrono::system_clock::now();
  //float flops = 0.0;

  // size the return vector
  //size_t ntarg  = targ.get_n();
  std::vector<S> rhs;
  rhs.resize(0);

  //auto end = std::chrono::system_clock::now();
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
  const std::array<Vector<S>,Dimensions>& tt = targ.get_tang();
  const std::array<Vector<S>,Dimensions>& tn = targ.get_norm();
  //const Vector<S>&                        tb = targ.get_tang_bcs();

  //std::cout << "tu[0].size() is " << tu[0].size() << std::endl;
  //std::cout << "tx[0].size() is " << tx[0].size() << std::endl;
  //std::cout << "ti.size() is " << ti.size() << std::endl;

  assert(tu[0].size() == tt[0].size() && "Input array sizes do not match");
  assert(tu[0].size() == tn[0].size() && "Input array sizes do not match");
  //assert(tx[0].size() == tb.size() && "Input array sizes do not match");

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
      //rhs[i] -= tb[i];

      //std::cout << "  elem " << i << " vel is " << tu[0][i] << " " << tu[1][i] << std::endl;
      //std::cout << "       " << " tan vec is " << tt[0][i] << " " << tt[1][i] << std::endl;
      //std::cout << "       " << " rhs is " << rhs[i] << std::endl;
    }

  } else {
    // we have unknown vortex and source strengths - use tangential and normal boundary conditions
    for (size_t i=0; i<ntarg; i++) {
      // dot product of normalized tangent with local velocity
      rhs[2*i]   = -(tu[0][i]*tt[0][i] + tu[1][i]*tt[1][i]);
      rhs[2*i+1] = -(tu[0][i]*tn[0][i] + tu[1][i]*tn[1][i]);
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
struct RHSVisitor {
  // source collection, target collection
  std::vector<float> operator()(Points<float> const& targ)   { return vels_to_rhs_points<float>(targ); } 
  std::vector<float> operator()(Surfaces<float> const& targ) { return vels_to_rhs_panels<float>(targ); } 
  std::vector<float> operator()(Volumes<float> const& targ)  { return vels_to_rhs_elems<float>(targ); } 
};

