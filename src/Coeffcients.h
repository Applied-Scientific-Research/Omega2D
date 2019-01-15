/*
 * Coefficients.h - Non-class influence coefficients calculations
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "NewKernels.h"
#include "Points.h"
#include "Surfaces.h"

#include <algorithm>	// for std::transform
#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>	// for M_PI


template <class S>
Vector<S> points_on_points_coeff (Points<S> const& src, Points<S>& targ) {
  auto start = std::chrono::system_clock::now();

  // when we need this, copy it from Influence.h
  Vector<S> coeffs;
  float flops = 0.0;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_on_points_coeff: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return coeffs;
}

template <class S>
Vector<S> panels_on_points_coeff (Surfaces<S> const& src, Points<S>& targ) {
  std::cout << "    1_0 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;

  Vector<S> coeffs;
  return coeffs;
/*
  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  //const Vector<S>&                      sr = src.get_rad();
  const std::vector<uint16_t>&            si = src.get_idx();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  //const Vector<S>&                      tr = targ.get_rad();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_n(); ++i) {
    //std::array<A,3> accum = {0.0};
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    for (size_t j=0; j<src.get_n(); ++j) {
      const size_t jp0 = si[2*j];
      const size_t jp1 = si[2*j+1];
      //kernel_1_0v<S,A>(&sx[2*si[2*j]], &sx[2*si[2*j+1]], ss[j],
      //                &tx[2*i], accum.data());
      kernel_1_0s<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                       sx[0][jp1], sx[1][jp1], sx[2][jp1],
                       ss[0][j], ss[1][j], ss[2][j],
                       tx[0][i], tx[1][i], tx[2][i],
                       //accum.data());
                       &accumu, &accumv, &accumw);
    }
    //tu[0][i] += accum[0];
    //tu[1][i] += accum[1];
    //tu[2][i] += accum[2];
    tu[0][i] += accumu;
    tu[1][i] += accumv;
    tu[2][i] += accumw;
  }
*/
}


template <class S>
Vector<S> points_on_panels_coeff (Points<S> const& src, Surfaces<S>& targ) {
  std::cout << "    0_1 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;

  Vector<S> coeffs;
  return coeffs;
/*
  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<uint16_t>&            ti = targ.get_idx();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_n(); ++i) {
    //std::array<A,3> accum = {0.0};
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    const size_t ip0 = ti[2*i];
    const size_t ip1 = ti[2*i+1];
    for (size_t j=0; j<src.get_n(); ++j) {
      // note that this is the same kernel as panels_on_points_coeff!
      //kernel_1_0v<S,A>(&tx[2*ti[2*i]], &tx[2*ti[2*i+1]], ss[j],
      //                 &sx[2*j], accum.data());
      kernel_1_0s<S,A>(tx[0][ip0], tx[1][ip0], tx[2][ip0],
                       tx[0][ip1], tx[1][ip1], tx[2][ip1],
                       ss[0][j], ss[1][j], ss[2][j],
                       sx[0][j], sx[1][j], sx[2][j],
                       //accum.data());
                       &accumu, &accumv, &accumw);
    }
    // we use it backwards, so the resulting velocities are negative
    //tu[0][i] -= accum[0];
    //tu[1][i] -= accum[1];
    //tu[2][i] -= accum[2];
    tu[0][i] -= accumu;
    tu[1][i] -= accumv;
    tu[2][i] -= accumw;
  }
*/
}

template <class S>
Vector<S> panels_on_panels_coeff (Surfaces<S> const& src, Surfaces<S>& targ) {
  std::cout << "    1_1 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // how large of a problem do we have?
  // this assumes one unknown per panel - not generally true!!!
  size_t nsrc  = src.get_npanels();
  size_t ntarg = targ.get_npanels();

  // pull references to the element arrays
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::vector<Int>&                 si = src.get_idx();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<Int>&                 ti = targ.get_idx();

  // allocate space for the output array
  Vector<S> coeffs;
  coeffs.resize(nsrc*ntarg);

  // run a panels-on-points algorithm - THIS CAN BE MORE EFFICIENT
  #pragma omp parallel for
  for (size_t j=0; j<nsrc; j++) {
    size_t iptr = ntarg * j;
    const Int sfirst  = si[2*j];
    const Int ssecond = si[2*j+1];
    const S sx0 = sx[0][sfirst];
    const S sy0 = sx[1][sfirst];
    const S sx1 = sx[0][ssecond];
    const S sy1 = sx[1][ssecond];

    for (size_t i=0; i<ntarg; i++) {

      const Int tfirst  = ti[2*i];
      const Int tsecond = ti[2*i+1];
      const S tx0 = tx[0][tfirst];
      const S ty0 = tx[1][tfirst];
      const S tx1 = tx[0][tsecond];
      const S ty1 = tx[1][tsecond];

      // collocation point for panel i
      const S xi = 0.5 * (tx1 + tx0);
      const S yi = 0.5 * (ty1 + ty0);

      // influence of vortex panel j with unit circulation on center of panel i
      auto vel = vortex_panel_affects_point<S,S>(sx0, sy0, sx1, sy1,
                                                 1.0, xi, yi);

      // target panel vector
      const S panelx = tx1 - tx0;
      const S panely = ty1 - ty0;
      const S panell = std::sqrt(panelx*panelx + panely*panely);

      // new way
      // dot product with tangent vector, applying normalization here
      coeffs[iptr++] = (vel[0]*panelx + vel[1]*panely) / panell;

      // old way
      // (un-normalized) surface normal pointing into fluid
      //S normx = -panely;
      //S normy = panelx;
      // dot product with normal vector, applying normalization here
      //A(i,j) = (vel[0]*normx + vel[1]*normy) / panell;
    }

    // special case: self-influence
    coeffs[j*ntarg+j] = M_PI;
  }

  // scale all influences by the constant
  const S fac = 1.0 / (2.0 * M_PI);
  std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
                 [fac](S elem) { return elem * fac; });

  return coeffs;
}


// helper struct for dispatching through a variant
struct CoefficientVisitor {
  // source collection, target collection
  Vector<float> operator()(Points<float> const& src,   Points<float>& targ)   { return points_on_points_coeff<float>(src, targ); } 
  Vector<float> operator()(Surfaces<float> const& src, Points<float>& targ)   { return panels_on_points_coeff<float>(src, targ); } 
  Vector<float> operator()(Points<float> const& src,   Surfaces<float>& targ) { return points_on_panels_coeff<float>(src, targ); } 
  Vector<float> operator()(Surfaces<float> const& src, Surfaces<float>& targ) { return panels_on_panels_coeff<float>(src, targ); } 
};

