/*
 * Coefficients.h - Non-class influence coefficients calculations
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "Kernels.h"
#include "Points.h"
#include "Surfaces.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

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
  const size_t nsrc  = src.get_npanels();
  const size_t ntarg = targ.get_npanels();

  // pull references to the element arrays
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::vector<Int>&                 si = src.get_idx();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<Int>&                 ti = targ.get_idx();

#ifdef USE_VC
  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
#endif

  // allocate space for the output array
  Vector<S> coeffs;
  coeffs.resize(nsrc*ntarg);

  // run a panels-on-points algorithm - THIS CAN BE MORE EFFICIENT
  #pragma omp parallel for
  for (size_t j=0; j<nsrc; j++) {
    size_t iptr = ntarg * j;
    const Int sfirst  = si[2*j];
    const Int ssecond = si[2*j+1];

#ifdef USE_VC
    const StoreVec sx0 = sx[0][sfirst];
    const StoreVec sy0 = sx[1][sfirst];
    const StoreVec sx1 = sx[0][ssecond];
    const StoreVec sy1 = sx[1][ssecond];

    const size_t ntargvec = 1 + (ntarg-1) / StoreVec::size();

    for (size_t i=0; i<ntargvec; i++) {

      // fill a 4- or 8-wide vector with the target coordinates
      StoreVec tx0, ty0, tx1, ty1;
      for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
        const size_t idx = i*StoreVec::size() + ii;
        const Int tfirst  = ti[2*idx];
        const Int tsecond = ti[2*idx+1];
        tx0[ii] = tx[0][tfirst];
        ty0[ii] = tx[1][tfirst];
        tx1[ii] = tx[0][tsecond];
        ty1[ii] = tx[1][tsecond];
      }

      // collocation point for panel i
      const StoreVec xi = StoreVec(0.5) * (tx1 + tx0);
      const StoreVec yi = StoreVec(0.5) * (ty1 + ty0);

      // influence of vortex panel j with unit circulation on center of panel i
      StoreVec resultu, resultv;
      kernel_1_0v<StoreVec,StoreVec>(sx0, sy0, sx1, sy1, StoreVec(1.0),
                                     xi, yi, &resultu, &resultv);

      // target panel vector
      const StoreVec panelx = tx1 - tx0;
      const StoreVec panely = ty1 - ty0;
      const StoreVec panell = Vc::sqrt(panelx*panelx + panely*panely);

      // dot product with tangent vector, applying normalization here
      const StoreVec newcoeffs = (resultu*panelx + resultv*panely) / panell;

      // spread the results from a vector register back to the primary array
      for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
        coeffs[iptr++] = newcoeffs[ii];
      }
    }
#else	// no Vc
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
      S resultu, resultv;
      kernel_1_0v<S,S>(sx0, sy0, sx1, sy1, 1.0, xi, yi, &resultu, &resultv);

      // target panel vector
      const S panelx = tx1 - tx0;
      const S panely = ty1 - ty0;
      const S panell = std::sqrt(panelx*panelx + panely*panely);

      // dot product with tangent vector, applying normalization here
      coeffs[iptr++] = (resultu*panelx + resultv*panely) / panell;
    }
#endif

    // special case: self-influence
    if (&src == &targ) coeffs[j*ntarg+j] = M_PI;
  }

  // scale all influences by the constant
  const S fac = 1.0 / (2.0 * M_PI);
  std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
                 [fac](S elem) { return elem * fac; });


  //
  // now we augment this "matrix" with an optional new row and column
  //

  const size_t nrows = ntarg + (targ.get_body_ptr() ? 1 : 0);
  const size_t ncols = nsrc  + ( src.get_body_ptr() ? 1 : 0);
  std::cout << "    augmenting the " << ntarg << " x " << nsrc << " block to " << nrows << " x " << ncols << std::endl;

  // make a new 1-D vector to contain the coefficients
  Vector<S> augcoeff;
  augcoeff.resize(nrows*ncols);

  // first, copy the old vector into the new and generate the new bottom row

  // original vector is column-major (src/column index changes slowest), so let's keep that
  for (size_t j=0; j<nsrc; j++) {

    // iterators in each vector
    auto c_iter = coeffs.begin() + j*ntarg;
    auto a_iter = augcoeff.begin() + j*nrows;

    // copy the next nsrc numbers into the new vector
    std::copy(c_iter, c_iter+ntarg, a_iter);

    if (src.get_body_ptr()) {
      a_iter += nsrc;
      if (&src == &targ) {
        // then write the last value in this column - the length of this panel
        const Int tfirst  = ti[2*j];
        const Int tsecond = ti[2*j+1];
        // target panel vector
        const S panelx = tx[0][tsecond] - tx[0][tfirst];
        const S panely = tx[1][tsecond] - tx[1][tfirst];
        const S panell = std::sqrt(panelx*panelx + panely*panely);
        // coefficient in matrix is the panel length
        *a_iter = panell;
      } else {
        *a_iter = 0.0;
      }
    }
  }

  // no longer need coeffs
  coeffs.clear();

  // then add the last column, if necessary, if target rotates
  if (targ.get_body_ptr()) {
    auto a_iter = augcoeff.begin() + nsrc*nrows;

    // get the rotational center of this Surface
    std::array<S,2> rc = targ.get_geom_center();

    // fill in the entries - the influence of the source body's panels, when vort and source terms are set
    //   such that the integration results in the flow imposed by the body's volume of vorticity,
    //   on this target panel.
    for (size_t i=0; i<ntarg; ++i) {

      // find target point - just above the panel
      // yes, I know I am calling these "source"
      const Int tfirst  = ti[2*i];
      const Int tsecond = ti[2*i+1];
      const S panelx = tx[0][tsecond] - tx[0][tfirst];
      const S panely = tx[1][tsecond] - tx[1][tfirst];
      const S panell = std::sqrt(panelx*panelx + panely*panely);
      //std::cout << "targ panel " << i << " at " << tx[0][tfirst] << " " << tx[1][tfirst] << std::endl;

#ifdef USE_VC
      // make an 8-wide vector of the target point
      const StoreVec tpx = 0.5 * (tx[0][tsecond] + tx[0][tfirst]) - 0.01*panely;
      const StoreVec tpy = 0.5 * (tx[1][tsecond] + tx[1][tfirst]) + 0.01*panelx;

      // running velocity sum
      StoreVec velx = 0.0;
      StoreVec vely = 0.0;

      const size_t nsrcvec = 1 + (nsrc-1) / StoreVec::size();
      for (size_t j=0; j<nsrcvec; ++j) {

        // fill a 4- or 8-wide vector with the source coordinates
        StoreVec sx0 = StoreVec(-1.0);
        StoreVec sy0 = StoreVec(0.0);
        StoreVec sx1 = StoreVec(1.0);
        StoreVec sy1 = StoreVec(0.0);
        for (size_t jj=0; jj<StoreVec::size() && j*StoreVec::size()+jj<nsrc; ++jj) {
          const size_t idx = j*StoreVec::size() + jj;
          const Int sfirst  = si[2*idx];
          const Int ssecond = si[2*idx+1];
          sx0[jj] = sx[0][sfirst];
          sy0[jj] = sx[1][sfirst];
          sx1[jj] = sx[0][ssecond];
          sy1[jj] = sx[1][ssecond];
        }

        // set the vortex and source strengths for this panel
        // this is the vector from the geometric center to the center of this panel
        const StoreVec spx = StoreVec(0.5) * (sx0 + sx1) - StoreVec(rc[0]);
        const StoreVec spy = StoreVec(0.5) * (sy0 + sy1) - StoreVec(rc[1]);

        // velocity of panel center under unit rotational rate
        const StoreVec upc = -spy;
        const StoreVec vpc =  spx;
        // and this is the tangential vector along this panel (easy to find normal)
        StoreVec spanelx = sx1 - sx0;
        StoreVec spanely = sy1 - sy0;
        const StoreVec spanell = Vc::rsqrt(spanelx*spanelx + spanely*spanely);
        spanelx *= spanell;
        spanely *= spanell;

        // finally can compute the strengths - making sure that the unfilled portion of the vector is 0.0
        StoreVec str_vor = -upc*spanelx - vpc*spanely;
        StoreVec str_src = -upc*spanely + vpc*spanelx;
        for (size_t jj=nsrc-j*StoreVec::size(); jj<StoreVec::size(); ++jj) {
          str_vor[jj] = 0.0;
          str_src[jj] = 0.0;
        }

        // influence of vortex panel j with unit circulation on center of panel i
        StoreVec resultu, resultv;
        kernel_1_0vs<StoreVec,StoreVec>(sx0, sy0, sx1, sy1, str_vor, str_src, tpx, tpy, &resultu, &resultv);
        //std::cout << "  src vec " << i << " uvel " << resultu << std::endl;
        velx += resultu;
        vely += resultv;

        // NOTE: I need to account for the source term if this is to work on non-circular bodies
      }

      // and find the component of that velocity along the "target" panel
      *a_iter = 0.0;
      // spread the results from a vector register back to the primary array
      for (size_t jj=0; jj<StoreVec::size() && jj<nsrc; ++jj) {
        *a_iter += fac * (velx[jj]*panelx + vely[jj]*panely) / panell;
      }
      ++a_iter;

#else	// no Vc
      const S tpx = 0.5 * (tx[0][tsecond] + tx[0][tfirst]) - 0.01*panely;
      const S tpy = 0.5 * (tx[1][tsecond] + tx[1][tfirst]) + 0.01*panelx;

      // running velocity sum
      std::array<S,2> vel = {0.0, 0.0};

      // integrate over all panels on the other body, with source and vortex terms set
      for (size_t j=0; j<nsrc; ++j) {

        const Int sfirst  = si[2*j];
        const Int ssecond = si[2*j+1];
        const S sx0 = sx[0][sfirst];
        const S sy0 = sx[1][sfirst];
        const S sx1 = sx[0][ssecond];
        const S sy1 = sx[1][ssecond];

        // set the vortex and source strengths for this panel
        // this is the vector from the geometric center to the center of this panel
        const S spx = 0.5 * (sx0 + sx1) - rc[0];
        const S spy = 0.5 * (sy0 + sy1) - rc[1];
        // velocity of panel center under unit rotational rate
        const S upc = -spy;
        const S vpc =  spx;
        // and this is the tangential vector along this panel (easy to find normal)
        S spanelx = sx1 - sx0;
        S spanely = sy1 - sy0;
        const S spanell = 1.0/std::sqrt(spanelx*spanelx + spanely*spanely);
        spanelx *= spanell;
        spanely *= spanell;

        //const S str_vor = 0.5;	// for D=1 circle case, approx.
        const S str_vor = -upc*spanelx - vpc*spanely;
        const S str_src = -upc*spanely + vpc*spanelx;

        // idea - have Surface::find_rotation_strengths() return a std::array<std::vector<S>,2> of strengths

        // influence of vortex panel j with given vortex strength on center of panel i
        S resultu, resultv;
        kernel_1_0vs<S,S>(sx0, sy0, sx1, sy1, str_vor, str_src, tpx, tpy, &resultu, &resultv);
        vel[0] += resultu;
        vel[1] += resultv;

        // NOTE: I need to account for the source term if this is to work on non-circular bodies
      }

      // and find the component of that velocity along the target panel
      *a_iter = fac * (vel[0]*panelx + vel[1]*panely) / panell;
       ++a_iter;
#endif	// no Vc
    }

    // finally, the bottom corner is the circulation at unit rotation of the body
    if (&src == &targ) {
      *a_iter = 2.0 * src.get_area();
    } else {
      *a_iter = 0.0;
    }
  } else {
    // if there's no body pointer, then what?
  }

  // debug print the top-left and bottom-right corners
  if (true) {
    //std::cout << "Top-left corner of influence matrix:" << std::endl;
    //for (size_t i=0; i<6; ++i) {
    //  for (size_t j=0; j<6; ++j) {
    //    std::cout << " \t" << augcoeff[nrows*j+i];
    //  }
    //  std::cout << std::endl;
    //}
    std::cout << "Top-right corner of influence matrix:" << std::endl;
    for (size_t i=0; i<6; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << augcoeff[nrows*j+i];
      }
      std::cout << std::endl;
    }
    std::cout << "Bottom-right corner of influence matrix:" << std::endl;
    for (size_t i=nrows-6; i<nrows; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << augcoeff[nrows*j+i];
      }
      std::cout << std::endl;
    }
  }

  return augcoeff;
}


// helper struct for dispatching through a variant
struct CoefficientVisitor {
  // source collection, target collection
  Vector<float> operator()(Points<float> const& src,   Points<float>& targ)   { return points_on_points_coeff<float>(src, targ); } 
  Vector<float> operator()(Surfaces<float> const& src, Points<float>& targ)   { return panels_on_points_coeff<float>(src, targ); } 
  Vector<float> operator()(Points<float> const& src,   Surfaces<float>& targ) { return points_on_panels_coeff<float>(src, targ); } 
  Vector<float> operator()(Surfaces<float> const& src, Surfaces<float>& targ) { return panels_on_panels_coeff<float>(src, targ); } 
};

