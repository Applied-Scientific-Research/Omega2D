/*
 * Coefficients.h - Non-class influence coefficients calculations
 *
 * (c)2017-21 Applied Scientific Research, Inc.
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
      //kernelu_1v_0p<S,A>(&sx[2*si[2*j]], &sx[2*si[2*j+1]], ss[j],
      //                   &tx[2*i], accum.data());
      kernelu_1s_0p<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
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
Vector<S> bricks_on_points_coeff (Volumes<S> const& src, Points<S>& targ) {
  std::cout << "    2_0 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  Vector<S> coeffs;
  return coeffs;
}

// ===========================================================================================================

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
      //kernelu_1v_0p<S,A>(&tx[2*ti[2*i]], &tx[2*ti[2*i+1]], ss[j],
      //                   &sx[2*j], accum.data());
      kernelu_1s_0p<S,A>(tx[0][ip0], tx[1][ip0], tx[2][ip0],
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

  const bool use_two_way = true;

  // use floats to prevent overruns
  float flops = 0.0;
  //auto start = std::chrono::system_clock::now();

  // how large of a problem do we have?
  const size_t nsrc  = src.get_npanels();
  const size_t ntarg = targ.get_npanels();

  // and how many rows and cols does that mean?
  const size_t oldncols =  nsrc * src.num_unknowns_per_panel();
  const size_t oldnrows = ntarg * targ.num_unknowns_per_panel();

  // pull references to the element arrays
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::vector<Int>&                 si = src.get_idx();
  const bool                    src_have_src = (src.num_unknowns_per_panel() == 2);
  const Vector<S>&                        sa = src.get_area();

  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<Int>&                 ti = targ.get_idx();
  const bool                   targ_have_src = (targ.num_unknowns_per_panel() == 2);
  const std::array<Vector<S>,Dimensions>& tt = targ.get_tang();
  const std::array<Vector<S>,Dimensions>& tn = targ.get_norm();
  const Vector<S>&                        ta = targ.get_area();

#ifdef USE_VC
  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
#endif

  // allocate space for the output array
  Vector<S> coeffs;
  coeffs.resize(oldncols*oldnrows);

  // run a panels-on-points algorithm - THIS CAN BE MORE EFFICIENT
  #pragma omp parallel
  {
  Vector<S> col1,col2;
  col1.resize(oldnrows);
  if (src_have_src) col2.resize(oldnrows);

  #pragma omp for
  for (int32_t j=0; j<(int32_t)nsrc; j++) {

    size_t rptr = 0;
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
      StoreVec tx0, ty0, tx1, ty1, ttx, tty, tnx, tny, tva;
      for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
        const size_t idx = i*StoreVec::size() + ii;
        const Int tfirst  = ti[2*idx];
        const Int tsecond = ti[2*idx+1];
        tx0[ii] = tx[0][tfirst];
        ty0[ii] = tx[1][tfirst];
        tx1[ii] = tx[0][tsecond];
        ty1[ii] = tx[1][tsecond];
        ttx[ii] = tt[0][idx];
        tty[ii] = tt[1][idx];
        tnx[ii] = tn[0][idx];
        tny[ii] = tn[1][idx];
        tva[ii] = ta[idx];
      }

      // collocation point for panel i
      const StoreVec xi = StoreVec(0.5) * (tx1 + tx0);
      const StoreVec yi = StoreVec(0.5) * (ty1 + ty0);

      // need to be more sophisticated about this
      if (src_have_src and targ_have_src) {
        // vortex and source strengths
        StoreVec vortu, vortv, srcu, srcv;
        // influence of vortex panel j with unit circulation on center of panel i
        kernelu_1vos_0p<StoreVec,StoreVec>(sx0, sy0, sx1, sy1,
                                           StoreVec(1.0), StoreVec(1.0),
                                           xi, yi,
                                           &vortu, &vortv, &srcu, &srcv);

        // dot product of vortex influence with tangent vector
        StoreVec c1e1 = vortu*ttx + vortv*tty;
        // dot product of vortex influence with normal vector
        StoreVec c1e2 = vortu*tnx + vortv*tny;
        // dot product of source influence with tangent vector
        StoreVec c2e1 = srcu*ttx + srcv*tty;
        // dot product of source influence with normal vector
        StoreVec c2e2 = srcu*tnx + srcv*tny;

        // average this with the point-affects-panel influence
        if (use_two_way) {
          // flipping source and target returns negative of desired influence
          kernelu_1vos_0p<StoreVec,StoreVec>(tx0, ty0, tx1, ty1,
                                             StoreVec(1.0), StoreVec(1.0),
                                             StoreVec(0.5)*(sx0+sx1), StoreVec(0.5)*(sy0+sy1),
                                             &vortu, &vortv, &srcu, &srcv);
          // subtract off the sum
          const StoreVec fac = sa[j]/tva;
          c1e1 -= (vortu*ttx + vortv*tty) * fac;
          c1e2 -= (vortu*tnx + vortv*tny) * fac;
          c2e1 -= ( srcu*ttx +  srcv*tty) * fac;
          c2e2 -= ( srcu*tnx +  srcv*tny) * fac;

          // and take the average
          c1e1 *= StoreVec(0.5);
          c1e2 *= StoreVec(0.5);
          c2e1 *= StoreVec(0.5);
          c2e2 *= StoreVec(0.5);
        }

        // spread the results from a vector register back to the primary array
        for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
          col1[rptr] = c1e1[ii];
          col1[rptr+1] = c1e2[ii];
          col2[rptr] = c2e1[ii];
          col2[rptr+1] = c2e2[ii];
          rptr += 2;
        }

      } else {
        // influence of vortex panel j with unit circulation on center of panel i
        StoreVec vortu, vortv;
        kernelu_1v_0p<StoreVec,StoreVec>(sx0, sy0, sx1, sy1, StoreVec(1.0),
                                         xi, yi, &vortu, &vortv);

        // dot product with tangent vector
        StoreVec newcoeffs = vortu*ttx + vortv*tty;

        // average this with the point-affects-panel influence
        if (use_two_way) {
          // flipping source and target returns negative of desired influence
          kernelu_1v_0p<StoreVec,StoreVec>(tx0, ty0, tx1, ty1, StoreVec(1.0),
                                           StoreVec(0.5)*(sx0+sx1), StoreVec(0.5)*(sy0+sy1),
                                           &vortu, &vortv);
          newcoeffs -= (vortu*ttx + vortv*tty) * (sa[j]/tva);
          // take the average
          newcoeffs *= StoreVec(0.5);
        }

        // spread the results from a vector register back to the primary array
        for (size_t ii=0; ii<StoreVec::size() && i*StoreVec::size()+ii<ntarg; ++ii) {
          col1[rptr++] = newcoeffs[ii];
        }
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

      if (src_have_src and targ_have_src) {
        // vortex and source strengths
        S vortu, vortv, srcu, srcv;
        // influence of vortex panel j with unit circulation on center of panel i
        kernelu_1vos_0p<S,S>(sx0, sy0, sx1, sy1, 1.0, 1.0, xi, yi, &vortu, &vortv, &srcu, &srcv);

        // dot product of vortex influence with tangent vector
        col1[rptr] = vortu*tt[0][i] + vortv*tt[1][i];
        // dot product of vortex influence with normal vector
        col1[rptr+1] = vortu*tn[0][i] + vortv*tn[1][i];
        // dot product of source influence with tangent vector
        col2[rptr] = srcu*tt[0][i] + srcv*tt[1][i];
        // dot product of source influence with normal vector
        col2[rptr+1] = srcu*tn[0][i] + srcv*tn[1][i];

        // average this with the point-affects-panel influence
        if (use_two_way) {
          kernelu_1vos_0p<S,S>(tx0, ty0, tx1, ty1, 1.0, 1.0, 0.5*(sx0+sx1), 0.5*(sy0+sy1), &vortu, &vortv, &srcu, &srcv);
          const S fact = sa[j]/ta[i];
          col1[rptr] -= fact*(vortu*tt[0][i] + vortv*tt[1][i]);
          col1[rptr+1] -= fact*(vortu*tn[0][i] + vortv*tn[1][i]);
          col2[rptr] -= fact*(srcu*tt[0][i] + srcv*tt[1][i]);
          col2[rptr+1] -= fact*(srcu*tn[0][i] + srcv*tn[1][i]);
          // take the average
          col1[rptr] *= 0.5;
          col1[rptr+1] *= 0.5;
          col2[rptr] *= 0.5;
          col2[rptr+1] *= 0.5;
        }

        rptr += 2;

      } else {
        // vortex-strengths only
        // influence of vortex panel j with unit circulation on center of panel i
        S vortu, vortv;
        kernelu_1v_0p<S,S>(sx0, sy0, sx1, sy1, 1.0, xi, yi, &vortu, &vortv);

        // dot product with tangent vector
        col1[rptr] = vortu*tt[0][i] + vortv*tt[1][i];

        // average this with the point-affects-panel influence
        if (use_two_way) {
          // flipping source and target returns negative of desired influence
          kernelu_1v_0p<S,S>(tx0, ty0, tx1, ty1, 1.0, 0.5*(sx0+sx1), 0.5*(sy0+sy1), &vortu, &vortv);
          col1[rptr] -= (vortu*tt[0][i] + vortv*tt[1][i]) * (sa[j]/ta[i]);
          // take the average
          col1[rptr] *= 0.5;
        }

        ++rptr;
      }
    }
#endif

    // special case: self-influence
    if (&src == &targ) {
      if (src_have_src and targ_have_src) {
        // vortex and source strengths
        col1[2*j]   = M_PI;
        col1[2*j+1] = 0.0;
        col2[2*j]   = 0.0;
        col2[2*j+1] = M_PI;
      } else {
        // vortex-strengths only
        col1[j] = M_PI;
      }
    }

    // copy these cols to the matrix
    //for (size_t i=0; i<ntarg; i++) {
      //std::cout << j << " " << i << " " << coeffs[j*ntarg+i] << " " << col1[i] << std::endl;
    //}
    if (src_have_src and targ_have_src) {
      std::copy(col1.begin(), col1.end(), coeffs.begin()+4*j*ntarg);
      std::copy(col2.begin(), col2.end(), coeffs.begin()+4*j*ntarg+2*ntarg);
    } else {
      std::copy(col1.begin(), col1.end(), coeffs.begin()+j*ntarg);
    }

  } // end omp for
  } // end omp parallel

  // update the flop count
  flops += (float)nsrc * (float)ntarg * (use_two_way ? 2.0 : 1.0) *
              (4.0 +
               (float)(src_have_src ? flopsu_1vos_0p<S,S>() : flopsu_1v_0p<S,S>()) +
               (targ_have_src ? 12.0 : 3.0));

  // scale all influences by the constant
  const S fac = 1.0 / (2.0 * M_PI);
  std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
                 [fac](S elem) { return elem * fac; });
  flops += 2.0 + (float)coeffs.size();

  // skip out if we don't augment
  if (not targ.is_augmented() and not src.is_augmented()) return coeffs;


  //
  // now we augment this "matrix" with an optional new row and column
  //

  const size_t nrows = oldnrows + (targ.is_augmented() ? 1 : 0);
  const size_t ncols = oldncols + ( src.is_augmented() ? 1 : 0);
  if (targ.is_augmented() or src.is_augmented()) {
    std::cout << "    augmenting the " << ntarg << " x " << nsrc << " block to " << nrows << " x " << ncols << std::endl;
  }

  bool debug = false;
  //if (targ.is_augmented() or src.is_augmented()) debug = true;

  // debug print the bottom-right corner
  if (debug) {
    std::cout << "Bottom-right corner of influence matrix:" << std::endl;
    for (size_t i=oldnrows-6; i<oldnrows; ++i) {
      for (size_t j=oldncols-6; j<oldncols; ++j) {
        std::cout << " \t" << coeffs[oldnrows*j+i];
      }
      std::cout << std::endl;
    }
  }

  // make a new 1-D vector to contain the coefficients
  Vector<S> augcoeff;
  augcoeff.resize(nrows*ncols);

  // first, copy the old vector into the new and generate the new bottom row

  // original vector is column-major (src/column index changes slowest), so let's keep that
  for (size_t j=0; j<oldncols; j++) {

    // iterators in each vector
    auto c_iter = coeffs.begin() + j*oldnrows;
    auto a_iter = augcoeff.begin() + j*nrows;

    // copy the next ntarg (*1 or *2) numbers into the new vector
    std::copy(c_iter, c_iter+oldnrows, a_iter);

    // and add the bottom value to this column
    if (targ.is_augmented()) {
      // always include the panel lengths of the source body
      a_iter += oldnrows;
      // then write the last value in this column - the length of this panel
      // coefficient in matrix is the panel length
      if (&src == &targ) {
        if (targ.num_unknowns_per_panel() == 1) {
          *a_iter = sa[j];
        } else if (j%2 == 0) {
          // but only if an even number if V+S
          *a_iter = sa[j/2];
        } else {
          *a_iter = 0.0;
        }
      } else {
        *a_iter = 0.0;
      }
    }
  }

  // no longer need coeffs
  coeffs.clear();

  // debug print the bottom-right corner
  if (debug) {
    std::cout << "Bottom-right corner of influence matrix:" << std::endl;
    for (size_t i=nrows-6; i<nrows; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << augcoeff[nrows*j+i];
      }
      std::cout << std::endl;
    }
  }

  // then add the last column, if necessary, if target rotates
  if (src.is_augmented()) {
    auto a_iter = augcoeff.begin() + oldncols*nrows;

    // this is the velocity influence from the source body with unit rotational rate on these target panels
    std::array<Vector<S>,Dimensions> const& vel = targ.get_vel();

    // fill in the entries - the influence of the source body's panels, when vort and source terms are set
    //   such that the integration results in the flow imposed by the body's volume of vorticity,
    //   on this target panel.
    for (size_t i=0; i<ntarg; ++i) {

      // find target point - just above the panel
      // yes, I know I am calling these "source"
      // and find the component of that velocity along the "target" panel
      *a_iter = vel[0][i]*tt[0][i] + vel[1][i]*tt[1][i];
      ++a_iter;

      if (targ_have_src) {
        *a_iter = vel[0][i]*tn[0][i] + vel[1][i]*tn[1][i];
        ++a_iter;
      }
    }

    // finally, the bottom corner is the circulation at unit rotation of the body
    if (&src == &targ) {
      *a_iter = 2.0 * src.get_vol();
    } else {
      *a_iter = 0.0;
    }
  } else {
    // if there's no body pointer, then what?
  }

  // debug print the top-left and bottom-right corners
  if (debug) {
    //std::cout << "Top-left corner of influence matrix:" << std::endl;
    //for (size_t i=0; i<6; ++i) {
    //  for (size_t j=0; j<6; ++j) {
    //    std::cout << " \t" << augcoeff[nrows*j+i];
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << "Top-right corner of influence matrix:" << std::endl;
    //for (size_t i=0; i<6; ++i) {
    //  for (size_t j=ncols-6; j<ncols; ++j) {
    //    std::cout << " \t" << augcoeff[nrows*j+i];
    //  }
    //  std::cout << std::endl;
    //}
    std::cout << "Bottom-right corner of influence matrix:" << std::endl;
    for (size_t i=nrows-6; i<nrows; ++i) {
      for (size_t j=ncols-6; j<ncols; ++j) {
        std::cout << " \t" << augcoeff[nrows*j+i];
      }
      std::cout << std::endl;
    }
  }

  //auto end = std::chrono::system_clock::now();
  //std::chrono::duration<double> elapsed_seconds = end-start;
  //const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  //printf("    matrix block:\t[%.4f] cpu seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return augcoeff;
}

template <class S>
Vector<S> bricks_on_panels_coeff (Volumes<S> const& src, Surfaces<S>& targ) {
  std::cout << "    2_1 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  Vector<S> coeffs;
  return coeffs;
}

// ===========================================================================================================

// bricks (volume elements) do not participate in BEM yet
template <class S>
Vector<S> points_on_bricks_coeff (Points<S> const& src, Volumes<S>& targ) {
  std::cout << "    0_2 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  Vector<S> coeffs;
  return coeffs;
}

template <class S>
Vector<S> panels_on_bricks_coeff (Surfaces<S> const& src, Volumes<S>& targ) {
  std::cout << "    1_2 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  Vector<S> coeffs;
  return coeffs;
}

template <class S>
Vector<S> bricks_on_bricks_coeff (Volumes<S> const& src, Volumes<S>& targ) {
  std::cout << "    2_2 compute coefficients of" << src.to_string() << " on" << targ.to_string() << std::endl;
  Vector<S> coeffs;
  return coeffs;
}


// ===========================================================================================================


// helper struct for dispatching through a variant
template <class S>
struct CoefficientVisitor {
  // source collection, target collection
  Vector<S> operator()(Points<S> const& src,   Points<S>& targ)   { return points_on_points_coeff<S>(src, targ); } 
  Vector<S> operator()(Surfaces<S> const& src, Points<S>& targ)   { return panels_on_points_coeff<S>(src, targ); } 
  Vector<S> operator()(Volumes<S> const& src,  Points<S>& targ)   { return bricks_on_points_coeff<S>(src, targ); } 
  Vector<S> operator()(Points<S> const& src,   Surfaces<S>& targ) { return points_on_panels_coeff<S>(src, targ); } 
  Vector<S> operator()(Surfaces<S> const& src, Surfaces<S>& targ) { return panels_on_panels_coeff<S>(src, targ); } 
  Vector<S> operator()(Volumes<S> const& src,  Surfaces<S>& targ) { return bricks_on_panels_coeff<S>(src, targ); } 
  Vector<S> operator()(Points<S> const& src,   Volumes<S>& targ)  { return points_on_bricks_coeff<S>(src, targ); } 
  Vector<S> operator()(Surfaces<S> const& src, Volumes<S>& targ)  { return panels_on_bricks_coeff<S>(src, targ); } 
  Vector<S> operator()(Volumes<S> const& src,  Volumes<S>& targ)  { return bricks_on_bricks_coeff<S>(src, targ); } 
};

