/*
 * Kernels.h - Non-class influence calculation kernels
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>

//
// Particle stuff
//

//
// non-class function for nbody inner-most loop
//
// sources have 4 values: x, y, str, rad; targets have same
template <class S, class A>
static inline void nbody_kernel_44 (const S* __restrict__ sx, const S* __restrict__ tx, A* __restrict__ tax) {
  // 15 flops
  const S dx = sx[0] - tx[0];
  const S dy = sx[1] - tx[1];
  S r2 = dx*dx + dy*dy + sx[3]*sx[3] + tx[3]*tx[3];
  r2 = sx[2]/r2;
  tax[0] += r2 * dy;
  tax[1] -= r2 * dx;
}

// sources have 4 values: x, y, str, rad; targets have only x, y (radius assumed to be zero)
template <class S, class A>
static inline void nbody_kernel_42 (const S* __restrict__ sx, const S* __restrict__ tx, A* __restrict__ tax) {
  // 13 flops
  const S dx = sx[0] - tx[0];
  const S dy = sx[1] - tx[1];
  S r2 = dx*dx + dy*dy + sx[3]*sx[3];
  r2 = sx[2]/r2;
  tax[0] += r2 * dy;
  tax[1] -= r2 * dx;
}

#ifdef USE_VC
// sources have 4 values: x, y, str, rad; targets have same
template <class S, class A>
static inline void nbody_kernel_Vc_44 (const S& __restrict__ sx, const S& __restrict__ sy,
                                       const S& __restrict__ ss, const S& __restrict__ sr,
                                       const S& __restrict__ tx, const S& __restrict__ ty, const S& __restrict__ tr,
                                       A* __restrict__ tu, A* __restrict__ tv) {
  // 15 flops
  const S dx = sx - tx;
  const S dy = sy - ty;
  S r2 = dx*dx + dy*dy + sr*sr * tr*tr;
  r2 = ss/r2;
  // must specialize this routine if S != A, because then we have 8-dbl and 8-float things here
  (*tu) += r2 * dy;
  (*tv) -= r2 * dx;
}

// sources have 4 values: x, y, str, rad; targets have only x, y (radius assumed to be zero)
template <class S, class A>
static inline void nbody_kernel_Vc_42 (const S& __restrict__ sx, const S& __restrict__ sy,
                                       const S& __restrict__ ss, const S& __restrict__ sr,
                                       const S& __restrict__ tx, const S& __restrict__ ty,
                                       A* __restrict__ tu, A* __restrict__ tv) {
  // 13 flops
  const S dx = sx - tx;
  const S dy = sy - ty;
  S r2 = dx*dx + dy*dy + sr*sr;
  r2 = ss/r2;
  (*tu) += r2 * dy;
  (*tv) -= r2 * dx;
}
#endif


//
// Panel stuff
//

//
// analytic influence of 2d linear constant-strength vortex panel on target point
//   ignoring the 1/2pi factor, which will be multiplied later
// 38 flops (assuming one M_PI check triggers)
//
template <class S, class A>
std::array<A,2> vortex_panel_affects_point(const S sx0, const S sy0,
                                           const S sx1, const S sy1,
                                           const S str,
                                           const S tx, const S ty) {

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(tx-sx0,2) + std::pow(ty-sy0,2);
  //const S rij   = std::sqrt(rij2);
  const S rij12 = std::pow(tx-sx1,2) + std::pow(ty-sy1,2);
  //const S rij1  = std::sqrt(rij12);
  //std::cout << "rij is " << rij << " and rijp1 is " << rij1 << std::endl;
  //const S vstar = std::log(rij/rij1);
  const S vstar = 0.5 * std::log(rij2/rij12);

  S ustar = std::atan2(tx-sx1, ty-sy1) - std::atan2(tx-sx0, ty-sy0);
  //std::cout << "ustar started off as " << ustar << std::endl;
  if (ustar < -M_PI) ustar += 2.*M_PI;
  if (ustar > M_PI) ustar -= 2.*M_PI;
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  const S px    = sx1-sx0;
  const S py    = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  const S mult  = str / std::sqrt(std::pow(px,2) + std::pow(py,2));
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  return std::array<A,2>{{mult*velx, mult*vely}};
}

#ifdef USE_VC
template <class S, class A>
static inline void vortex_panel_affects_point_Vc(const S& sx0, const S& sy0,
                                                 const S& sx1, const S& sy1,
                                                 const S& str,
                                                 const S& tx, const S& ty,
                                                 A* tu, A* tv) {

  // side lengths of the triangle s0, s1, t
  const S dx0   = tx - sx0;
  const S dy0   = ty - sy0;
  const S dx1   = tx - sx1;
  const S dy1   = ty - sy1;
  //std::cout << "dx0 is " << dx0 << std::endl;
  //std::cout << "  dy0 is " << dy0 << std::endl;
  const S rij2  = dx0*dx0 + dy0*dy0;
  //const S rij   = Vc::sqrt(rij2);
  const S rij12 = dx1*dx1 + dy1*dy1;
  //const S rij1  = Vc::sqrt(rij12);
  //const S vstar = Vc::log(rij/rij1);
  const S vstar = S(0.5) * Vc::log(rij2/rij12);

  S ustar = std::atan2(dx1, dy1) - std::atan2(dx0, dy0);
  //std::cout << "ustar started off as " << ustar << std::endl;
  //const Vc::Mask<float> toolow = (ustar < -M_PI);
  //std::cout << "  ustar was " << ustar << std::endl;
  Vc::where(ustar < S(-M_PI)) | ustar += S(2.*M_PI);
  Vc::where(ustar > S(M_PI)) | ustar -= S(2.*M_PI);
  //std::cout << "  ustar is " << ustar << std::endl;
  //std::cout << "  vstar is " << vstar << std::endl;

  const S px    = sx1-sx0;
  const S py    = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  const S mult  = str * Vc::rsqrt(px*px + py*py);
  //std::cout << "  finalx is " << (mult*velx) << std::endl;
  //std::cout << "  finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  (*tu) = mult*velx;
  (*tv) = mult*vely;
  // 37 flops
}
#endif
