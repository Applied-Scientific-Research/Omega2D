/*
 * Kernels.h - Non-class inner kernels for influence calculations
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#include "CoreFunc.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>


//
// velocity influence functions
//
// here is the naming system:
//   kernel_NS_MT
//     N is the number of dimensions of the source element (0=point, 2=surface)
//     M is the number of dimensions of the target element
//     S is the type of the source element ('v'=vortex, 's'=source, 'vs'=vortex and source)
//     T is the type of the target element
//         first character is 'p' for a singular point, 'b' for a vortex blob
//         second character is 'g' if gradients must be returned
//

// thick-cored particle on thick-cored point, no gradients
template <class S, class A> size_t flops_0v_0v () { return 10 + flops_tv_nograds<S>(); }
template <class S, class A>
static inline void kernel_0v_0v (const S sx, const S sy, const S sr, const S ss,
                                 const S tx, const S ty, const S tr,
                                 A* const __restrict__ tu, A* const __restrict__ tv) {
  // 18 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S r2 = ss * core_func<S>(dx*dx + dy*dy, sr, tr);
  *tu -= r2 * dy;
  *tv += r2 * dx;
}

// thick-cored particle on singular point, no gradients
template <class S, class A> size_t flops_0v_0p () { return 10 + flops_tp_nograds<S>(); }
template <class S, class A>
static inline void kernel_0v_0p (const S sx, const S sy, const S sr, const S ss,
                                 const S tx, const S ty,
                                 A* const __restrict__ tu, A* const __restrict__ tv) {
  // 16 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S r2 = ss * core_func<S>(dx*dx + dy*dy, sr);
  *tu -= r2 * dy;
  *tv += r2 * dx;
}


//
// analytic influence of 2d linear constant-strength vortex panel on target point
//   ignoring the 1/2pi factor, which will be multiplied later
//   35 flops average
//
template <class S, class A> size_t flops_1_0v () { return 35; }
template <class S, class A>
static inline void kernel_1_0v (const S sx0, const S sy0,
                                const S sx1, const S sy1,
                                const S str,
                                const S tx, const S ty,
                                A* const __restrict__ tu, A* const __restrict__ tv) {

  // segment vector
  const S dx0   = tx - sx0;
  const S dy0   = ty - sy0;
  const S dx1   = tx - sx1;
  const S dy1   = ty - sy1;

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
#ifdef USE_VC
  const S vstar = S(0.5f) * Vc::log(rij2/rij12);
  S ustar = Vc::atan2(dx1, dy1) - Vc::atan2(dx0, dy0);
  //Vc::where(Vc::isnan(ustar)) | ustar = S(0.0f);
#else
  const S vstar = 0.5f * std::log(rij2/rij12);
  S ustar = std::atan2(dx1, dy1) - std::atan2(dx0, dy0);
#endif
  //std::cout << "ustar started off as " << ustar << std::endl;
#ifdef USE_VC
  Vc::where(ustar < S(-M_PI)) | ustar += S(2.0f*M_PI);
  Vc::where(ustar > S(M_PI)) | ustar -= S(2.0f*M_PI);
#else
  if (ustar < -M_PI) ustar += 2.0f*M_PI;
  if (ustar > M_PI) ustar -= 2.0f*M_PI;
#endif
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  const S px    = sx1-sx0;
  const S py    = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
#ifdef USE_VC
  const S mult  = str * Vc::rsqrt(px*px + py*py);
#else
  const S mult  = str / std::sqrt(px*px + py*py);
#endif
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  *tu = mult*velx;
  *tv = mult*vely;
}


//
// analytic influence of 2d linear constant-strength vortex AND source panel
//   on target point ignoring the 1/2pi factor, which will be multiplied later
//   47 flops average
//
template <class S, class A> size_t flops_1_0vs () { return 47; }
template <class S, class A>
static inline void kernel_1_0vs (const S sx0, const S sy0,
                                 const S sx1, const S sy1,
                                 const S vs, const S ss,
                                 const S tx, const S ty,
                                 A* const __restrict__ tu, A* const __restrict__ tv) {

  // segment vector
  const S dx0   = tx - sx0;
  const S dy0   = ty - sy0;
  const S dx1   = tx - sx1;
  const S dy1   = ty - sy1;

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
#ifdef USE_VC
  const S vstar = S(0.5f) * Vc::log(rij2/rij12);
  S ustar = Vc::atan2(dx1, dy1) - Vc::atan2(dx0, dy0);
#else
  const S vstar = 0.5f * std::log(rij2/rij12);
  S ustar = std::atan2(dx1, dy1) - std::atan2(dx0, dy0);
#endif
  //std::cout << "ustar started off as " << ustar << std::endl;
#ifdef USE_VC
  Vc::where(ustar < S(-M_PI)) | ustar += S(2.0f*M_PI);
  Vc::where(ustar > S(M_PI)) | ustar -= S(2.0f*M_PI);
#else
  if (ustar < -M_PI) ustar += 2.0f*M_PI;
  if (ustar > M_PI) ustar -= 2.0f*M_PI;
#endif
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  S px = sx1-sx0;
  S py = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;
#ifdef USE_VC
  const S mult  = Vc::rsqrt(px*px + py*py);
#else
  const S mult  = 1.0f / std::sqrt(px*px + py*py);
#endif
  px *= mult;
  py *= mult;

  // finally, rotate back into global coordinates
  // and multiply by vortex sheet strength
  *tu = vs * (ustar*px - vstar*py);
  *tv = vs * (ustar*py + vstar*px);
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and now the source strength
  *tu += ss * (ustar*py + vstar*px);
  *tv += ss * (vstar*py - ustar*px);
}


//
// analytic influence of 2d linear constant-strength vortex AND source panel
//   on target point ignoring the 1/2pi factor, and separating out each velocity
//   45 flops average
//
template <class S, class A> size_t flops_1_0vps () { return 45; }
template <class S, class A>
static inline void kernel_1_0vps (const S sx0, const S sy0,
                                  const S sx1, const S sy1,
                                  const S vs, const S ss,
                                  const S tx, const S ty,
                                  A* const __restrict__ vu, A* const __restrict__ vv,
                                  A* const __restrict__ su, A* const __restrict__ sv) {

  // segment vector
  const S dx0   = tx - sx0;
  const S dy0   = ty - sy0;
  const S dx1   = tx - sx1;
  const S dy1   = ty - sy1;

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
#ifdef USE_VC
  const S vstar = S(0.5f) * Vc::log(rij2/rij12);
  S ustar = Vc::atan2(dx1, dy1) - Vc::atan2(dx0, dy0);
#else
  const S vstar = 0.5f * std::log(rij2/rij12);
  S ustar = std::atan2(dx1, dy1) - std::atan2(dx0, dy0);
#endif
  //std::cout << "ustar started off as " << ustar << std::endl;
#ifdef USE_VC
  Vc::where(ustar < S(-M_PI)) | ustar += S(2.0f*M_PI);
  Vc::where(ustar > S(M_PI)) | ustar -= S(2.0f*M_PI);
#else
  if (ustar < -M_PI) ustar += 2.0f*M_PI;
  if (ustar > M_PI) ustar -= 2.0f*M_PI;
#endif
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  S px = sx1-sx0;
  S py = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;
#ifdef USE_VC
  const S mult  = Vc::rsqrt(px*px + py*py);
#else
  const S mult  = 1.0f / std::sqrt(px*px + py*py);
#endif
  px *= mult;
  py *= mult;

  // finally, rotate back into global coordinates
  // and multiply by vortex sheet strength
  *vu = vs * (ustar*px - vstar*py);
  *vv = vs * (ustar*py + vstar*px);
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and now the source strength
  *su = ss * (ustar*py + vstar*px);
  *sv = ss * (vstar*py - ustar*px);
}

