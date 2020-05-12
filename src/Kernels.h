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

/*
// gradient kernels here
template <class S, class A> size_t flops_0v_0vg () { return 25 + flops_tv_grads<S>(); }
template <class S, class A>
static inline void kernel_0v_0vg (const S sx, const S sy, const S sr, const S ss,
                                  const S tx, const S ty, const S tr,
                                  A* const __restrict__ tu, A* const __restrict__ tv,
                                  A* const __restrict__ tux, A* const __restrict__ tvx,
                                  A* const __restrict__ tuy, A* const __restrict__ tvy) {
  // 25 flops without core_func
  const S dx = tx - sx;
  const S dy = ty - sy;
  S r2, bbb;
  core_func<S>(dx*dx + dy*dy, sr, tr);
  r2 *= ss;
  bbb *= ss;
  *tu -= r2 * dy;
  *tv += r2 * dx;
  // and the grads
  *tux += -bbb*dx*dy;
  *tuy += -bbb*dy*dy - r2;
  *tvx +=  bbb*dx*dx + r2;
  *tvy +=  bbb*dx*dy;
}
*/


//
// analytic influence of 2d linear constant-strength vortex panel on target point
//   ignoring the 1/2pi factor, which will be multiplied later
//   35 flops average
//
// first, the flops count
template <class S, class A> size_t flops_1_0v () { return 35; }

// then some useful inlines, to pull out all of the Vc-specific language
#ifdef USE_VC
template <class S>
static inline S get_vstar (const S rij2, const S rij12) {
  return S(0.5f) * Vc::log(rij2/rij12);
}
#endif
template <> inline float get_vstar (const float rij2, const float rij12) {
  return 0.5f * std::log(rij2/rij12);
}
template <> inline double get_vstar (const double rij2, const double rij12) {
  return 0.5 * std::log(rij2/rij12);
}

#ifdef USE_VC
template <class S>
static inline S get_ustar (const S dx0, const S dy0, const S dx1, const S dy1) {
  S ustar = Vc::atan2(dx1, dy1) - Vc::atan2(dx0, dy0);
  Vc::where(ustar < S(-M_PI)) | ustar += S(2.0f*M_PI);
  Vc::where(ustar > S(M_PI)) | ustar -= S(2.0f*M_PI);
  return ustar;
}
#endif
template <> inline float get_ustar (const float dx0, const float dy0, const float dx1, const float dy1) {
  float ustar = std::atan2(dx1, dy1) - std::atan2(dx0, dy0);
  if (ustar < -M_PI) ustar += 2.0f*M_PI;
  if (ustar > M_PI) ustar -= 2.0f*M_PI;
  return ustar;
}
template <> inline double get_ustar (const double dx0, const double dy0, const double dx1, const double dy1) {
  double ustar = std::atan2(dx1, dy1) - std::atan2(dx0, dy0);
  if (ustar < -M_PI) ustar += 2.0*M_PI;
  if (ustar > M_PI) ustar -= 2.0*M_PI;
  return ustar;
}

#ifdef USE_VC
template <class S>
static inline S get_rsqrt (const S _in) {
  return Vc::rsqrt(_in);
}
#endif
template <> inline float get_rsqrt (const float _in) {
  return 1.0f / std::sqrt(_in);
}
template <> inline double get_rsqrt (const double _in) {
  return 1.0 / std::sqrt(_in);
}

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
  const S ustar = get_ustar<S>(dx0, dy0, dx1, dy1);

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
  const S vstar = get_vstar<S>(rij2, rij12);
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  const S px    = sx1-sx0;
  const S py    = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  const S mult  = str * get_rsqrt<S>(px*px + py*py);
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
  const S ustar = get_ustar<S>(dx0, dy0, dx1, dy1);

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
  const S vstar = get_vstar<S>(rij2, rij12);
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  S px = sx1-sx0;
  S py = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;
  const S mult  = get_rsqrt<S>(px*px + py*py);
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
  const S ustar = get_ustar<S>(dx0, dy0, dx1, dy1);

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
  const S vstar = get_vstar<S>(rij2, rij12);
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  S px = sx1-sx0;
  S py = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;
  const S mult  = get_rsqrt<S>(px*px + py*py);
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

