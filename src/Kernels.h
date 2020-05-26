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

// then some useful inlines, to pull out all of the Vc-specific language
#ifdef USE_VC
template <class S>
static inline S get_vstar (const S rij2, const S rij12) {
  return S(0.5f) * Vc::log(rij2/rij12);
}
template <> inline float get_vstar (const float rij2, const float rij12) {
  return 0.5f * std::log(rij2/rij12);
}
template <> inline double get_vstar (const double rij2, const double rij12) {
  return 0.5 * std::log(rij2/rij12);
}
#else
template <class S>
static inline S get_vstar (const S rij2, const S rij12) {
  return 0.5f * std::log(rij2/rij12);
}
#endif

#ifdef USE_VC
template <class S>
static inline S get_ustar (const S dx0, const S dy0, const S dx1, const S dy1) {
  S ustar = Vc::atan2(dx1, dy1) - Vc::atan2(dx0, dy0);
  Vc::where(ustar < S(-M_PI)) | ustar += S(2.0f*M_PI);
  Vc::where(ustar > S(M_PI)) | ustar -= S(2.0f*M_PI);
  return ustar;
}
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
#else
template <class S>
static inline S get_ustar (const S dx0, const S dy0, const S dx1, const S dy1) {
  S ustar = std::atan2(dx1, dy1) - std::atan2(dx0, dy0);
  if (ustar < -M_PI) ustar += 2.0f*M_PI;
  if (ustar > M_PI) ustar -= 2.0f*M_PI;
  return ustar;
}
#endif

// from https://developer.download.nvidia.com/cg/acos.html
// this is 16 flops
#ifdef USE_VC
template <class S>
inline S my_acos(const S _x) {
  S negate = S(0.0f);
  Vc::where(_x < S(0.0f)) | negate = S(1.0f);
  //const S negate = (_x < 0.0f) ? 1.0f : 0.0f;
  S x = Vc::abs(_x);
  Vc::where(x > S(1.0f)) | x = S(1.0f);
  S ret = S(-0.0187293f);
  ret *= x;
  ret += S(0.0742610f);
  ret *= x;
  ret -= S(0.2121144f);
  ret *= x;
  ret += S(1.5707288f);    // NOT pi/2
  ret *= Vc::sqrt(S(1.0f)-x);
  ret -= S(2.0f) * ret * negate;
  return negate * S(M_PI) + ret;
}
template <>
inline float my_acos(const float _x) {
  return std::acos(_x);
}
template <>
inline double my_acos(const double _x) {
  return std::acos(_x);
}
#else
template <class S>
inline S my_acos(const S _x) {
  return std::acos(_x);
}
#endif

#ifdef USE_VC
// this is flops for the Vc version
template <class S> size_t flops_usf () { return 7+17; }
template <class S>
static inline S get_ustar_fast (const S a2, const S b2, const S c2, const S norm) {
  const S numer = b2 + c2 - a2;
  const S denom = S(0.5) * Vc::rsqrt(b2*c2);
  //S ustar = Vc::acos(numer * denom);	// there is no Vc::acos
  S ustar = -my_acos<S>(numer * denom);
  Vc::where(norm < S(0.0)) | ustar = -ustar;
  return ustar;
}
template <>
inline float get_ustar_fast (const float a2, const float b2, const float c2, const float norm) {
  const float numer = b2 + c2 - a2;
  const float denom = 0.5f / std::sqrt(b2*c2);
  float ustar = std::acos(numer * denom);
  return std::copysign(ustar, -norm);
}
#else
template <class S> size_t flops_usf () { return 8; }
template <class S>
static inline S get_ustar_fast (const S a2, const S b2, const S c2, const S norm) {
  const S numer = b2 + c2 - a2;
  assert(b2*c2 > 0 && "Can't take the square root of a negative number; Can't divide by 0");
  const S denom = 0.5f / std::sqrt(b2*c2);
  /*if (abs(numer * denom) > 1) {
      std::cout << "a2: " << a2 << " b2: " << b2 << "c2: " << c2 << std::endl;
  }*/
  assert(abs(numer * denom) <= 1 && "acos takes values in [-1, 1]");
  float ustar = std::acos(numer * denom);
  return std::copysign(ustar, -norm);
  //float ustar = -std::acos(numer * denom);
  //if (norm < 0.0f) ustar = -ustar;
  //return ustar;
}
#endif

#ifdef USE_VC
template <class S>
static inline S get_rsqrt (const S _in) {
  return Vc::rsqrt(_in);
}
template <> inline float get_rsqrt (const float _in) {
  return 1.0f / std::sqrt(_in);
}
template <> inline double get_rsqrt (const double _in) {
  return 1.0 / std::sqrt(_in);
}
#else
template <class S>
static inline S get_rsqrt (const S _in) {
    assert(_in > 0); // Can't take the square root of a negative number; Can't divide by 0
  return 1.0f / std::sqrt(_in);
}
#endif

// first, the flops count (star is 3, rsqrt counts as 1)
template <class S, class A> size_t flops_1_0v () { return 28 + 3 + flops_usf<S>(); }
// then the first kernel
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
  const S px    = sx1- sx0;
  const S py    = sy1- sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
  const S panl  = px*px   + py*py;

  // compute u* and v*
  const S ustar = get_ustar<S>(dx0, dy0, dx1, dy1);
  //const S ustar = get_ustar_fast<S>(panl, rij2, rij12, px*dy0-py*dx0);
  const S vstar = get_vstar<S>(rij2, rij12);
  //std::cout << "ustar is " << ustar << std::endl;
  //std::cout << "ustarf is " << ustarf << std::endl;
  //std::cout << "norm is " << (px*dy0-py*dx0) << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  // assert(panl > 0); // Can't take the square root of a negative number; Can't divide by 0
  const S mult  = str * get_rsqrt<S>(panl);
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
  S px = sx1 - sx0;
  S py = sy1 - sy0;

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
  const S panl  = px*px   + py*py;

  //std::cout << "px is " << px << " and py is " << py << std::endl;
  // assert(panl > 0); // Can't take the square root of a negative number; Can't divide by 0
  const S mult  = get_rsqrt<S>(panl);
  px *= mult;
  py *= mult;

  // compute u* and v*
  const S ustar = get_ustar<S>(dx0, dy0, dx1, dy1);
  //const S ustar = get_ustar_fast<S>(panl, rij2, rij12, px*dy0-py*dx0);
  const S vstar = get_vstar<S>(rij2, rij12);
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

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
  S px = sx1 - sx0;
  S py = sy1 - sy0;

  // side lengths of the triangle s0, s1, t
  const S rij2  = dx0*dx0 + dy0*dy0;
  const S rij12 = dx1*dx1 + dy1*dy1;
  const S panl  = px*px   + py*py;
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  //std::cout << "px is " << px << " and py is " << py << std::endl;
  // assert(panl > 0); // Can't take the square root of a negative number; Can't divide by 0
  const S mult  = get_rsqrt<S>(panl);
  px *= mult;
  py *= mult;

  // compute u* and v*
  const S ustar = get_ustar<S>(dx0, dy0, dx1, dy1);
  //const S ustar = get_ustar_fast<S>(panl, rij2, rij12, px*dy0-py*dx0);
  const S vstar = get_vstar<S>(rij2, rij12);

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

