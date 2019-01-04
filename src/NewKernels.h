/*
 * Kernels.h - Non-class inner kernels for influence calculations
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

//#define _USE_MATH_DEFINES
//#include <cmath>


// velocity influence functions

template <class S, class A>
static inline void kernel_0v_0v (const S sx, const S sy, const S sr, const S ss,
                                 const S tx, const S ty, const S tr,
                                 A* __restrict__ tu, A* __restrict__ tv) {
  // 15 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  S r2 = dx*dx + dy*dy + sr*sr + tr*tr;
#ifdef USE_VC
  r2 = ss * Vc::reciprocal(r2);
#else
  r2 = ss / r2;
#endif
  *tu -= r2 * dy;
  *tv += r2 * dx;
}


template <class S, class A>
static inline void kernel_0v_0p (const S sx, const S sy, const S sr, const S ss,
                                 const S tx, const S ty,
                                 A* __restrict__ tu, A* __restrict__ tv) {
  // 13 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  S r2 = dx*dx + dy*dy + sr*sr;
#ifdef USE_VC
  r2 = ss * Vc::reciprocal(r2);
#else
  r2 = ss / r2;
#endif
  *tu -= r2 * dy;
  *tv += r2 * dx;
}


/*
//
// analytic influence of 2d linear constant-strength vortex panel on target point
//   ignoring the 1/2pi factor, which will be multiplied later
//   40 flops average
//
template <class S, class A>
static inline void kernel_1_0v (const S* __restrict__ sx0, const S* __restrict__ sx1, const S str,
                                const S* __restrict__ tx, A* __restrict__ tu) {

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(tx[0]-sx0[0],2) + std::pow(tx[1]-sx0[1],2);
  const S rij   = std::sqrt(rij2);
  const S rij12 = std::pow(tx[0]-sx1[0],2) + std::pow(tx[1]-sx1[1],2);
  const S rij1  = std::sqrt(rij12);
  //std::cout << "rij is " << rij << " and rijp1 is " << rij1 << std::endl;
  const S vstar = std::log(rij/rij1);
  S ustar = std::atan2(tx[0]-sx1[0], tx[1]-sx1[1]) - std::atan2(tx[0]-sx0[0], tx[1]-sx0[1]);
  //std::cout << "ustar started off as " << ustar << std::endl;
  if (ustar < -M_PI) ustar += 2.*M_PI;
  if (ustar > M_PI) ustar -= 2.*M_PI;
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  const S px    = sx1[0]-sx0[0];
  const S py    = sx1[1]-sx0[1];
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  const S mult  = str / std::sqrt(std::pow(px,2) + std::pow(py,2));
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  tu[0] += mult*velx;
  tu[1] += mult*vely;
}

template <class S, class A>
static inline void kernel_1_0s (const S sx0, const S sy0, const S sz0,
                                const S sx1, const S sy1, const S sz1,
                                const S ssx, const S ssy, const S ssz,
                                const S tx, const S ty, const S tz,
                                A* __restrict__ tu, A* __restrict__ tv, A* __restrict__ tw) {

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(tx-sx0,2) + std::pow(ty-sy0,2);
  const S rij   = std::sqrt(rij2);
  const S rij12 = std::pow(tx-sx1,2) + std::pow(ty-sy1,2);
  const S rij1  = std::sqrt(rij12);
  //std::cout << "rij is " << rij << " and rijp1 is " << rij1 << std::endl;
  const S vstar = std::log(rij/rij1);
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
  const S mult  = ssx / std::sqrt(std::pow(px,2) + std::pow(py,2));
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  *tu += mult*velx;
  *tv += mult*vely;
}
*/
