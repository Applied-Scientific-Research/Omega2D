/*
 * CoreFunc.h - Non-class core function inlines for influence calculations
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#include "MathHelper.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>

//#define USE_RM_KERNEL
//#define USE_EXPONENTIAL_KERNEL
#define USE_WL_KERNEL
//#define USE_V2_KERNEL
//#define USE_V3_KERNEL


#ifdef USE_RM_KERNEL
//
// Rosenhead-Moore velocity-only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = distsq + sr*sr;
  return my_recip(r2);
}
template <class S> size_t flops_tp_nograds () { return 3; }

// and the one for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = distsq + sr*sr + tr*tr;
  return my_recip(r2);
}
template <class S> size_t flops_tv_nograds () { return 5; }

//
// Rosenhead-Moore with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S td2 = distsq + sr*sr;
  *r2 = my_recip(td2);
  *bbb = S(-2.0) * (*r2) * (*r2);
}
template <class S> size_t flops_tp_grads () { return 9; }

// and the one for non-singular targets
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S td2 = distsq + sr*sr + tr*tr;
  *r2 = my_recip(td2);
  *bbb = S(-2.0) * (*r2) * (*r2);
}
template <class S> size_t flops_tv_grads () { return 9; }
#endif


#ifdef USE_EXPONENTIAL_KERNEL
//
// exponential core - velocity only
//
#ifdef USE_VC
template <class S>
static inline S exp_cond (const S ood2, const S corefac, const S reld2) {
  S returnval = ood2;
  returnval(reld2 < S(16.0)) = ood2 * (S(1.0) - Vc::exp(-reld2));
  returnval(reld2 < S(0.001)) = corefac;
  return returnval;
}
template <>
inline float exp_cond (const float ood2, const float corefac, const float reld2) {
  if (reld2 > 16.0f) {
    return ood2;
  } else if (reld2 < 0.001f) {
    return corefac;
  } else {
    return ood2 * (1.0f - std::exp(-reld2));
  }
}
template <>
inline double exp_cond (const double ood2, const double corefac, const double reld2) {
  if (reld2 > 16.0) {
    return ood2;
  } else if (reld2 < 0.001) {
    return corefac;
  } else {
    return ood2 * (1.0 - std::exp(-reld2));
  }
}
#else
template <class S>
static inline S exp_cond (const S ood2, const S corefac, const S reld2) {
  if (reld2 > 16.0f) {
    return ood2;
    // 1 flop (comparison)
  } else if (reld2 < 0.001f) {
    return corefac;
    // 2 flops
  } else {
    return ood2 * (1.0f - std::exp(-reld2));
    // 3 flops
  }
}
#endif

template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S ood2 = my_recip(distsq);
  const S corefac = my_recip(sr*sr);
  const S reld2 = corefac / ood2;
  // 4 flops to here
  return exp_cond(ood2, corefac, reld2);
}
template <class S> size_t flops_tp_nograds () { return 7; }

// non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S ood2 = my_recip(distsq);
  const S corefac = my_recip(sr*sr + tr*tr);
  const S reld2 = corefac / ood2;
  return exp_cond(ood2, corefac, reld2);
}
template <class S> size_t flops_tv_nograds () { return 9; }

//
// exponential core - with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S ood2 = my_recip(distsq);
  const S corefac = my_recip(sr*sr);
  const S reld2 = corefac / ood2;
  // 4 flops to here
  *r2 = exp_cond(ood2, corefac, reld2);
  *bbb = S(-2.0) * (*r2) * ood2
       + S(2.0) * corefac * ood2 * my_exp(-reld2);
}
template <class S> size_t flops_tp_grads () { return 14; }

// non-singular targets
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S ood2 = my_recip(distsq);
  const S corefac = my_recip(sr*sr + tr*tr);
  const S reld2 = corefac / ood2;
  *r2 = exp_cond(ood2, corefac, reld2);
  *bbb = S(-2.0) * (*r2) * ood2
       + S(2.0) * corefac * ood2 * my_exp(-reld2);
}
template <class S> size_t flops_tv_grads () { return 16; }
#endif


#ifdef USE_WL_KERNEL
//
// Winckelmans-Leonard - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S s2 = sr*sr;
  const S d2 = distsq + s2;
  return (d2 + s2) / (d2*d2);
}
template <class S> size_t flops_tp_nograds () { return 5; }

// and the one for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S s2 = sr*sr + tr*tr;
  const S d2 = distsq + s2;
  return (d2 + s2) / (d2*d2);
}
template <class S> size_t flops_tv_nograds () { return 7; }

//
// Winckelmans-Leonard - with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S d2 = distsq + s2;
  const S ood2s = my_recip(d2*d2);
  *r2 = (d2 + s2) * ood2s;
  *bbb = (S(2.0) - S(4.0) * (*r2) * d2) * ood2s;
}
template <class S> size_t flops_tp_grads () { return 10; }

// and the one for non-singular targets
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S s2 = sr*sr + tr*tr;
  const S d2 = distsq + s2;
  const S ood2s = my_recip(d2*d2);
  *r2 = (d2 + s2) * ood2s;
  *bbb = (S(2.0) - S(4.0) * (*r2) * d2) * ood2s;
}
template <class S> size_t flops_tv_grads () { return 12; }
#endif


#ifdef USE_V2_KERNEL
//
// Vatistas n=2 - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S s2 = sr*sr;
  return my_rsqrt(distsq*distsq + s2*s2);
}
template <class S> size_t flops_tp_nograds () { return 6; }

// and for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S s2 = sr*sr;
  const S t2 = tr*tr;
  return my_rsqrt(distsq*distsq + s2*s2 + t2*t2);
}
template <class S> size_t flops_tv_nograds () { return 9; }

//
// Vatistas n=2 - with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  *r2 = my_rsqrt(distsq*distsq + s2*s2);
  *bbb = S(-2.0) * (*r2) * (*r2) * (*r2) * distsq;
}
template <class S> size_t flops_tp_grads () { return 10; }

// and the one for non-singular targets
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S t2 = tr*tr;
  *r2 = my_rsqrt(distsq*distsq + s2*s2 + t2*t2);
  *bbb = S(-2.0) * (*r2) * (*r2) * (*r2) * distsq;
}
template <class S> size_t flops_tv_grads () { return 13; }
#endif


#ifdef USE_V3_KERNEL
//
// Vatistas n=3 - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S s2 = sr*sr;
  return my_rcbrt(distsq*distsq*distsq + s2*s2*s2);
}
template <class S> size_t flops_tp_nograds () { return 8; }

// and for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S s2 = sr*sr;
  const S t2 = tr*tr;
  return my_rcbrt(distsq*distsq*distsq + s2*s2*s2 + t2*t2*t2);
}
template <class S> size_t flops_tv_nograds () { return 12; }

//
// Vatistas n=3 - with gradients
//
template <class S>
static inline void core_func (const S distsq, const S sr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S ds6 = distsq*distsq*distsq + s2*s2*s2;
  *r2 = my_rcbrt(ds6);
  *bbb = S(-2.0) * distsq*distsq * (*r2) / ds6;
}
template <class S> size_t flops_tp_grads () { return 12; }

// and the one for non-singular targets
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S s2 = sr*sr;
  const S t2 = tr*tr;
  const S ds6 = distsq*distsq*distsq + s2*s2*s2 + t2*t2*t2;
  *r2 = my_rcbrt(ds6);
  *bbb = S(-2.0) * distsq*distsq * (*r2) / ds6;
}
template <class S> size_t flops_tv_grads () { return 16; }
#endif

