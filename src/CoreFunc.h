/*
 * CoreFunc.h - Non-class core function inlines for influence calculations
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>

//#define USE_RM_KERNEL
//#define USE_EXPONENTIAL_KERNEL
//#define USE_WL_KERNEL
#define USE_V2_KERNEL
//#define USE_V3_KERNEL


// helper functions: recip, rsqrt, rcbrt

#ifdef USE_VC
template <class S>
static inline S my_recip(const S _in) {
  return Vc::reciprocal(_in);
}
template <>
inline float my_recip(const float _in) {
  return 1.0f / _in;
}
template <>
inline double my_recip(const double _in) {
  return 1.0 / _in;
}
#else
template <class S>
static inline S my_recip(const S _in) {
  return S(1.0) / _in;
}
#endif

#ifdef USE_VC
template <class S>
static inline S my_rsqrt(const S _in) {
  return Vc::rsqrt(_in);
}
template <>
inline float my_rsqrt(const float _in) {
  return 1.0f / std::sqrt(_in);
}
template <>
inline double my_rsqrt(const double _in) {
  return 1.0 / std::sqrt(_in);
}
#else
template <class S>
static inline S my_rsqrt(const S _in) {
  return S(1.0) / std::sqrt(_in);
}
#endif

#ifdef USE_VC
template <class S>
static inline S my_rcbrt(const S _in) {
  return Vc::exp(S(-0.3333333)*Vc::log(_in));
}
template <>
inline float my_rcbrt(const float _in) {
  return 1.0f / std::cbrt(_in);
}
template <>
inline double my_rcbrt(const double _in) {
  return 1.0 / std::cbrt(_in);
}
#else
template <class S>
static inline S my_rcbrt(const S _in) {
  return S(1.0) / std::cbrt(_in);
}
#endif


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
  const S ood2 = Vc::reciprocal(distsq);
  const S corefac = Vc::reciprocal(sr*sr + tr*tr);
  const S reld2 = corefac / ood2;
  return exp_cond(ood2, corefac, reld2);
}
template <class S> size_t flops_tv_nograds () { return 9; }

//
// exponential core - with gradients
//
// not done
#endif


#ifdef USE_WL_KERNEL
//
// Winckelmans–Leonard - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  const S d2 = distsq + r2;
  return (d2 + r2) / (d2*d2);
}
template <class S> size_t flops_tp_nograds () { return 5; }

// and the one for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr + tr*tr;
  const S d2 = distsq + r2;
  return (d2 + r2) / (d2*d2);
}
template <class S> size_t flops_tv_nograds () { return 7; }
#endif


#ifdef USE_V2_KERNEL
//
// Vatistas n=2 - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  return my_rsqrt(distsq*distsq + r2*r2);
}
template <class S> size_t flops_tp_nograds () { return 6; }

// and for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr;
  const S o2 = tr*tr;
  return my_rsqrt(distsq*distsq + r2*r2 + o2*o2);
}
template <class S> size_t flops_tv_nograds () { return 9; }
#endif


#ifdef USE_V3_KERNEL
//
// Vatistas n=3 - velocity only
//
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  return my_rcbrt(distsq*distsq*distsq + r2*r2*r2);
}
template <class S> size_t flops_tp_nograds () { return 8; }

// and for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr;
  const S o2 = tr*tr;
  return my_rcbrt(distsq*distsq*distsq + r2*r2*r2 + o2*o2*o2);
}
template <class S> size_t flops_tv_nograds () { return 12; }
#endif

