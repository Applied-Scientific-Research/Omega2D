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


#ifdef USE_RM_KERNEL
//
// Rosenhead-Moore velocity-only
//
#ifdef USE_VC
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = distsq + sr*sr;
  return Vc::reciprocal(r2);
}
#endif
template <>
inline float core_func (const float distsq, const float sr) {
  const float r2 = distsq + sr*sr;
  return 1.0f / r2;
}
template <>
inline double core_func (const double distsq, const double sr) {
  const double r2 = distsq + sr*sr;
  return 1.0 / r2;
}
template <class S> size_t flops_tp_nograds () { return 3; }

// and the one for non-singular targets
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  return core_func(distsq + tr*tr, sr);
}
template <class S> size_t flops_tv_nograds () { return 5; }

//
// Rosenhead-Moore with gradients
//
#ifdef USE_VC
template <class S>
static inline void core_func (const S distsq, const S sr, const S tr,
                              S* const __restrict__ r2, S* const __restrict__ bbb) {
  const S td2 = distsq + sr*sr + tr*tr;
  *r2 = Vc::reciprocal(td2);
  *bbb = S(-2.0) * (*r2) * (*r2);
}
#endif
template <>
inline void core_func (const float distsq, const float sr, const float tr,
                       float* const __restrict__ r2, float* const __restrict__ bbb) {
  const float td2 = distsq + sr*sr + tr*tr;
  *r2 = 1.0f / td2;
  *bbb = -2.0f * (*r2) * (*r2);
}
template <>
inline void core_func (const double distsq, const double sr, const double tr,
                       double* const __restrict__ r2, double* const __restrict__ bbb) {
  const double td2 = distsq + sr*sr + tr*tr;
  *r2 = 1.0 / td2;
  *bbb = -2.0 * (*r2) * (*r2);
}
template <class S> size_t flops_tv_grads () { return 9; }
#endif


#ifdef USE_EXPONENTIAL_KERNEL
//
// exponential core - velocity only
//
#ifdef USE_VC
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S ood2 = Vc::reciprocal(distsq);
  const S corefac = Vc::reciprocal(sr*sr);
  const S reld2 = corefac / ood2;
  // 4 flops to here
  S returnval = ood2;
  returnval(reld2 < S(16.0)) = ood2 * (S(1.0) - Vc::exp(-reld2));
  returnval(reld2 < S(0.001)) = corefac;
  return returnval;
}
#endif
template <>
inline float core_func (const float distsq, const float sr) {
  const float ood2 = 1.0f / distsq;
  const float corefac = 1.0f / std::pow(sr,2);
  const float reld2 = corefac / ood2;
  // 4 flops to here
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
template <>
inline double core_func (const double distsq, const double sr) {
  const double ood2 = 1.0 / distsq;
  const double corefac = 1.0 / std::pow(sr,2);
  const double reld2 = corefac / ood2;
  // 4 flops to here
  if (reld2 > 16.0) {
    return ood2;
    // 1 flop (comparison)
  } else if (reld2 < 0.001) {
    return corefac;
    // 2 flops
  } else {
    return ood2 * (1.0 - std::exp(-reld2));
    // 3 flops
  }
}
template <class S> size_t flops_tp_nograds () { return 7; }

// non-singular targets
#ifdef USE_VC
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S ood2 = Vc::reciprocal(distsq);
  const S corefac = Vc::reciprocal(sr*sr + tr*tr);
  const S reld2 = corefac / ood2;
  S returnval = ood2;
  returnval(reld2 < S(16.0)) = ood2 * (S(1.0) - Vc::exp(-reld2));
  returnval(reld2 < S(0.001)) = corefac;
  return returnval;
}
#endif
template <>
inline float core_func (const float distsq, const float sr, const float tr) {
  const float ood2 = 1.0f / distsq;
  const float corefac = 1.0f / (std::pow(sr,2) + std::pow(tr,2));
  const float reld2 = corefac / ood2;
  if (reld2 > 16.0f) {
    return ood2;
  } else if (reld2 < 0.001f) {
    return corefac;
  } else {
    return ood2 * (1.0f - std::exp(-reld2));
  }
}
template <>
inline double core_func (const double distsq, const double sr, const double tr) {
  const double ood2 = 1.0 / distsq;
  const double corefac = 1.0 / (std::pow(sr,2) + std::pow(tr,2));
  const double reld2 = corefac / ood2;
  if (reld2 > 16.0) {
    return ood2;
  } else if (reld2 < 0.001) {
    return corefac;
  } else {
    return ood2 * (1.0 - std::exp(-reld2));
  }
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
#ifdef USE_VC
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  return Vc::rsqrt(distsq*distsq + r2*r2);
}
template <>
inline float core_func (const float distsq, const float sr) {
  const float r2 = sr*sr;
  return 1.0f / std::sqrt(distsq*distsq + r2*r2);
}
template <>
inline double core_func (const double distsq, const double sr) {
  const double r2 = sr*sr;
  // std::hypot was slower
  return 1.0 / std::sqrt(distsq*distsq + r2*r2);
}
#else
template <class S>
inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  return S(1.0) / std::sqrt(distsq*distsq + r2*r2);
}
#endif
template <class S> size_t flops_tp_nograds () { return 6; }

// and for non-singular targets
#ifdef USE_VC
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr;
  const S o2 = tr*tr;
  return Vc::rsqrt(distsq*distsq + r2*r2 + o2*o2);
}
template <>
inline float core_func (const float distsq, const float sr, const float tr) {
  const float r2 = sr*sr;
  const float o2 = tr*tr;
  return 1.0f / std::sqrt(distsq*distsq + r2*r2 + o2*o2);
}
template <>
inline double core_func (const double distsq, const double sr, const double tr) {
  const double r2 = sr*sr;
  const double o2 = tr*tr;
  return 1.0 / std::sqrt(distsq*distsq + r2*r2 + o2*o2);
}
#else
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr;
  const S o2 = tr*tr;
  return 1.0f / std::sqrt(distsq*distsq + r2*r2 + o2*o2);
}
#endif
template <class S> size_t flops_tv_nograds () { return 9; }
#endif


#ifdef USE_V3_KERNEL
//
// Vatistas n=3 - velocity only
//
#ifdef USE_VC
template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  return Vc::exp(S(-0.3333333)*Vc::log(distsq*distsq*distsq + r2*r2*r2));
}
#endif
template <>
inline float core_func (const float distsq, const float sr) {
  const float r2 = sr*sr;
  return 1.0f / std::cbrt(distsq*distsq*distsq + r2*r2*r2);
}
template <>
inline double core_func (const double distsq, const double sr) {
  const double r2 = sr*sr;
  return 1.0 / std::cbrt(distsq*distsq*distsq + r2*r2*r2);
}
template <class S> size_t flops_tp_nograds () { return 8; }

// and for non-singular targets
#ifdef USE_VC
template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr;
  const S o2 = tr*tr;
  return Vc::exp(S(-0.3333333)*Vc::log(distsq*distsq*distsq + r2*r2*r2 + o2*o2*o2));
}
#endif
template <>
inline float core_func (const float distsq, const float sr, const float tr) {
  const float r2 = sr*sr;
  const float o2 = tr*tr;
  return 1.0f / std::cbrt(distsq*distsq*distsq + r2*r2*r2 + o2*o2*o2);
}
template <>
inline double core_func (const double distsq, const double sr, const double tr) {
  const double r2 = sr*sr;
  const double o2 = tr*tr;
  return 1.0 / std::cbrt(distsq*distsq*distsq + r2*r2*r2 + o2*o2*o2);
}
template <class S> size_t flops_tv_nograds () { return 12; }
#endif

