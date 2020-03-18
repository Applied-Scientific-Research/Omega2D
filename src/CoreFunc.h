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
#define USE_WL_KERNEL
//#define USE_V2_KERNEL


#ifdef USE_RM_KERNEL
//
// core functions - Rosenhead-Moore
//

template <class S> size_t flops_tv_nograds () { return 5; }

template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = distsq + sr*sr + tr*tr;
#ifdef USE_VC
  return Vc::reciprocal(r2);
#else
  return S(1.0) / r2;
#endif
}

template <class S> size_t flops_tp_nograds () { return 3; }

template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = distsq + sr*sr;
#ifdef USE_VC
  return Vc::reciprocal(r2);
#else
  return S(1.0) / r2;
#endif
}
#endif


#ifdef USE_EXPONENTIAL_KERNEL
//
// core functions - exponential
//

template <class S> size_t flops_tv_nograds () { return 9; }

template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
#ifdef USE_VC
  const S ood2 = Vc::reciprocal(distsq);
  const S corefac = Vc::reciprocal(sr*sr + tr*tr);
#else
  const S ood2 = S(1.0) / distsq;
  const S corefac = S(1.0) / (std::pow(sr,2) + std::pow(tr,2));
#endif
  const S reld2 = corefac / ood2;
#ifdef USE_VC
  S returnval = ood2;
  returnval(reld2 < S(16.0)) = ood2 * (S(1.0) - Vc::exp(-reld2));
  returnval(reld2 < S(0.001)) = corefac;
  return returnval;
#else
  if (reld2 > S(16.0)) {
    return ood2;
  } else if (reld2 < S(0.001)) {
    return corefac;
  } else {
    return ood2 * (S(1.0) - std::exp(-reld2));
  }
#endif
}

template <class S> size_t flops_tp_nograds () { return 7; }

template <class S>
static inline S core_func (const S distsq, const S sr) {
#ifdef USE_VC
  const S ood2 = Vc::reciprocal(distsq);
  const S corefac = Vc::reciprocal(sr*sr);
#else
  const S ood2 = S(1.0) / distsq;
  const S corefac = S(1.0) / std::pow(sr,2);
#endif
  const S reld2 = corefac / ood2;
  // 4 flops to here
#ifdef USE_VC
  S returnval = ood2;
  returnval(reld2 < S(16.0)) = ood2 * (S(1.0) - Vc::exp(-reld2));
  returnval(reld2 < S(0.001)) = corefac;
  //std::cout << std::endl << dist << std::endl << reld2 << std::endl << returnval << std::endl;
  return returnval;
#else
  if (reld2 > S(16.0)) {
    return ood2;
    // 1 flop (comparison)
  } else if (reld2 < S(0.001)) {
    return corefac;
    // 2 flops
  } else {
    return ood2 * (S(1.0) - std::exp(-reld2));
    // 3 flops
  }
#endif
}
#endif


#ifdef USE_WL_KERNEL
//
// core functions - Winckelmans–Leonard
//

template <class S> size_t flops_tv_nograds () { return 8; }

template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr + tr*tr;
  const S d2 = distsq + r2;
  return (distsq + S(2.0)*r2) / (d2*d2);
}

template <class S> size_t flops_tp_nograds () { return 6; }

template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S r2 = sr*sr;
  const S d2 = distsq + r2;
  return (distsq + S(2.0)*r2) / (d2*d2);
}
#endif


#ifdef USE_V2_KERNEL
//
// core functions - Vatistas n=2
//

template <class S> size_t flops_tv_nograds () { return 7; }

template <class S>
static inline S core_func (const S distsq, const S sr, const S tr) {
  const S d2 = distsq + sr*sr + tr*tr;
#ifdef USE_VC
  return Vc::rsqrt(d2*d2);
#else
  return S(1.0) / std::sqrt(d2*d2);
#endif
}

template <class S> size_t flops_tp_nograds () { return 5; }

template <class S>
static inline S core_func (const S distsq, const S sr) {
  const S d2 = distsq + sr*sr;
#ifdef USE_VC
  return Vc::rsqrt(d2*d2);
#else
  return S(1.0) / std::sqrt(d2*d2);
#endif
}
#endif

