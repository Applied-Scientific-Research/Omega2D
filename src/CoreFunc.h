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

//
// core functions - Rosenhead-Moore
//

// this is always 5 flops
template <class S>
static inline S core_rm (const S distsq, const S sr, const S tr) {
  const S r2 = distsq + sr*sr + tr*tr;
#ifdef USE_VC
  return Vc::reciprocal(r2);
#else
  return S(1.0) / r2;
#endif
}

// this is always 3 flops
template <class S>
static inline S core_rm (const S distsq, const S sr) {
  const S r2 = distsq + sr*sr;
#ifdef USE_VC
  return Vc::reciprocal(r2);
#else
  return S(1.0) / r2;
#endif
}

//
// core functions - exponential
//

// this probably averages out to 9 flops
template <class S>
static inline S core_exp (const S distsq, const S sr, const S tr) {
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

// this probably averages out to 7 flops
template <class S>
static inline S core_exp (const S distsq, const S sr) {
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

//
// core functions - Winckelmans–Leonard
//

// this is always 8 flops
template <class S>
static inline S core_wl (const S distsq, const S sr, const S tr) {
  const S r2 = sr*sr + tr*tr;
  const S d2 = distsq + r2;
  return (distsq + S(2.0)*r2) / (d2*d2);
}

// this is always 6 flops
template <class S>
static inline S core_wl (const S distsq, const S sr) {
  const S r2 = sr*sr;
  const S d2 = distsq + r2;
  return (distsq + S(2.0)*r2) / (d2*d2);
}

//
// core functions - Vatistas n=2
//

// this is always 7 flops
template <class S>
static inline S core_v2 (const S distsq, const S sr, const S tr) {
  const S d2 = distsq + sr*sr + tr*tr;
#ifdef USE_VC
  return Vc::rsqrt(d2*d2);
#else
  return S(1.0) / std::sqrt(d2*d2);
#endif
}

// this is always 5 flops
template <class S>
static inline S core_v2 (const S distsq, const S sr) {
  const S d2 = distsq + sr*sr;
#ifdef USE_VC
  return Vc::rsqrt(d2*d2);
#else
  return S(1.0) / std::sqrt(d2*d2);
#endif
}

