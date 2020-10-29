/*
 * MathHelper.h - Non-class helper inlines for core functions, kernels, influences
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#ifdef _WIN32
#define __restrict__ __restrict
#endif

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <cmath>

// helper functions: recip, rsqrt, rcbrt, acos

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
template <class S>
static inline S my_halflog (const S _x) {
  return S(0.5f) * Vc::log(_x);
}
template <> inline float my_halflog (const float _x) {
  return 0.5f * std::log(_x);
}
template <> inline double my_halflog (const double _x) {
  return 0.5 * std::log(_x);
}
#else
template <class S>
static inline S my_halflog (const S _x) {
  return 0.5f * std::log(_x);
}
#endif
