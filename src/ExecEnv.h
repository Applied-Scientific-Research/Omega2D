/*
 * ExecEnv.h - Execution environment class
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <string>

// solver type/order
enum summation_t {
  direct    = 1,
  barneshut = 2,
  vic       = 3,	// unsupported internally
  fmm       = 4		// unsupported internally
};

// solver acceleration
enum accel_t {
  cpu_x86    = 1,
  cpu_vc     = 2,
  gpu_opengl = 3,	// unsupported internally
  gpu_cuda   = 4	// unsupported internally
};


//
// Class for the execution environment
//
class ExecEnv {
public:
  // primary constructor
  ExecEnv(const bool _internal,
          const summation_t _sumtype,
          const accel_t _acceltype)
    : m_internal(_internal),
      m_summ(_sumtype),
      m_accel(_acceltype)
    {}

  // default (delegating) ctor
  ExecEnv()
#ifdef EXTERNAL_VEL_SOLVE
  #ifdef USE_VC
    : ExecEnv(true, direct, cpu_vc)
  #else
    : ExecEnv(true, direct, cpu_x86)
  #endif
#else
  #ifdef USE_VC
    : ExecEnv(false, direct, cpu_vc)
  #else
    : ExecEnv(false, direct, cpu_x86)
  #endif
#endif
    {}

  void use_internal() { m_internal = true; };
  void use_external() { m_internal = false; };
  bool is_internal() const { return m_internal; };
  void set_summation(const summation_t _newsumm) { m_summ = _newsumm; };
  void set_instrs(const accel_t _newaccel) { m_accel = _newaccel; };
  accel_t get_instrs() const { return m_accel; };

  std::string to_string() const {
    std::string mystr;
    if (m_internal) {
      mystr += " external solver";
    } else {
      if (m_accel == cpu_x86) {
        mystr += " native";
      } else if (m_accel == cpu_vc) {
        mystr += " Vc-accelerated";
      } else {
        mystr += " unknown acceleration";
      }
      if (m_summ == direct) {
        mystr += " direct sums";
      } else if (m_summ == barneshut) {
        mystr += " treecode";
      } else {
        mystr += " unknown algorithm";
      }
    }
    return mystr;
  }

protected:
  bool m_internal;
  summation_t m_summ;
  accel_t m_accel;
};

