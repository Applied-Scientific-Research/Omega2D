/*
 * ExecEnv.h - Execution environment class
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
    : ExecEnv(false, direct, cpu_vc)
  #else
    : ExecEnv(false, direct, cpu_x86)
  #endif
#else
  #ifdef USE_VC
    : ExecEnv(true, direct, cpu_vc)
  #else
    : ExecEnv(true, direct, cpu_x86)
  #endif
#endif
    {}

  void use_internal() { m_internal = true; };
  void use_external() { m_internal = false; };
  void set_internal(const bool _isint) { m_internal = _isint; };
  bool is_internal() const { return m_internal; };
  void set_summation(const summation_t _newsumm) { m_summ = _newsumm; };
  void set_instrs(const accel_t _newaccel) { m_accel = _newaccel; };
  accel_t get_instrs() const { return m_accel; };

  std::string to_string() const {
    std::string mystr;
    if (m_internal) {
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
    } else {
      mystr += " external solver";
    }
    return mystr;
  }

protected:
  bool m_internal;
  summation_t m_summ;
  accel_t m_accel;
};

