/*
 * ExecEnv.h - Execution environment class
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once


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

  // delegating ctor
  ExecEnv()
    : ExecEnv(true, direct, cpu_x86)
    {
#ifdef EXTERNAL_VEL_SOLVE
      m_internal = false;
#endif
#ifdef USE_VC
      m_accel = cpu_vc;
#endif
    }

  void use_internal() { m_internal = true; };
  void use_external() { m_internal = false; };
  bool is_internal() const { return m_internal; };
  void set_summation(const summation_t _newsumm) { m_summ = _newsumm; };
  void set_instrs(const accel_t _newaccel) { m_accel = _newaccel; };
  accel_t get_instrs() const { return m_accel; };

protected:
  bool m_internal;
  summation_t m_summ;
  accel_t m_accel;
};

