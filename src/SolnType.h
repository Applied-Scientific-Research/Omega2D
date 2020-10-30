/*
 * SolnType.h - Solution types from influence calculations
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <string>

// solver type/order
enum solution_t {
  velonly    = 1,
  velandgrad = 2,
  psionly    = 3,
  velandvort = 4
};

//
// Class for the execution environment
//
class SolnType {
public:
  // primary constructor
  SolnType(const solution_t _type)
    : m_stype(_type)
    {}

  // default (delegating) ctor, solve for vels only
  SolnType()
    : SolnType(velonly)
    {}

  // accept list of results types, pick appropriate enum value
  SolnType(const bool _psi,
           const bool _vel,
           const bool _grad,
           const bool _vort) {
    if (_vel and not _psi and not _grad and not _vort) m_stype = velonly;
    else if (_vel and not _psi and _grad and not _vort) m_stype = velandgrad;
    else if (_vel and not _psi and not _grad and not _vort) m_stype = velonly;
    else if (not _vel and _psi and not _grad and not _vort) m_stype = psionly;
    else if (_vel and not _psi and not _grad and _vort) m_stype = velandvort;
    else {
      // what do we do here?
    }
  }

  // streamfunction (1-component)
  const bool compute_psi() const { return m_stype == psionly; }
  // velocities (2-component)
  const bool compute_vel() const { return m_stype != psionly; }
  // velocity gradients (2x2 matrix)
  const bool compute_grad() const { return m_stype == velandgrad; }
  // vorticity (1-component)
  const bool compute_vort() const { return m_stype == velandvort; }
  const solution_t get_soln_type() const { return m_stype; }

  std::string to_string() const {
    return to_string(m_stype);
  }

  std::string to_string(const solution_t _type) {
    std::string mystr;
    if (_type == psionly) mystr += " for psi";
    else if (_type == velonly) mystr += " for vel";
    else if (_type == velandgrad) mystr += " for vel and grad";
    else if (_type == velandvort) mystr += " for vel and vort";
    return mystr;
  }

protected:
  solution_t m_stype;
};

