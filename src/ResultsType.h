/*
 * ResultsType.h - Solution types from influence calculations
 *
 * (c)2020-1 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <string>

// solver type/order
enum results_t {
  velonly     = 1,
  velandgrad  = 2,
  psionly     = 3,
  velandvort  = 4,
  velandshear = 5
};

//
// Class for the execution environment
//
class ResultsType {
public:
  // primary constructor
  ResultsType(const results_t _type)
    : m_rtype(_type)
    {}

  // default (delegating) ctor, solve for vels only
  ResultsType()
    : ResultsType(velonly)
    {}

  // accept list of results types, pick appropriate enum value
  ResultsType(const bool _psi,
              const bool _vel,
              const bool _grad,
              const bool _vort) {
    if (_vel and not _psi and not _grad and not _vort) m_rtype = velonly;
    else if (_vel and not _psi and _grad and not _vort) m_rtype = velandgrad;
    else if (_vel and not _psi and not _grad and not _vort) m_rtype = velonly;
    else if (not _vel and _psi and not _grad and not _vort) m_rtype = psionly;
    else if (_vel and not _psi and not _grad and _vort) m_rtype = velandvort;
    else {
      // what do we do here?
    }
  }

  // streamfunction (1-component)
  const bool compute_psi() const { return m_rtype == psionly; }
  // velocities (2-component)
  const bool compute_vel() const { return m_rtype != psionly; }
  // velocity gradients (2x2 matrix)
  const bool compute_grad() const { return m_rtype == velandgrad; }
  // vorticity (1-component)
  const bool compute_vort() const { return m_rtype == velandvort; }
  // shear rate (1-component)
  const bool compute_shear() const { return m_rtype == velandshear; }

  const results_t get_type() const { return m_rtype; }

  std::string to_string() {
    return to_string(m_rtype);
  }

  std::string to_string(results_t _type) {
    std::string mystr;
    if (_type == psionly) mystr += " for psi";
    else if (_type == velonly) mystr += " for vel";
    else if (_type == velandgrad) mystr += " for vel and grad";
    else if (_type == velandvort) mystr += " for vel and vort";
    else if (_type == velandshear) mystr += " for vel and shear";
    return mystr;
  }

protected:
  results_t m_rtype;
};

