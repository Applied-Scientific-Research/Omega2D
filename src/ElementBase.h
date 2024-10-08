/*
 * ElementBase.h - abstract class for arrays of any computational elements
 *
 * (c)2018-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

#include "VectorHelper.h"
#include "Omega2D.h"
#include "Body.h"
#include "ElementPacket.h"

#include <iostream>
#include <vector>
#include <memory>
#include <cassert>
#include <optional>
#include <variant>
#include <algorithm>
#include <cmath>


// the superclass

template <class S>
class ElementBase {
public:
  ElementBase<S>(const size_t _n,
                 const elem_t _e,
                 const move_t _m,
                 std::shared_ptr<Body> _bp) :
      E(_e), M(_m), B(_bp), n(_n) {
  }

  size_t get_n() const { return n; }
  bool is_inert() const { return E==inert; }
  elem_t                                  get_elemt() const { return E; }
  move_t                                  get_movet() const { return M; }
  const std::shared_ptr<Body>             get_body_ptr() const { return B; }
  std::shared_ptr<Body>                   get_body_ptr()   { return B; }
  const std::array<Vector<S>,Dimensions>& get_pos() const  { return x; }
  std::array<Vector<S>,Dimensions>&       get_pos()        { return x; }
  const Vector<S>&                        get_str() const  { return *s; }
  Vector<S>&                              get_str()        { return *s; }
  const std::array<Vector<S>,Dimensions>& get_vel() const  { return u; }
  std::array<Vector<S>,Dimensions>&       get_vel()        { return u; }

  const bool has_vort() const { return (bool)w; }
  const Vector<S>& get_vort() const {
    // can't change this object, so return what we have
    return *w;
  }
  Vector<S>& get_vort() {
    // allowed to change this object, to resize or create the Vector
    if (w) {
      // check size before returning
      if (w->size() != n) w->resize(n);
    } else {
      // allocate and zero before returning
      Vector<S> new_vort;
      new_vort.resize(n);
      w = std::move(new_vort);
      std::fill((*w).begin(), (*w).end(), 0.0);
    }
    return *w;
  }
  void set_vort(Vector<S> _in) {
    assert(_in.size() == n && "ERROR (ElementBase::set_vort) input vector length mismatch");
    if (w) {
      if (w->size() != _in.size()) w->resize(n);
    } else {
      // allocate and move
      Vector<S> new_vort;
      new_vort.resize(n);
      w = std::move(new_vort);
    }
    *w = _in;
    //std::cout << "Received vorticity on " << n << " nodes, starting with " << (*w)[0] << std::endl;
  }

  const bool has_shear() const { return (bool)e; }
  const Vector<S>& get_shear() const {
    // can't change this object, so return what we have
    return *e;
  }
  Vector<S>& get_shear() {
    // allowed to change this object, to resize or create the Vector
    if (e) {
      // check size before returning
      if (e->size() != n) e->resize(n);
    } else {
      // allocate and zero before returning
      Vector<S> new_shear;
      new_shear.resize(n);
      e = std::move(new_shear);
      std::fill((*e).begin(), (*e).end(), 0.0);
    }
    return *e;
  }

  void set_str(const size_t ioffset, const size_t icnt, Vector<S> _in) {
    assert(s && "Strength array does not exist");
    assert(_in.size() == (*s).size() && "Set strength array size does not match");

    // copy over the strengths
    *s = _in;
  }

  // child class calls here to add nodes and other properties
  void add_new(const ElementPacket<float>& _in) {

    // check inputs
    if (_in.x.size() == 0 or _in.nelem == 0) return;
    assert(_in.x.size() % Dimensions == 0 && "Input ElementPacket does not have correct number of node coords");

    // new number of nodes (not elements)
    const size_t nnew = _in.x.size()/Dimensions;

    // add node coordinates
    for (size_t d=0; d<Dimensions; ++d) {
      // extend with more space for new values
      x[d].resize(n+nnew);
      // copy new values to end of vector
      for (size_t i=0; i<nnew; ++i) {
        x[d][n+i] = _in.x[Dimensions*i+d];
      }
    }

    // strength - ASSUMPTION: the strength is the first of each val entry
    if (s) {
      assert(_in.val.size() >= nnew && "Input ElementPacket does not have enough values in val");
      const size_t nper = _in.val.size() / nnew;
      // must dereference s to get the actual vector
      (*s).resize(n+nnew);
      for (size_t i=0; i<nnew; ++i) {
        (*s)[n+i] = _in.val[nper*i+0];
      }
    }

    // extend the other vectors as well
    for (size_t d=0; d<Dimensions; ++d) {
      u[d].resize(n+nnew);
    }

    // finally, update n
    n += nnew;
  }


  // up-size all arrays to the new size, filling with sane values
  // this only happens right after diffusion
  void resize(const size_t _nnew) {
    const size_t currn = n;
    //std::cout << "  inside ElementBase::resize with " << currn << " " << _nnew << std::endl;
    if (_nnew == currn) return;

    // positions first
    for (size_t d=0; d<Dimensions; ++d) {
      const size_t thisn = x[d].size();
      x[d].resize(_nnew);
      for (size_t i=thisn; i<_nnew; ++i) {
        x[d][i] = 0.0;
      }
    }

    // strength
    if (s) {
      const size_t thisn = (*s).size();
      (*s).resize(_nnew);
      for (size_t i=thisn; i<_nnew; ++i) {
        (*s)[i] = 0.0;
      }
    }

    // and finally velocity (no need to set it)
    for (size_t d=0; d<Dimensions; ++d) {
      u[d].resize(_nnew);
    }

    // lastly, update n
    n = _nnew;
  }

  // should rename these zero_results
  void zero_vels() {
    for (size_t d=0; d<Dimensions; ++d) {
      std::fill(u[d].begin(), u[d].end(), 0.0);
    }
    if (w) {
      if (w->size() != n) w->resize(n);
      std::fill((*w).begin(), (*w).end(), 0.0);
    }
    if (e) {
      if (e->size() != n) e->resize(n);
      std::fill((*e).begin(), (*e).end(), 0.0);
    }
  }

  // should rename these finalize_results
  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    const double factor = 0.5/M_PI;
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<get_n(); ++i) {
        u[d][i] = _fs[d] + u[d][i] * factor;
      }
    }
    if (w) {
      assert(w->size() == n && "Vorticity vector not sized properly in finalize_vels");
      for (size_t i=0; i<get_n(); ++i) {
        (*w)[i] *= factor;
      }
    }
    if (e) {
      assert(e->size() == n && "Shear rate vector not sized properly in finalize_vels");
      for (size_t i=0; i<get_n(); ++i) {
        (*e)[i] *= factor;
      }
    }
  }

  void add_body_motion(const S factor, const double _time) {
    // do nothing here
  }

  void zero_strengths() {
    if (s) {
      std::fill((*s).begin(), (*s).end(), 0.0);
    }
  }

  // do nothing here
  void add_unit_rot_strengths() {}
  void add_solved_rot_strengths(const S _factor) {}

  void transform(const double _time) {
    // reset positions according to prescribed motion
    if (B and M == bodybound) {

      // tell the Body to compute and save its position, vel, angular pos and angular vel
      B->transform(_time);

      // for the no-rotation case, we can just transform here
      std::array<double,Dimensions> thispos = B->get_pos();
      const double theta = B->get_orient();
      const S st = std::sin(theta);
      const S ct = std::cos(theta);

      std::cout << "    transforming body at time " << (S)_time << " to " << (S)thispos[0] << " " << (S)thispos[1]
                << " and theta " << theta << " omega " << B->get_rotvel() << std::endl;

      // and do the transform
      for (size_t i=0; i<get_n(); ++i) {
        // rotate and translate
        x[0][i] = (S)thispos[0] + (*ux)[0][i]*ct - (*ux)[1][i]*st;
        x[1][i] = (S)thispos[1] + (*ux)[0][i]*st + (*ux)[1][i]*ct;
      }
    }
  }

  // time is the starting time, time+dt is the ending time
  void move(const double _time, const double _dt,
            const double _wt1, ElementBase<S> const & _u1) {

    if (M == lagrangian) {
      std::cout << "  Moving" << to_string() << std::endl;

      // we do not need to update vels because _u1 is the same as this

      // update positions
      for (size_t d=0; d<Dimensions; ++d) {
        for (size_t i=0; i<n; ++i) {
          x[d][i] += (S)_dt * _wt1*_u1.u[d][i];
        }
      }

      // update strengths (in derived class)

    } else if (B and M == bodybound) {
      transform(_time+_dt);
    }
  }

  // time is the starting time, time+dt is the ending time
  void move(const double _time, const double _dt,
            const double _wt1, ElementBase<S> const & _u1,
            const double _wt2, ElementBase<S> const & _u2) {

    // must confirm that incoming time derivates include velocity
    // if this has vels, then lets advect it
    if (M == lagrangian) {
      std::cout << "  Moving" << to_string() << std::endl;

      // update vels, note that _u1 is the same as this
      for (size_t d=0; d<Dimensions; ++d) {
        for (size_t i=0; i<n; ++i) {
          u[d][i] = _wt1*_u1.u[d][i] + _wt2*_u2.u[d][i];
        }
      }

      // update positions
      for (size_t d=0; d<Dimensions; ++d) {
        for (size_t i=0; i<n; ++i) {
          x[d][i] += (S)_dt * u[d][i];
        }
      }

      // update strengths (in derived class)

    } else if (B and M == bodybound) {
      transform(_time+_dt);
    }
  }

  // time is the starting time, time+dt is the ending time
  void move(const double _time, const double _dt,
            const double _wt1, ElementBase<S> const & _u1,
            const double _wt2, ElementBase<S> const & _u2,
            const double _wt3, ElementBase<S> const & _u3) {

    // must confirm that incoming time derivates include velocity
    // if this has vels, then lets advect it
    if (M == lagrangian) {
      std::cout << "  Moving" << to_string() << std::endl;

      // update vels, note that _u1 is the same as this
      for (size_t d=0; d<Dimensions; ++d) {
        for (size_t i=0; i<n; ++i) {
          u[d][i] = _wt1*_u1.u[d][i] + _wt2*_u2.u[d][i] + _wt3*_u3.u[d][i];
        }
      }

      // update positions
      for (size_t d=0; d<Dimensions; ++d) {
        for (size_t i=0; i<n; ++i) {
          x[d][i] += (S)_dt * u[d][i];
        }
      }

      // update strengths (in derived class)

    } else if (B and M == bodybound) {
      transform(_time+_dt);
    }
  }

  // find the new peak strength magnitude
  S get_max_str() {
    if (s) {
      // we have strengths, go through and check them
      const S this_max = *std::max_element(std::begin(*s), std::end(*s));
      const S this_min = *std::min_element(std::begin(*s), std::end(*s));
      return std::max(this_max, -this_min);
    } else {
      return 1.0;
    }
  }

  // add and return the total circulation of all elements
  S get_total_circ(const double _time) {
    S circ = 0.0;

    if (s) {
      if (s->size() < 40000) {
        // we have strengths, add them up
        // this is the c++17 way
        //return std::reduce(std::execution::par, s->begin(), s->end());
        // this is the c++11 way
        circ = std::accumulate(s->begin(), s->end(), 0.0);
      } else {
        // do it in double precision instead
        Vector<double> dblstr(s->begin(), s->end());
        double dblcirc = std::accumulate(dblstr.begin(), dblstr.end(), 0.0);
        circ = (float)dblcirc;
      }
    }

    return circ;
  }

  // add and return the total circulation of the rotating body volume
  //   will be specialized by Surface
  S get_body_circ(const double _time) {
    return 0.0;
  }

  std::string to_string() const {
    std::string mystr;
    if (E == active) {
      mystr += " Active";
    } else if (E == reactive) {
      mystr += " Reactive";
    } else {
      mystr += " Inert";
    }
    if (M == lagrangian) {
      mystr += " Lagrangian";
    } else if (M == bodybound) {
      mystr += " Body-fixed";
    } else {
      mystr += " Fixed";
    }
    return mystr;
  }

protected:
  // if you add anything here, you need to wipe out all build files and run cmake again!

  // active, reactive, or inert?
  elem_t E;
  // how does it move? use move_t and Body*
  move_t M;
  // if attached to a body, which one?
  std::shared_ptr<Body> B;

  // common arrays for all derived types
  size_t n;						// number of nodes

  // state vector
  std::array<Vector<S>,Dimensions> x;                   // position of nodes
  std::optional<Vector<S>> s;                           // strength at nodes

  // time derivative of state vector
  std::array<Vector<S>,Dimensions> u;                   // velocity at nodes
  std::optional<Vector<S>> w;                           // vorticity at nodes
  std::optional<Vector<S>> e;                           // shear rate at nodes
  //std::optional<std::array<Vector<S>,Dimensions>> dsdt; // strength change

  // for objects moving with a body
  std::optional<std::array<Vector<S>,Dimensions>> ux;   // untransformed position of nodes
};

