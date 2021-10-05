/*
 * Surfaces.h - Specialized class for surfaces in 2D
 *
 * (c)2019-21 Applied Scientific Research, Inc.
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

#include "Omega2D.h"
#include "VectorHelper.h"
#include "ElementBase.h"

#ifdef USE_GL
#include "GlState.h"
#include "RenderParams.h"
#include "OglHelper.h"
#include "ShaderHelper.h"
#include "glad/glad.h"
#endif

#include <iostream>
#include <iomanip> // for setfill and setw
#include <vector>
#include <array>
#include <algorithm> // for max_element
#include <optional>
#include <cassert>


// useful structure for panel strengths and BCs
// Strength<S> ps;
//   ps[0] means that the tangential (vortex) strength is present
//   ps[1] means that the normal (source) strength is present
template <class S> using Strength = std::array<std::optional<Vector<S>>,2>;

// useful structure for basis vectors
// Basis<S> b;
//   b[0] is the pair of arrays of the tangential normalized vectors, b[0][0] for x, b[0][1] for y
//   b[1] is the same for the normal vectors (+ is into fluid), b[1][0] for x, b[1][1] for y
template <class S> using Basis = std::array<std::array<Vector<S>,Dimensions>,Dimensions>;


// 1-D elements
template <class S>
class Surfaces: public ElementBase<S> {
public:
  // constructor - accepts standard ElementPacket now
  //               last parameter (_val) is either fixed strength or BC for each panel
  Surfaces(const ElementPacket<float>& _elems,
           const elem_t _e,
           const move_t _m,
           std::shared_ptr<Body> _bp)
    : ElementBase<S>(0, _e, _m, _bp),
      np(0),
      source_str_is_unknown(true),
      vol(-1.0),
      solved_omega(0.0),
      omega_error(0.0),
      this_omega(0.0),
      reabsorbed_gamma(0.0),
      max_strength(-1.0) {

    // pull out the vectors
    const std::vector<float>& _x   = _elems.x;
    const std::vector<Int>&   _idx = _elems.idx;
    const std::vector<float>& _val = _elems.val;

    // make sure input arrays are correctly-sized
    assert(_elems.ndim == 1 && "Input elements are not (1D) surfaces");
    const size_t nsurfs = _elems.nelem;

    // always initialize the ps panel strength optionals
    if (this->E != inert) {
      for (size_t d=0; d<2; ++d) {
        Vector<S> new_s;
        ps[d] = std::move(new_s);
       }
    }

    // if no surfs, quit out now
    if (nsurfs == 0) {
      // but still initialize ux before we go (in case first bfeature is not enabled)
      if (_bp) this->ux = this->x;
      return;
    }

    std::cout << "Surfaces sizes: _x=" << _x.size() << " _idx=" << _idx.size() << " _val=" << _val.size() << std::endl;
    assert(_val.size() % nsurfs == 0 && "Value array is not an even multiple of panel count");
    assert(_x.size() % Dimensions == 0 && "Position array is not an even multiple of dimensions");
    const size_t nnodes = _x.size() / Dimensions;

    std::cout << "  new collection with " << nsurfs << " panels and " << nnodes << " nodes" << std::endl;

    // pull out the node locations, they go in the base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(nnodes);
      for (size_t i=0; i<nnodes; ++i) {
        this->x[d][i] = _x[Dimensions*i+d];
      }
    }

    // save them as untransformed if we are given a Body pointer
    if (_bp) {
      this->ux = this->x;
    }

    // copy over the node indices (with a possible type change)
    bool idx_are_all_good = true;
    idx.resize(_idx.size());
    const size_t nper = _idx.size() / nsurfs;	// we can have >2 nodes per elem now
    for (size_t i=0; i<nper*nsurfs; ++i) {
      // make sure it exists in the nodes array
      if (_idx[i] >= nnodes) idx_are_all_good = false;
      idx[i] = _idx[i];
    }
    assert(idx_are_all_good && "Some indicies are bad");

    // compute all basis vectors and panel areas
    compute_bases(nsurfs);

    // are strengths/values on a per-node or per-panel basis? - per panel now

    // now, depending on the element type, put the value somewhere
    const size_t nvalper = _val.size() / nsurfs;
    if (this->E == active) {
      // value is a fixed strength for the segment
      // first the tangential/vortex component
      Vector<S> new_s(nsurfs);
      ps[0] = std::move(new_s);
      // and the normal/source component
      Vector<S> new_s1(nsurfs);
      ps[1] = std::move(new_s1);

      // now assign these
      if (nvalper == 1) {
        std::copy(_val.begin(), _val.end(), ps[0]->begin());
        std::fill(ps[1]->begin(), ps[1]->end(), 0.0);
      } else if (nvalper == 2) {
        for (size_t i=0; i<nsurfs; ++i) {
          (*ps[0])[i] = _val[2*i];
          (*ps[1])[i] = _val[2*i+1];
        }
      } else {
        std::fill(ps[0]->begin(), ps[0]->end(), 0.0);
        std::fill(ps[1]->begin(), ps[1]->end(), 0.0);
      }

    } else if (this->E == reactive or this->E == hybrid) {
      // value is a boundary condition
      // bc[0] is the tangential BC, typically 0.0
      Vector<S> new_bc(nsurfs);
      bc[0] = std::move(new_bc);
      // bc[1] is the normal BC, also typically 0.0
      Vector<S> new_bc1(nsurfs);
      bc[1] = std::move(new_bc1);

      // now assign these
      if (nvalper == 1) {
        std::copy(_val.begin(), _val.end(), bc[0]->begin());
        std::fill(bc[1]->begin(), bc[1]->end(), 0.0);
      } else if (nvalper == 2) {
        for (size_t i=0; i<nsurfs; ++i) {
          (*bc[0])[i] = _val[2*i];
          (*bc[1])[i] = _val[2*i+1];
        }
      } else {
        std::fill(bc[0]->begin(), bc[0]->end(), 0.0);
        std::fill(bc[1]->begin(), bc[1]->end(), 0.0);
      }

      // make space for the unknown sheet strengths (only if reactive)
      if (this->E == reactive) {
        for (size_t d=0; d<2; ++d) {
          if (ps[d]) {
            ps[d]->resize(nsurfs);
            std::fill(ps[d]->begin(), ps[d]->end(), 0.0);
          }
        }
      }

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is per node, in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnodes);
    }

    // but panel velocity is per panel
    for (size_t d=0; d<Dimensions; ++d) {
      pu[d].resize(nsurfs);
    }

    // debug print
    if (false) {
      std::cout << "Nodes" << std::endl;
      for (size_t i=0; i<std::min((size_t)100,nnodes); ++i) {
        std::cout << "  " << i << " " << this->x[0][i] << " " << this->x[1][i] << std::endl;
      }
      std::cout << "Segments" << std::endl;
      for (size_t i=0; i<std::min((size_t)100,nsurfs); ++i) {
        std::cout << "  " << i << " " << idx[2*i] << " " << idx[2*i+1] << std::endl;
      }
    }

    // need to reset the base class n and the local np
    this->n = nnodes;
    np = nsurfs;

    // find geometric center
    if (this->M == bodybound) {
      set_geom_center();
    }
  }

  size_t                            get_npanels()     const { return np; }
  const S                           get_vol()         const { return vol; }
  const std::array<S,Dimensions>    get_geom_center() const { return tc; }

  // panel geometry
  const std::vector<Int>&                  get_idx()  const { return idx; }
  const std::array<Vector<S>,Dimensions>&  get_tang() const { return b[0]; }
  const std::array<Vector<S>,Dimensions>&  get_norm() const { return b[1]; }
  const Vector<S>&                         get_area() const { return area; }

  // override the ElementBase versions and send the panel-center vels
  const std::array<Vector<S>,Dimensions>&   get_vel() const { return pu; }
  std::array<Vector<S>,Dimensions>&         get_vel()       { return pu; }

  // fixed or unknown surface strengths, or those due to rotation
  const Vector<S>&                          get_str() const { return *ps[0]; }
  Vector<S>&                                get_str()       { return *ps[0]; }
  const Vector<S>&                     get_vort_str() const { return *ps[0]; }
  Vector<S>&                           get_vort_str()       { return *ps[0]; }
  const bool                           have_src_str() const { return (bool)ps[1]; }
  const Vector<S>&                      get_src_str() const { return *ps[1]; }
  Vector<S>&                            get_src_str()       { return *ps[1]; }

  // and (reactive only) boundary conditions
  const Vector<S>&                     get_tang_bcs() const { return *bc[0]; }
  const Vector<S>&                     get_norm_bcs() const { return *bc[1]; }

  // find out the next row index in the BEM after this collection
  void set_first_row(const Int _i) { istart = _i; }
  const Int num_unknowns_per_panel() const { return (source_str_is_unknown ? 2 : 1); }
  const Int get_first_row() const { return istart; }
  const Int get_num_rows()  const { return (get_npanels()*num_unknowns_per_panel() + (is_augmented() ? 1 : 0)); }
  const Int get_next_row()  const { return istart+get_num_rows(); }

  // assign the new strengths from BEM - do not let base class do this
  void set_str(const size_t ioffset, const size_t icnt, Vector<S> _in) {

    assert(ps[0] && "Strength array does not exist");

    // pop off the "unknown" rotation rate and save it
    if (is_augmented()) {
      solved_omega = _in.back();
      std::cout << "    solved rotation rate is " << solved_omega << std::endl;
      omega_error = solved_omega - this->B->get_rotvel();
      std::cout << "    error in rotation rate is " << omega_error << std::endl;
      _in.pop_back();
    }

    assert(_in.size() == (*ps[0]).size()*num_unknowns_per_panel() && "Set strength array size does not match");
    //assert(ioffset == 0 && "Offset is not zero");

    // copy the BEM-solved strengths into the panel-strength data structures
    if (source_str_is_unknown) {
      // there are source strengths
      for (size_t i=0; i<get_npanels(); ++i) {
        (*ps[0])[i] = _in[2*i];
        (*ps[1])[i] = _in[2*i+1];
        //Int id0 = idx[2*i];
        //Int id1 = idx[2*i+1];
        //std::cout << i << " " << (*ps[0])[i] << " " << this->x[0][id0] << " " << this->x[1][id0] << " " << this->x[0][id1] << " " << this->x[1][id1] << " " << (*ps[1])[i] << std::endl;
      }
    } else {
      // only vortex strengths
      *ps[0] = _in;
      //for (size_t i=0; i<get_npanels(); ++i) {
        //Int id0 = idx[2*i];
        //Int id1 = idx[2*i+1];
        //std::cout << i << " " << (*ps[0])[i] << " " << this->x[0][id0] << " " << this->x[1][id0] << " " << this->x[0][id1] << " " << this->x[1][id1] << std::endl;
      //}
    }
  }

  // a little logic to see if we should augment the BEM equations for this object
  const bool is_augmented() const {
    bool augment = true;

    // don't augment the ground body, or the boundary to an internal flow
    if (this->B) {
      // is the body pointer ground?
      if (std::string("ground").compare(this->B->get_name()) == 0) {
        // now, does the object bound internal flow?
        if (vol < 0.0) augment = false;
      }
    } else {
      // nullptr for Body? no augment (old way of turning it off)
      augment = false;
    }
    // and only need to augment reactive surfaces (ones participating in BEM)
    if (this->E != reactive) augment = false;

    // force no augmentation at all
    //augment = false;

    return augment;
  }

  const float get_max_bc_value() const {
    if (this->E == reactive) {
      const S vort_max = *std::max_element(std::begin(*bc[0]), std::end(*bc[0]));
      const S vort_min = *std::min_element(std::begin(*bc[0]), std::end(*bc[0]));
      S this_max = std::max(vort_max, -vort_min);
      const S src_max = *std::max_element(std::begin(*bc[1]), std::end(*bc[1]));
      const S src_min = *std::min_element(std::begin(*bc[1]), std::end(*bc[1]));
      this_max = std::max(this_max, std::max(src_max, -src_min));
      return (float)this_max;
    } else {
      return (float)0.0;
    }
  }

  // append nodes and panels to this collection
  void add_new(const ElementPacket<float>& _in) {

    // ensure that this packet really is Surfaces
    assert(_in.idx.size() != 0 && "Input ElementPacket is not Surfaces");
    assert(_in.ndim == 1 && "Input ElementPacket is not Surfaces");

    // if no surfs, quit out now
    const size_t nsurfs = _in.nelem;
    if (nsurfs == 0) return;

    // and that it has the right number of values per elem
    if (this->E == inert) assert(_in.val.size() == 0 && "Input ElementPacket with inert Surfaces has nonzero val array");

    // make sure input arrays are correctly-sized
    assert(_in.x.size() % Dimensions == 0 && "Position array is not an even multiple of dimensions");
    const size_t nnodes = _in.x.size() / Dimensions;
    assert(_in.idx.size() % nsurfs == 0 && "Index array is not an even multiple of panel count");
    const size_t nper = _in.idx.size() / nsurfs;	// we can have >2 nodes per elem now
    assert(_in.idx.size() % nper == 0 && "Index array is not an even multiple of nper");
    assert(_in.val.size() % nsurfs == 0 && "Value array is not an even multiple of panel count");

    std::cout << "  adding " << nsurfs << " new surface panels and " << nnodes << " new points to collection..." << std::endl;

    // remember old sizes of nodes and element arrays
    const size_t nnold = this->n;
    const size_t neold = get_npanels();
    const size_t nperold = idx.size() / neold;
    assert(nper == nperold && "Cannot change element order in a single Collection");

    // DON'T call the method in the base class, because we do things differently here
    //ElementBase<S>::add_new(_in);

    // pull out the node locations, they are base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(nnold+nnodes);
      for (size_t i=0; i<nnodes; ++i) {
        this->x[d][nnold+i] = _in.x[Dimensions*i+d];
      }
    }

    // save them as untransformed if we have a Body pointer
    if (this->B) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*this->ux)[d].resize(nnold+nnodes);
        for (size_t i=nnold; i<nnold+nnodes; ++i) {
          (*this->ux)[d][i] = this->x[d][i];
        }
      }
    }

    // copy over the node indices, taking care to offset into the new array
    bool idx_are_all_good = true;
    idx.resize(nper*neold + _in.idx.size());
    for (size_t i=0; i<nper*nsurfs; ++i) {
      // make sure it exists in the nodes array
      if (_in.idx[i] >= nnold+nnodes) idx_are_all_good = false;
      idx[nper*neold+i] = nnold + _in.idx[i];
    }
    assert(idx_are_all_good && "Some indicies are bad");

    // compute all basis vectors and panel areas
    compute_bases(neold+nsurfs);

    // now, depending on the element type, put the value somewhere - but panel-wise, so here
    const size_t nvalper = _in.val.size() / nsurfs;
    if (this->E == active) {
      // value is a fixed strength for the element
      ps[0]->resize(neold+nsurfs); 
      ps[1]->resize(neold+nsurfs); 

      // now assign these
      if (nvalper == 1) {
        std::copy(_in.val.begin(), _in.val.end(), ps[0]->begin()+neold);
        std::fill(ps[1]->begin()+neold, ps[1]->end(), 0.0);
      } else if (nvalper == 2) {
        for (size_t i=0; i<nsurfs; ++i) {
          const size_t iat = neold + i;
          (*ps[0])[iat] = _in.val[2*i];
          (*ps[1])[iat] = _in.val[2*i+1];
        }
      } else {
        std::fill(ps[0]->begin()+neold, ps[0]->end(), 0.0);
        std::fill(ps[1]->begin()+neold, ps[1]->end(), 0.0);
      }

    } else if (this->E == reactive) {
      // value is a boundary condition
      bc[0]->resize(neold+nsurfs); 
      bc[1]->resize(neold+nsurfs); 

      // now assign these
      if (nvalper == 1) {
        std::copy(_in.val.begin(), _in.val.end(), bc[0]->begin()+neold);
        std::fill(bc[1]->begin()+neold, bc[1]->end(), 0.0);
      } else if (nvalper == 2) {
        for (size_t i=0; i<nsurfs; ++i) {
          const size_t iat = neold + i;
          (*bc[0])[iat] = _in.val[2*i];
          (*bc[1])[iat] = _in.val[2*i+1];
        }
      } else {
        std::fill(bc[0]->begin()+neold, bc[0]->end(), 0.0);
        std::fill(bc[1]->begin()+neold, bc[1]->end(), 0.0);
      }

      // upsize vortex sheet and raw strength arrays, too
      for (size_t d=0; d<2; ++d) {
        if (ps[d]) {
          ps[d]->resize(neold+nsurfs);
          std::fill(ps[d]->begin()+neold, ps[d]->end(), 0.0);
        }
      }

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnold+nnodes);
    }

    // panel velocity is here
    for (size_t d=0; d<Dimensions; ++d) {
      pu[d].resize(neold+nsurfs);
    }

    // debug print
    if (false) {
      std::cout << "Nodes" << std::endl;
      for (size_t i=0; i<nnold+nnodes; ++i) {
        std::cout << "  " << i << " " << this->x[0][i] << " " << this->x[1][i] << std::endl;
      }
      std::cout << "Segments" << std::endl;
      for (size_t i=0; i<neold+nsurfs; ++i) {
        std::cout << "  " << i << " " << idx[2*i] << " " << idx[2*i+1] << std::endl;
      }
    }

    // need to reset the base class n
    this->n += nnodes;
    np += nsurfs;

    // re-find geometric center
    if (this->M == bodybound) {
      set_geom_center();
    }
  }


  void add_body_motion(const S _factor, const double _time) {
    // no need to call base class now
    //ElementBase<S>::add_body_motion(_factor);

    // if no body pointer, or attached to ground, return
    if (not this->B) return;
    if (std::string("ground").compare(this->B->get_name()) == 0) return;

    // make sure we've calculated transformed center (we do this when we do volume)
    assert(vol > 0.0 && "Have not calculated transformed center, or volume is negative");
    // and we trust that we've transformed utc to tc

    // do this for all nodes - what about panels?
    for (size_t i=0; i<get_npanels(); ++i) {

      // apply the translational velocity
      const std::array<double,Dimensions> thisvel = this->B->get_vel(_time);
      for (size_t d=0; d<Dimensions; ++d) {
        pu[d][i] += _factor * (float)thisvel[d];
      }

      // now compute the rotational velocity with respect to the geometric center
      const double thisrotvel = this->B->get_rotvel(_time);
      // center of this panel
      Int id0 = idx[2*i];
      Int id1 = idx[2*i+1];
      // panel center
      const S xc = 0.5 * (this->x[0][id1] + this->x[0][id0]);
      const S yc = 0.5 * (this->x[1][id1] + this->x[1][id0]);
      // add rotational velocity
      pu[0][i] -= _factor * (float)thisrotvel * (yc - tc[1]);
      pu[1][i] += _factor * (float)thisrotvel * (xc - tc[0]);
    }
  }
 
  void zero_strengths() {
    // call base class first
    ElementBase<S>::zero_strengths();

    // and reset any panel vortex or source strengths
    for (size_t d=0; d<2; ++d) {
      if (ps[d]) {
        std::fill(ps[d]->begin(), ps[d]->end(), 0.0);
      }
    }
  }

  // three ways to add source and vortex rotational strengths to the surface
  // first: as a multiple of the current, defined rotation rate
  void add_rot_strengths(const S _factor) {
    // if no parent Body, forget it
    if (not this->B) return;
    // get current rotation rate
    const S rotvel = (S)this->B->get_rotvel();
    // call parent
    add_rot_strengths_base(_factor * rotvel);
  }

  // second: assuming unit rotation rate (for finding the BEM influence matrix)
  void add_unit_rot_strengths() {
    add_rot_strengths_base(1.0);
  }

  // third: as a multiple of the rotation rate
  void add_solved_rot_strengths(const S _factor) {
    if (is_augmented()) {
      // use the augmented-BEM result for rotation rate
      add_rot_strengths_base(_factor * solved_omega);
    } else {
      // use the predefined rotationrate
      add_rot_strengths(_factor);
    }
  }

  // augment the strengths with a value equal to that which accounts for
  //   the solid-body rotation of the object
  // NOTE: this needs to provide both the vortex AND source strengths!
  // AND we don't have the time - assume bodies have been transformed
  void add_rot_strengths_base(const S _rotvel) {

    // No - ALWAYS allow this to run
    // if no rotation, strengths, or no parent Body, or attached to ground, then no problem!
    //if (not this->B) return;
    //if (not ps[0]) return;
    //if (std::string("ground").compare(this->B->get_name()) == 0) return;

    //if (std::abs(_rotvel) < std::numeric_limits<float>::epsilon()) return;

    // make sure we've calculated transformed center (we do this when we do volume)
    // NO - we can't check this because internal flow setups will have negative volume
    //assert(vol > 0.0 && "Have not calculated transformed center, or volume is negative");
    // and we trust that we've transformed utc to tc

    assert(ps[0]->size() == get_npanels() && "Strength array is not the same as panel count");
    assert(ps[1]->size() == get_npanels() && "Strength array is not the same as panel count");

    // still here? let's do it. use the untransformed coordinates
    for (size_t i=0; i<get_npanels(); i++) {
      const size_t j   = idx[2*i];
      const size_t jp1 = idx[2*i+1];
      // vector from object geometric center to panel center
      const S dx = 0.5 * ((*this->ux)[0][j] + (*this->ux)[0][jp1]) - utc[0];
      const S dy = 0.5 * ((*this->ux)[1][j] + (*this->ux)[1][jp1]) - utc[1];
      // velocity of the panel center
      const S ui = -_rotvel * dy;
      const S vi =  _rotvel * dx;

      // panel tangential vector, fluid to the left, body to the right
      S panelx = (*this->ux)[0][jp1] - (*this->ux)[0][j];
      S panely = (*this->ux)[1][jp1] - (*this->ux)[1][j];
      const S oopanell = 1.0 / area[i];
      panelx *= oopanell;
      panely *= oopanell;

      // the vortex strength - we ADD to the existing
      const S new_vort = -1.0 * (ui*panelx + vi*panely);
      (*ps[0])[i] += new_vort;

      // the source strength
      const S new_src = -1.0 * (ui*panely - vi*panelx);
      (*ps[1])[i] += new_src;

      // debug print
      if (std::abs(_rotvel) > 0.0 and false) {
        std::cout << "  panel " << i << " at " << dx << " " << dy << " adds to vortex str "
                  << new_vort << " and source str " << new_src << std::endl;
      }
    }
  }

  // calculate the geometric center of all geometry in this object
  void set_geom_center() {

    // we must have an attached body and a set of untransformed coordinates
    assert(this->B && "Body pointer has not been set");
    assert(this->ux && "Untransformed positions have not been set");

    std::cout << "  inside Surfaces::set_geom_center with " << get_npanels() << " panels" << std::endl;

    // iterate over panels, accumulating area and CM in double precision
    using A = double;
    A asum = 0.0;
    A xsum = 0.0;
    A ysum = 0.0;
    for (size_t i=0; i<get_npanels(); i++) {
      const size_t j   = idx[2*i];
      const size_t jp1 = idx[2*i+1];
      //std::cout << "elem " << i << " from " << (*this->ux)[0][j] << " " << (*this->ux)[1][j] << " to " << (*this->ux)[0][jp1] << " " << (*this->ux)[1][jp1] << std::endl;

      // assume a triangle from 0,0 to two ends of each panel
      const A xc = (0.0 + (*this->ux)[0][j] + (*this->ux)[0][jp1]) / 3.0;
      const A yc = (0.0 + (*this->ux)[1][j] + (*this->ux)[1][jp1]) / 3.0;
      const A panelx = (*this->ux)[0][jp1] - (*this->ux)[0][j];
      const A panely = (*this->ux)[1][jp1] - (*this->ux)[1][j];

      // and the side lengths
      const A a = std::sqrt(std::pow((*this->ux)[0][j],2)+std::pow((*this->ux)[1][j],2));
      const A b = std::sqrt(std::pow(panelx,2)+std::pow(panely,2));
      const A c = std::sqrt(std::pow((*this->ux)[0][jp1],2)+std::pow((*this->ux)[1][jp1],2));
      //std::cout << "  panel " << i << " has side lens " << a << " " << b << " " << c << std::endl;

      // Heron's formula for the area
      const A hs = 0.5*(a+b+c);
      //std::cout << "  hs " << hs << " comps " << a << " " << b << " " << c << std::endl;
      const A asqrd = hs*(hs-a)*(hs-b)*(hs-c);
      A thisarea = 0.0;
      // asqrd can be minutely less than zero due to roundoff errors, that throws nans
      if (asqrd > 0.0) {
        thisarea = std::sqrt(asqrd);
      }
      // negate area if the winding is backwards
      if ((*this->ux)[1][j]*panelx - (*this->ux)[0][j]*panely < 0.0) thisarea = -thisarea;

      // add this to the running sums
      //std::cout << "    and area " << thisarea << " and center " << xc << " " << yc << std::endl;
      asum += thisarea;
      xsum += xc*thisarea;
      ysum += yc*thisarea;
      //std::cout << "    running sums " << asum << " " << xsum << " " << ysum << std::endl;
    }
    vol = (S)asum;
    utc[0] = (S)(xsum/asum);
    utc[1] = (S)(ysum/asum);

    std::cout << "    geom center is " << utc[0] << " " << utc[1] << " and area is " << vol << std::endl;
  }

  // need to maintain the 2x2 set of basis vectors for each panel
  // this also calculates the triangle areas
  // always recalculate everything!
  void compute_bases(const Int nnew) {

    // since we may be saving higher-order elements, this is no longer valid
    //assert(2*nnew == idx.size() && "Array size mismatch");

    // resize any vectors
    for (size_t i=0; i<Dimensions; ++i) {
      for (size_t j=0; j<Dimensions; ++j) {
        b[i][j].resize(nnew);
      }
    }
    area.resize(nnew);

    // we'll reuse these vectors
    std::array<S,Dimensions> x1, norm;

    // update everything
    for (size_t i=0; i<nnew; ++i) {
      const size_t id0 = idx[2*i];
      const size_t id1 = idx[2*i+1];
      //std::cout << "elem near " << this->x[0][id0] << " " << this->x[1][id0] << std::endl;

      // x1 vector is along direction from node 0 to node 1
      for (size_t j=0; j<Dimensions; ++j) x1[j] = this->x[j][id1] - this->x[j][id0];
      const S base = std::sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
      for (size_t j=0; j<Dimensions; ++j) x1[j] *= (1.0/base);
      //std::cout << "  has x1 " << x1[0] << " " << x1[1] << std::endl;

      // now we have the area
      area[i] = base;
      //std::cout << "elem near " << this->x[0][id0] << " " << this->x[1][id0] << " has area " << area[i] << std::endl;

      // the normal vector points into the fluid (to the left when walking from nodes 0 to 1)
      norm[0] = -x1[1];
      norm[1] =  x1[0];
      //std::cout << "  norm " << norm[0] << " " << norm[1] << std::endl;

      // and assign
      for (size_t j=0; j<Dimensions; ++j) b[0][j][i] = x1[j];
      for (size_t j=0; j<Dimensions; ++j) b[1][j][i] = norm[j];

      //std::cout << "elem near " << this->x[0][id0] << " " << this->x[1][id0] << " has norm " << b[1][0][i] << " " << b[1][1][i] << " and tang " << b[0][0][i] << " " << b[0][1][i] << std::endl;
    }
  }

  // when transforming a body-bound object to a new time, we must also transform the geometric center
  void transform(const double _time) {
    // must explicitly call the method in the base class
    ElementBase<S>::transform(_time);

    // and recalculate the basis vectors
    compute_bases(np);

    if (this->B and this->M == bodybound) {
    //if (this->B) {
      // prepare for the transform
      std::array<double,Dimensions> thispos = this->B->get_pos();
      const double theta = this->B->get_orient();
      const S st = std::sin(theta);
      const S ct = std::cos(theta);

      // transform the utc to tc here
      tc[0] = (S)thispos[0] + utc[0]*ct - utc[1]*st;
      tc[1] = (S)thispos[1] + utc[0]*st + utc[1]*ct;

    } else {
      // transform the utc to tc here
      tc[0] = utc[0];
      tc[1] = utc[1];
    }
  }


  void zero_vels() {
    // zero the local, panel-center vels
    for (size_t d=0; d<Dimensions; ++d) {
      std::fill(pu[d].begin(), pu[d].end(), 0.0);
    }
    // then explicitly call the method in the base class to zero theirs
    ElementBase<S>::zero_vels();
  }

  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    // finalize panel-center vels first
    const double factor = 0.5/M_PI;
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<get_npanels(); ++i) {
        pu[d][i] = _fs[d] + pu[d][i] * factor;
      }
    }
    // must explicitly call the method in the base class, too
    ElementBase<S>::finalize_vels(_fs);
  }


  //
  // 1st order Euler advection and stretch
  //
  void move(const double _time, const double _dt,
            const double _wt1, Surfaces<S> const & _u1) {

    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt, _wt1, _u1);

    // and update the max strength measure
    (void) update_max_str();
  }


  //
  // 2nd order RK advection and stretch
  //
  void move(const double _time, const double _dt,
            const double _wt1, Surfaces<S> const & _u1,
            const double _wt2, Surfaces<S> const & _u2) {

    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt, _wt1, _u1, _wt2, _u2);

    // must confirm that incoming time derivates include velocity (?)

    // and update the max strength measure
    (void) update_max_str();
  }


  //
  // 3rd order RK advection and stretch
  //
  void move(const double _time, const double _dt,
            const double _wt0, Surfaces<S> const & _u0,
            const double _wt1, Surfaces<S> const & _u1,
            const double _wt2, Surfaces<S> const & _u2) {

    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt, _wt0, _u0, _wt1, _u1, _wt2, _u2);

    // must confirm that incoming time derivates include velocity (?)

    // and update the max strength measure
    (void) update_max_str();
  }


  //
  // return a particle version of the panels (useful during Diffusion)
  // offset is in world units - NOT scaled
  //
  ElementPacket<float> represent_as_particles(const S _offset, const S _ips) {

    // prepare the vectors
    std::vector<float> _x;
    std::vector<Int> _idx;
    std::vector<float> _vals;

    // since there may be more points than panels, reserve space and emplace_back
    _x.reserve(Dimensions*get_npanels());
    _vals.reserve(get_npanels());

    // get basis vectors
    std::array<Vector<S>,2>& norm = b[1];

    const size_t nper = idx.size() / get_npanels();	// we can have >2 nodes per elem now

    // if caller does not know ips, it sends -1.0, which we interpret here
    S thisips = _ips;
    if (_ips < 0.0) thisips = area[0];

    // the fluid is to the left walking from one point to the next
    // so go CW around an external boundary starting at theta=0 (+x axis)

    for (size_t i=0; i<get_npanels(); ++i) {
      // panel end points
      const Int id0 = idx[nper*i];
      const Int id1 = idx[nper*i+1];
      // the panel strength is the solved strength plus the boundary condition
      float this_str = (*ps[0])[i];
      // add on the (vortex) bc value here
      if (this->E == reactive) this_str += (*bc[0])[i];

      // is the panel long enough to shed more than one particle?
      const size_t nshed = 1 + std::max(0, (int)(area[i] / (1.2*thisips)));

      for (size_t j=0; j<nshed; ++j) {
        const S wgt = (1+2*j) / (S)(2*nshed);
        // start on the panel and push out a fixed distance
        _x.emplace_back(wgt*this->x[0][id1] + (1.0-wgt)*this->x[0][id0] + _offset * norm[0][i]);
        _x.emplace_back(wgt*this->x[1][id1] + (1.0-wgt)*this->x[1][id0] + _offset * norm[1][i]);
        // complete the element with a strength
        _vals.emplace_back(this_str * area[i] / (S)nshed);
      }
    }

    return ElementPacket<float>({_x, _idx, _vals, _vals.size(), 0});
  }

  // find the new peak vortex strength magnitude
  S get_max_str() {
    if (ps[0]) {
      const S this_max = *std::max_element(std::begin(*ps[0]), std::end(*ps[0]));
      const S this_min = *std::min_element(std::begin(*ps[0]), std::end(*ps[0]));
      return std::max(this_max, -this_min);
    } else {
      return 1.0;
    }
  }

  // smooth the peak strength magnitude
  void update_max_str() {
    S thismax = get_max_str();

    // and slowly update the current saved value
    if (max_strength < 0.0) {
      max_strength = thismax;
    } else {
      max_strength = 0.1*thismax + 0.9*max_strength;
    }
  }

  // add and return the total circulation of all elements
  //   specializing the one in ElementBase because we need
  //   to scale by panel length here
  S get_total_circ(const double _time) {
    S circ = 0.0;

    if (ps[0]) {
      // we have strengths, add them up
      for (size_t i=0; i<get_npanels(); i++) {
        // total circulation is just vortex sheet strength times panel length
        circ += (*ps[0])[i] * area[i];
      }
    }

    return circ;
  }

  // add and return the total circulation of all elements
  S get_body_circ(const double _time) {
    S circ = 0.0;

    // do not call the parent
    if (this->B) {
      // we're attached to a body - great! what's the rotation rate?
      circ = 2.0 * vol * (S)this->B->get_rotvel(_time);
    } else {
      // we are fixed, thus not rotating
    }

    return circ;
  }

  // return the rotational speed from the previous step
  S get_last_body_circ() { return 2.0 * vol * (S)this_omega; }

  // compute the change in circulation of the rotating body since the last shedding event
/*
  S get_body_circ_change(const double _time) {
    const S old_circ = 2.0 * vol * (S)this_omega;
    S new_circ = 0.0;

    // do not call the parent
    if (this->B) {
      // we're attached to a body - great! what's the rotation rate?
      new_circ = 2.0 * vol * (S)this->B->get_rotvel(_time);
    } else {
      // we are fixed, thus not rotating
    }

    return (new_circ - old_circ);
  }
*/

  // add and return the total impulse of all elements
  std::array<S,Dimensions> get_total_impulse() {

    // here is the return vector
    std::array<S,Dimensions> imp;
    imp.fill(0.0);

    if (ps[0]) {
      // make this easy - represent as particles - do we count BCs here?!?
      ElementPacket<float> pts = represent_as_particles(0.0, -1.0);

      // now compute impulse of those
      for (size_t i=0; i<pts.nelem; ++i) {
        imp[0] -= (S)pts.val[i] * (S)pts.x[2*i];
        imp[1] += (S)pts.val[i] * (S)pts.x[2*i+1];
      }
    }

    return imp;
  }

  // total the arc length for this collection
  S get_total_arclength() {
    S outsum = S(0.0);
    for (size_t i=0; i<get_npanels(); ++i) {
      outsum += area[i];
    }
    return outsum;
  }

  // total the inlet volume flow rate for this collection
  S get_total_inflow() {
    if (this->E != reactive and this->E != hybrid) return S(0.0);

    S insum = S(0.0);
    for (size_t i=0; i<get_npanels(); ++i) {
      if ((*bc[1])[i] > S(0.0)) {
        // length times flow rate
        insum += (*bc[1])[i] * area[i];
      }
    }

    return insum;
  }

  // if no outlets were defined, must zero the inflow to allow sim to proceed
  void zero_inflow() {
    if (this->E != reactive and this->E != hybrid) return;
    for (size_t i=0; i<get_npanels(); ++i) {
      (*bc[1])[i] = 0.0;
    }
  }

  // total the inlet volume flow rate for this collection
  S get_total_outflow() {
    if (this->E != reactive and this->E != hybrid) return S(0.0);

    S outsum = S(0.0);
    for (size_t i=0; i<get_npanels(); ++i) {
      if ((*bc[1])[i] < S(0.0)) {
        // length times flow rate
        outsum += -(*bc[1])[i] * area[i];
      }
    }

    return outsum;
  }

  // linearly scale the outflow (normal BC vels) to enforce continuity
  void scale_outflow(const S factor) {
    if (this->E != reactive and this->E != hybrid) return;
    for (size_t i=0; i<get_npanels(); ++i) {
      if ((*bc[1])[i] < S(0.0)) {
        (*bc[1])[i] *= factor;
      }
    }
  }

  // outright set the outflow (normal BC vels) (ideally to enforce continuity)
  void set_outflow(const S factor) {
    if (this->E != reactive and this->E != hybrid) return;
    for (size_t i=0; i<get_npanels(); ++i) {
      (*bc[1])[i] = -factor;
    }
  }

  // reset the circulation counter and saved rotation rate - useful for augmented BEM
  void reset_augmentation_vars() {
    this_omega = this->B->get_rotvel();
    reabsorbed_gamma = 0.0;
  }

  S get_last_body_circ_error() {
    return (S)(2.0 * vol * omega_error);
  }

  // *add* the given circulation to the reabsorbed accumulator
  void add_to_reabsorbed(const S _circ) {
    reabsorbed_gamma += _circ;
  }

  // return that amount of reabsorbed circulation
  S get_reabsorbed() {
    return reabsorbed_gamma;
  }


#ifdef USE_GL
  //
  // OpenGL functions
  //

  // helper function to clean up initGL
  void prepare_opengl_buffer(GLuint _prog, GLuint _idx, const GLchar* _name) {
    glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[_idx]);
    const GLint position_attribute = glGetAttribLocation(_prog, _name);
    // Specify how the data for position can be accessed
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    // Enable the attribute
    glEnableVertexAttribArray(position_attribute);
  }

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor,
              float*              _defcolor) {

    //std::cout << "inside Surfaces.initGL" << std::endl;
    std::cout << "inside Surfaces.initGL with E=" << this->E << " and M=" << this->M << std::endl;

    // generate the opengl state object with space for 4 vbos and 1 shader program
    mgl = std::make_shared<GlState>(4,1);

    // Allocate space, but don't upload the data from CPU to GPU yet
    for (size_t i=0; i<Dimensions; ++i) {
      glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
      glBufferData(GL_ARRAY_BUFFER, 0, this->x[i].data(), GL_STATIC_DRAW);
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mgl->vbo[Dimensions]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, idx.data(), GL_STATIC_DRAW);

    if (this->s) {
      glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[3]);
      glBufferData(GL_ARRAY_BUFFER, 0, (*ps[0]).data(), GL_STATIC_DRAW);
    }

    // Load and create the blob-drawing shader program
    mgl->spo[0] = create_draw_surface_line_prog();

    // Now do the four arrays
    prepare_opengl_buffer(mgl->spo[0], 0, "px");
    prepare_opengl_buffer(mgl->spo[0], 1, "py");
    prepare_opengl_buffer(mgl->spo[0], 2, "rawstr");

    // and for the compute shaders! (later)

    // Get the location of the attributes that enters in the vertex shader
    mgl->projmat_attribute = glGetUniformLocation(mgl->spo[0], "Projection");

    // upload the projection matrix
    glUniformMatrix4fv(mgl->projmat_attribute, 1, GL_FALSE, _projmat.data());

    // locate where the colors and color scales go
    mgl->pos_color_attribute = glGetUniformLocation(mgl->spo[0], "pos_color");
    mgl->neg_color_attribute = glGetUniformLocation(mgl->spo[0], "neg_color");
    mgl->def_color_attribute = glGetUniformLocation(mgl->spo[0], "def_color");
    mgl->str_scale_attribute = glGetUniformLocation(mgl->spo[0], "str_scale");

    // send the current values
    glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_poscolor);
    glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_negcolor);
    glUniform4fv(mgl->def_color_attribute, 1, (const GLfloat *)_defcolor);
    glUniform1f (mgl->str_scale_attribute, (const GLfloat)1.0);
    //std::cout << "init pos color as " << _poscolor[0] << " " << _poscolor[1] << " " << _poscolor[2] << " " << _poscolor[3] << std::endl;

    // and indicate the fragment color output
    glBindFragDataLocation(mgl->spo[0], 0, "frag_color");

    glBindVertexArray(0);
  }

  // this gets done every time we change the size of the index array
  void updateGL() {
    //std::cout << "inside Surfaces.updateGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) return;
    if (glIsVertexArray(mgl->vao) == GL_FALSE) return;

    const size_t vlen = this->x[0].size()*sizeof(S);
    if (vlen > 0) {
      glBindVertexArray(mgl->vao);

      // Indicate and upload the data from CPU to GPU
      for (size_t i=0; i<Dimensions; ++i) {
        // the positions
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
        glBufferData(GL_ARRAY_BUFFER, vlen, this->x[i].data(), GL_DYNAMIC_DRAW);
      }

      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mgl->vbo[Dimensions]);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Int)*idx.size(), idx.data(), GL_DYNAMIC_DRAW);

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {
        // just don't upload strengths

      } else { // this->E is active or reactive
        // the strengths
        if (ps[0]) {
          const size_t slen = ps[0]->size()*sizeof(S);
          glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[3]);
          glBufferData(GL_ARRAY_BUFFER, slen, (*ps[0]).data(), GL_DYNAMIC_DRAW);
        }
      }

      glBindVertexArray(0);

      // must tell draw call how many elements are there - or, really, how many indices
      mgl->num_uploaded = idx.size();
    }
  }

  // OpenGL3 stuff to draw segments, called once per frame
  void drawGL(std::vector<float>& _projmat,
              RenderParams&       _rparams,
              const float         _vdelta) {

    //std::cout << "inside Surfaces.drawGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) {
      initGL(_projmat, _rparams.pos_circ_color,
                       _rparams.neg_circ_color,
                       _rparams.default_color);
      updateGL();
    }

    if (mgl->num_uploaded > 0) {
      glBindVertexArray(mgl->vao);

      // get blending ready
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      glLineWidth(2.0);

      // here is where we split on element type: active/reactive vs. inert
      //if (this->E == inert) {
      //} else { // this->E is active or reactive

      // draw as lines
      glUseProgram(mgl->spo[0]);

      // upload the current projection matrix
      glUniformMatrix4fv(mgl->projmat_attribute, 1, GL_FALSE, _projmat.data());

      // upload the current color values
      glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_rparams.pos_circ_color);
      glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_rparams.neg_circ_color);
      glUniform4fv(mgl->def_color_attribute, 1, (const GLfloat *)_rparams.default_color);
      glUniform1f (mgl->str_scale_attribute, (const GLfloat)max_strength);

      // the one draw call here
      glDrawElements(GL_LINES, mgl->num_uploaded, get_gl_type<Int>, 0);

      // return state
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glBindVertexArray(0);
    }
  }
#endif

  std::string to_string() const {
    std::string retstr = " " + std::to_string(get_npanels()) + ElementBase<S>::to_string() + " Panels";
    return retstr;
  }

  std::string write_vtk(const size_t _index, const size_t _frameno, const double _time) {
    assert(this->np > 0 && "Inside write_vtu_panels with no panels");
  
    const bool asbase64 = true;
    bool has_vort_str = false;
    bool has_src_str = false;
    std::string prefix = "panel_";
    if (this->E != inert) {
      has_vort_str = true;
      has_src_str = (bool)this->ps[1];
    }
  
    // generate file name
    std::stringstream vtkfn;
    vtkfn << prefix << std::setfill('0') << std::setw(2) << _index << "_" << std::setw(5) << _frameno << ".vtu";
    VtkXmlWriter panelWriter = VtkXmlWriter(vtkfn.str(), asbase64);
    // push comment with sim time?
 
    // include simulation time here
    panelWriter.addElement("FieldData");
    {
      std::map<std::string, std::string> attribs = {{"type",           "Float64"},
                                                    {"Name",           "TimeValue"},
                                                    {"NumberOfTuples", "1"}};
      panelWriter.addElement("DataArray", attribs);
      Vector<double> time_vec = {_time};
      panelWriter.writeDataArray(time_vec);
      panelWriter.closeElement();
    }
    // FieldData
    panelWriter.closeElement();
 
    {
      std::map<std::string, std::string> attribs = {{"NumberOfPoints", std::to_string(this->n).c_str()},
                                                    {"NumberOfCells", std::to_string(this->np).c_str()}};
      panelWriter.addElement("Piece", attribs);
    }

    panelWriter.addElement("Points");
    {
      std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                    {"Name",               "position"},
                                                    {"type",               "Float32"}}; 
      panelWriter.addElement("DataArray", attribs);
      Vector<float> pos = panelWriter.unpackArray(this->x);
      panelWriter.writeDataArray(pos);
      panelWriter.closeElement();
    }
    // Points
    panelWriter.closeElement();
 
    panelWriter.addElement("Cells");
    // again, all connectivities and offsets must be Int32!
    {
      std::map<std::string, std::string> attribs = {{"Name", "connectivity"},
                                                    {"type", "Int32"}};
      panelWriter.addElement("DataArray", attribs);
      std::vector<Int> const & idx = this->idx;
      Vector<int32_t> v(std::begin(idx), std::end(idx));
      panelWriter.writeDataArray(v);
      panelWriter.closeElement();
    }
 
    {
      std::map<std::string, std::string> attribs = {{"Name", "offsets"},
                                                    {"type", "Int32"}};
      panelWriter.addElement("DataArray", attribs);
      Vector<int32_t> v(this->np);
      std::iota(v.begin(), v.end(), 1);
      std::transform(v.begin(), v.end(), v.begin(),
                     std::bind(std::multiplies<int32_t>(), std::placeholders::_1, 2));
      panelWriter.writeDataArray(v);
      panelWriter.closeElement();
    }
 
    {
      std::map<std::string, std::string> attribs = {{"Name", "types"},
                                                    {"type", "UInt8"}};
      panelWriter.addElement("DataArray", attribs);
      Vector<uint8_t> v(this->np);
      std::fill(v.begin(), v.end(), 3);
      panelWriter.writeDataArray(v);
      panelWriter.closeElement();
    }
    // Cells
    panelWriter.closeElement();
 
    // do we have CellData?
    std::string scalar_list;
    std::string vector_list;
    if (has_vort_str) scalar_list.append("vortex sheet strength,");
    if (has_src_str) scalar_list.append("source sheet strength,");
    //if (has_radii) scalar_list.append("area,");

    if (scalar_list.size() + vector_list.size() > 0) {

      {
        std::map<std::string, std::string> attribs;
        if (scalar_list.size()>1) {
          scalar_list.pop_back();
          attribs.insert({"Scalars", scalar_list});
        }
        panelWriter.addElement("CellData", attribs);
      }
 
      if (has_vort_str) {
        std::map<std::string, std::string> attribs = {{"Name", "vortex sheet strength"},
                                                      {"type", "Float32"}};
        panelWriter.addElement("DataArray", attribs);
        panelWriter.writeDataArray(*this->ps[0]);
        panelWriter.closeElement();
      }
 
      if (has_src_str) {
        std::map<std::string, std::string> attribs = {{"Name", "source sheet strength"},
                                                      {"type", "Float32"}};
        panelWriter.addElement("DataArray", attribs);
        panelWriter.writeDataArray(*this->ps[1]);
        panelWriter.closeElement();
      }

      // CellData
      panelWriter.closeElement();
    }

    // do we have PointData?
    bool havePointData = false;

    if (havePointData) {

      {
        std::map<std::string, std::string> attribs = {{"Vectors", "velocity"}};
        panelWriter.addElement("CellData", attribs);
      }
 
      // should write vels on a node-wise basis if type is "active"
      {
        std::map<std::string, std::string> attribs = {{"NumberOfComponents", "3"},
                                                      {"Name",               "velocity"},
                                                      {"type",               "Float32"}}; 
        panelWriter.addElement("DataArray", attribs);
        Vector<float> vel = panelWriter.unpackArray(this->u);
        panelWriter.writeDataArray(vel);
        panelWriter.closeElement();
      }

      // PointData
      panelWriter.closeElement();
    }

    // Piece 
    panelWriter.closeElement();
  
    panelWriter.finish();
    std::cout << "Wrote " << this->np << " panels to " << vtkfn.str() << std::endl;
    return vtkfn.str();
  }

protected:
  // ElementBase.h has x, s, u, ux on the *nodes*

  size_t np;				// number of panels

  // element-wise variables special to triangular panels
  std::vector<Int>                 idx;	// indexes into the x array
  Vector<S>                       area; // panel areas
  Basis<S>                           b; // transformed basis vecs: tangent is b[0], normal is b[1], x norm is b[1][0]
  std::array<Vector<S>,Dimensions>  pu; // velocities on panel centers - "u" is node vels in ElementBase

  // strengths and BCs
  Strength<S>                       ps; // panel-wise strengths per unit length (for "active" and "reactive")
                                        // vortex sheet strengths are 0, source is in 1
  Strength<S>                       bc; // boundary condition for the elements (only when "reactive")
  bool           source_str_is_unknown; // should the BEM solve for source strengths?

  // parameters for the encompassing body
  Int                           istart; // index of first entry in RHS vector and A matrix
  S                                vol; // volume of the body - for augmented BEM solution
  std::array<S,Dimensions>         utc; // untransformed geometric center
  std::array<S,Dimensions>          tc; // transformed geometric center

  // augmented-BEM-related
  double                  solved_omega; // rotation rate returned from augmented row in BEM
  double                   omega_error; // error in rotation rate from last BEM
  double                    this_omega; // rotation rate at most recent Diffusion step
  S                   reabsorbed_gamma; // amount of circulation reabsorbed by this collection since last Diffusion step

private:
#ifdef USE_GL
  std::shared_ptr<GlState> mgl;
#endif
  float max_strength;
};

