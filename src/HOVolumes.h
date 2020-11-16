/*
 * HOVolumes.h - Specialized class for high-order volumes
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "Volumes.h"
#include "Surfaces.h"
#include "Points.h"

#include <iostream>
#include <vector>
#include <array>
#include <algorithm> // for max_element
#include <optional>
#include <cassert>


// 2D elements in 2D
template <class S>
class HOVolumes: public Volumes<S> {
public:

  // use parent constructor
  HOVolumes(const ElementPacket<S>& _elems,
            const elem_t _e,
            const move_t _m,
            std::shared_ptr<Body> _bp)
    : Volumes(_elems, _e, _m, _bp)
  { }

  // return a Collection of the nodes on the open boundary
  const Points<S>& get_bc_nodes(const S _time) const { return open_p; }
  Points<S>&       get_bc_nodes(const S _time)       { return open_p; }

  // return a Collection of the solution nodes
  const Points<S>& get_vol_nodes(const S _time) const { return soln_p; }
  Points<S>&       get_vol_nodes(const S _time)       { return soln_p; }


  // add more nodes and elements to this collection
  void add_new(const std::vector<S>&   _x,
               const std::vector<Int>& _idx,
               const std::vector<S>&   _val) {

    // remember old sizes of nodes and element arrays
    const size_t nnold = this->n;
    const size_t neold = get_nelems();

    // make sure input arrays are correctly-sized

    // assume all elements are 1st order quads (4 corner indices)
    const size_t nper = 4;
    assert(_idx.size() % nper == 0 && "Index array is not an even multiple of 4");
    const size_t nelems = _idx.size() / nper;
    // if no surfs, quit out now
    if (nelems == 0) return;

    assert(_x.size() % Dimensions == 0 && "Position array is not an even multiple of dimensions");
    const size_t nnodes = _x.size() / Dimensions;

    assert(_val.size() % nnodes == 0 && "Value array is not an even multiple of node count");

    std::cout << "  adding " << nelems << " new elements and " << nnodes << " new points to collection..." << std::endl;

    // DON'T call the method in the base class, because we do things differently here
    //ElementBase<S>::add_new(_in);

    // pull out the node locations, they are base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(nnold+nnodes);
      for (size_t i=0; i<nnodes; ++i) {
        this->x[d][nnold+i] = _x[Dimensions*i+d];
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
    idx.resize(4*neold + _idx.size());
    for (size_t i=0; i<nper*nelems; ++i) {
      // make sure it exists in the nodes array
      if (_idx[i] >= nnold+nnodes) idx_are_all_good = false;
      idx[nper*neold+i] = nnold + _idx[i];
    }
    assert(idx_are_all_good && "Some indicies are bad");

    // compute all basis vectors and element areas
    compute_bases(neold+nelems);

    // now, depending on the element type, put the value somewhere - but element-wise, so here
    if (this->E == active) {
/*
      // value is a fixed strength for the element
      ps[0]->reserve(neold+nelems); 
      ps[0]->insert(ps[0]->end(), _val.begin(), _val.end());
      // HACK - should use the size of _val to determine whether we have data here
      ps[1]->reserve(neold+nelems); 
      ps[1]->insert(ps[1]->end(), _val.begin(), _val.end());
*/

    } else if (this->E == reactive) {
/*
      // value is a boundary condition
      bc[0]->reserve(neold+nelems); 
      bc[0]->insert(bc[0]->end(), _val.begin(), _val.end());
      // HACK - should use the size of _val to determine whether we have data here
      bc[1]->reserve(neold+nelems); 
      bc[1]->insert(bc[1]->end(), _val.begin(), _val.end());

      // upsize vortex sheet and raw strength arrays, too
      for (size_t d=0; d<2; ++d) {
        if (ps[d]) {
          ps[d]->resize(neold+nelems);
          //std::fill(ps[d]->begin(), ps[d]->end(), 0.0);
        }
      }
*/

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnold+nnodes);
    }

    // element velocity is here
/*
    for (size_t d=0; d<Dimensions; ++d) {
      pu[d].resize(neold+nelems);
    }
*/

    // debug print
    if (false) {
      std::cout << "Nodes" << std::endl;
      for (size_t i=0; i<nnold+nnodes; ++i) {
        std::cout << "  " << i << " " << this->x[0][i] << " " << this->x[1][i] << std::endl;
      }
      std::cout << "Elems" << std::endl;
      for (size_t i=0; i<neold+nelems; ++i) {
        std::cout << "  " << i << " " << idx[4*i] << " " << idx[4*i+1] << " " << idx[4*i+2] << " " << idx[4*i+3] << std::endl;
      }
    }

    // need to reset the base class n
    this->n += nnodes;
    nb += nelems;

    // re-find geometric center
/*
    if (this->M == bodybound) {
      set_geom_center();
    }
*/
  }

  // append a Surfaces to the object
  void add_wall(const ElementPacket<float>& _in) {

    // ensure that this packet really is Surfaces
    assert(_in.idx.size() != 0 && "Input ElementPacket is not Volumes");
    assert(_in.ndim == 2 && "Input ElementPacket is not Volumes");

    assert("We are not yet checking for matching element orders!");

    wall_s = _in;

    // and that it has the right number of values per element/node
    //if (this->E == inert) assert(_in.val.size() == 0 && "Input ElementPacket with inert Volumes has nonzero val array");
    //else assert(_in.val.size() == _in.nelem && "Input ElementPacket with Volumes has bad val array size");
  }


/*
  // up-size all arrays to the new size, filling with sane values
  void resize(const size_t _nnew) {
    const size_t currn = this->n;
    //std::cout << "  inside Volumes::resize with " << currn << " " << _nnew << std::endl;

    // must explicitly call the method in the base class - this sets n
    ElementBase<S>::resize(_nnew);

    if (_nnew == currn) return;

    // radii here
    const size_t thisn = r.size();
    r.resize(_nnew);
    for (size_t i=thisn; i<_nnew; ++i) {
      r[i] = 1.0;
    }
  }
*/

  //
  // return a particle version of the elements (useful during Diffusion)
  // offset is in world units - NOT scaled
  //
/*
  std::vector<S> represent_elems_as_particles(const S _offset, const S _vdelta) {

    // how many elements?
    const size_t num_pts = get_nelems();

    // init the output vector (x, y, s, r)
    std::vector<S> px(num_pts*4);

    // get basis vectors
    std::array<Vector<S>,2>& norm = b[1];

    // the fluid is to the left walking from one point to the next
    // so go CW around an external boundary starting at theta=0 (+x axis)

    for (size_t i=0; i<num_pts; i++) {
      Int id0 = idx[2*i];
      Int id1 = idx[2*i+1];
      // start at center of element
      px[4*i+0] = 0.5 * (this->x[0][id1] + this->x[0][id0]);
      px[4*i+1] = 0.5 * (this->x[1][id1] + this->x[1][id0]);
      //std::cout << "  element center is " << px[4*i+0] << " " << px[4*i+1];
      // push out a fixed distance
      // this assumes properly resolved, vdelta and dt
      px[4*i+0] += _offset * norm[0][i];
      px[4*i+1] += _offset * norm[1][i];
      // the element strength is the solved strength plus the boundary condition
      float this_str = (*ps[0])[i];
      // add on the (vortex) bc value here
      if (this->E == reactive) this_str += (*bc[0])[i];
      // complete the element with a strength
      px[4*i+2] = this_str * area[i];
      // and the core size
      px[4*i+3] = _vdelta;
      //std::cout << "  new part is " << px[4*i+0] << " " << px[4*i+1] << " " << px[4*i+2] << " " << px[4*i+3] << std::endl;
    }

    return px;
  }
*/

protected:
  // ElementBase.h has x, s, u, ux on the *nodes*
  // Volumes.h has nb, idx, area

  // boundary elements (only used when this is an euler/hybrid Collection)
  Surfaces<S>             wall_s;   // wall boundary surface (from msh file)
  Surfaces<S>             open_s;   // open boundary surface (from msh file)
  Points<S>               open_p;   // solution nodes at open boundary (from HO Solver)
  Points<S>               soln_p;   // solution nodes (from HO Solver)

private:

};

