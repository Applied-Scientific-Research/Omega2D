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
  HOVolumes(const ElementPacket<S>& _velems,
            const ElementPacket<S>& _welems,
            const ElementPacket<S>& _oelems,
            const elem_t _e,
            const move_t _m,
            std::shared_ptr<Body> _bp)
    : Volumes<S>(_velems, _e, _m, _bp),
      wall_s(Surfaces<S>(_welems, _e, _m, _bp)),
      open_s(Surfaces<S>(_oelems, _e, _m, _bp)),
      open_p(Points<S>(ElementPacket<S>((uint8_t)0), _e, _m, _bp, 0.0)),
      soln_p(Points<S>(ElementPacket<S>((uint8_t)0), _e, _m, _bp, 0.0))
  { }

  // return a Collection of the nodes on the open boundary
  const Points<S>& get_bc_nodes(const S _time) const { return open_p; }
  Points<S>&       get_bc_nodes(const S _time)       { return open_p; }

  // return a Collection of the solution nodes
  const Points<S>& get_vol_nodes(const S _time) const { return soln_p; }
  Points<S>&       get_vol_nodes(const S _time)       { return soln_p; }

  // return a new vector of the node locations in the geometry
  std::vector<double> get_node_pos() {
    // from array of vectors
    std::array<Vector<S>,Dimensions>& nodegeom = this->get_pos();
    // to single vector
    std::vector<double> allpos(Dimensions*nodegeom[0].size());
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<nodegeom[d].size(); ++i) {
        allpos[Dimensions*i+d] = nodegeom[d][i];
      }
    }
    return allpos;
  }

  // return a new vector of the indices of the geometry
  std::vector<uint32_t> get_elem_idx() {
    const std::vector<Int>& elemgeom = this->get_idx();
    std::vector<uint32_t> elemidx(elemgeom.begin(), elemgeom.end());
    return elemidx;
  }
  std::vector<uint32_t> get_wall_idx() {
    const std::vector<Int>& wallgeom = wall_s.get_idx();
    std::vector<uint32_t> wallidx(wallgeom.begin(), wallgeom.end());
    return wallidx;
  }
  std::vector<uint32_t> get_open_idx() {
    const std::vector<Int>& opengeom = open_s.get_idx();
    std::vector<uint32_t> openidx(opengeom.begin(), opengeom.end());
    return openidx;
  }

  // retrieve solution points from solver, set them here
  void set_soln_pts(std::vector<double> _in) {
    assert(_in.size() > 0 && "ERROR (set_soln_pts): received zero solution points from solver");

    // convert input vector to local coords
    std::vector<S> x(_in.begin(), _in.end());

    // and dummy val array with zeros
    std::vector<S> val(x.size()/2);

    // ready to pass a packet to the ctor
    soln_p.add_new(ElementPacket<S>(x, std::vector<Int>(), val, x.size()/2, 0), 0.0);
  }

  // retrieve solution points on open boundary from solver, set them here
  void set_open_pts(std::vector<double> _in) {
    assert(_in.size() > 0 && "ERROR (set_open_pts): received zero open bc points from solver");

    // convert input vector to local coords
    std::vector<S> x(_in.begin(), _in.end());

    // and dummy val array with zeros
    std::vector<S> val(x.size()/2);

    // ready to pass a packet to the ctor
    open_p.add_new(ElementPacket<S>(x, std::vector<Int>(), val, x.size()/2, 0), 0.0);
  }

  // append a Surfaces to the object for the wall-boundary geometry
  void add_wall(const ElementPacket<float>& _in) {

    // ensure that this packet really is Surfaces
    assert(_in.idx.size() != 0 && "Input ElementPacket is empty");
    assert(_in.ndim == 1 && "Input ElementPacket is not Surfaces");

    assert("We are not yet checking for matching element orders!");

    wall_s = Surfaces<S>(_in, this->E, this->M, this->B);

    // and that it has the right number of values per element/node
    //if (this->E == inert) assert(_in.val.size() == 0 && "Input ElementPacket with inert Volumes has nonzero val array");
    //else assert(_in.val.size() == _in.nelem && "Input ElementPacket with Volumes has bad val array size");
  }

  // append a Surfaces to the object for the open-boundary geometry
  void add_open(const ElementPacket<float>& _in) {

    // ensure that this packet really is Surfaces
    assert(_in.idx.size() != 0 && "Input ElementPacket is empty");
    assert(_in.ndim == 1 && "Input ElementPacket is not Surfaces");

    assert("We are not yet checking for matching element orders!");

    open_s = Surfaces<S>(_in, this->E, this->M, this->B);

    // and that it has the right number of values per element/node
    //if (this->E == inert) assert(_in.val.size() == 0 && "Input ElementPacket with inert Volumes has nonzero val array");
    //else assert(_in.val.size() == _in.nelem && "Input ElementPacket with Volumes has bad val array size");
  }


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

