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

// versions of the HO solver
#ifdef HOFORTRAN
#include "hofortran_interface.h"
#else
#include "dummysolver.h"
#endif

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

  // useful getter for the areas of participating elements
  const Vector<S>&             get_soln_area() const { return solnarea; }
  const Vector<S>&              get_gtop_wgt() const { return gtop_wgt; }
  const Vector<S>&              get_ptog_wgt() const { return ptog_wgt; }

  // return a Collection of the nodes on the open boundary
  const Points<S>& get_bc_nodes(const S _time) const { return open_p; }
  Points<S>&       get_bc_nodes(const S _time)       { return open_p; }

  // return a Collection of the solution nodes
  const Points<S>& get_vol_nodes(const S _time) const { return soln_p; }
  Points<S>&       get_vol_nodes(const S _time)       { return soln_p; }

  // if anyone needs to know the order of the geometric elements
  const int32_t get_geom_elem_order() const {
    // what kinds of elements do we have?
    const size_t nper = this->idx.size() / this->nb;
    // return order, where Lnod = order+1
    if (nper == 4) return 1;
    else if (nper == 9) return 2;
    else if (nper == 16) return 3;
    else if (nper == 25) return 4;
    else if (nper == 36) return 5;
    else {
      assert(false and "ERROR (get_geom_elem_order): Should not get here");
      return -1;
    }
  }

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

  // set vorticity at solution points, to use for reprojecting to particles
  void set_soln_vort(std::vector<double> _in) {
    assert(_in.size() > 0 && "ERROR (set_soln_vort): received zero solution points from solver");

    // convert input vector to local type
    Vector<S> vort(_in.begin(), _in.end());

    // pass it to the Collection of Points for use later
    soln_p.set_vort(vort);
  }

  //
  // generate and save a measure of the area of each *solution* node
  //
  void set_soln_areas() {

    // return if we do not need to recalculate these (vdelta changes)

    // how many solution nodes per geometric element?
    //const size_t nper = soln_p.get_n() / this->nb;
    //assert(nper*this->nb == soln_p.get_n() && "ERROR (set_mask_area_ho) element - solution node mismatch");

    //std::cout << "In set_soln_areas with weight mask size " << nper << std::endl;
    std::cout << "In set_soln_areas with " << soln_p.get_n() << " solution nodes" << std::endl;

    // the temporary array of doubles
    std::vector<double> hoarea(soln_p.get_n());

    // get an array of weights from the HO solver for a HO element
    //std::vector<double> wgt(nper);
#ifdef HOFORTRAN
    getsolnareas_d((int32_t)hoarea.size(), hoarea.data());
    //get_hoquad_weights_d((int32_t)nper, wgt.data());
#else
    //std::fill(wgt.begin(), wgt.end(), 1.0/(double)nper);
#endif
    //std::cout << "  first row of weight mask is ";
    //for (size_t j=0; j<std::sqrt(nper); ++j) std::cout << " " << wgt[j];
    //std::cout << std::endl;

    // copying Jacobian from the HO solver
    solnarea.resize(soln_p.get_n());
    std::copy(hoarea.begin(), hoarea.end(), solnarea.begin());

    // HACK using area here
    //solnarea.resize(soln_p.get_n());
    //for (size_t i=0; i<this->nb; ++i) {
    //  for (size_t j=0; j<nper; ++j) {
        // get quad areas from the parent class
    //    solnarea[nper*i+j] = (S)(this->area[i] * wgt[j]);
        //std::cout << "  " << i << " " << j << "  " << hoarea[nper*i+j] << " " << solnarea[nper*i+j] << std::endl;
    //  }
    //}
  }

  //
  // generate and save a scalar on each *solution* node corresponding to
  //   weights for grid-to-particle operation and particle-to-grid op
  //
  void set_overlap_weights(const S gtop_center, const S gtop_width,
                           const S ptog_center, const S ptog_width) {

    // get this array so we can reference it more easily
    const std::array<Vector<S>,Dimensions>& ptx = soln_p.get_pos();

    // size the arrays properly
    gtop_wgt.resize(soln_p.get_n());
    ptog_wgt.resize(soln_p.get_n());

    // now, based on the location of the solution nodes, zero out any
    //   that are too close to the body - HACK - assume D=1 cylinder
    for (size_t i=0; i<soln_p.get_n(); ++i) {

      // HACK - assume geometry is annular from 0.5 to 1.0
      const S thisrad = std::sqrt(std::pow(ptx[0][i],2)+std::pow(ptx[1][i],2));

      // grid-to-particle weight - typically peak in the middle and fall off near boundaries
      if (std::abs(thisrad-gtop_center) < gtop_width) {
        gtop_wgt[i] = 0.5 + 0.5*std::cos((1./gtop_width)*M_PI*(thisrad-gtop_center));
      } else {
        gtop_wgt[i] = 0.0;
      }

      // particle-to-grid weight - typically peak at the open BC
      if (std::abs(thisrad-ptog_center) < ptog_width) {
        ptog_wgt[i] = 0.5 + 0.5*std::cos((1./ptog_width)*M_PI*(thisrad-ptog_center));
      } else {
        ptog_wgt[i] = 0.0;
      }
    }

    std::cout << "  set_overlap_weights on " << soln_p.get_n() << " solution nodes" << std::endl;
  }

  //
  // return a particle version of the elements - one particle per *solution* node
  //
  ElementPacket<S> get_equivalent_particles(const Vector<S>& _circ, const S _vd) {

    assert(_circ.size() == soln_p.get_n() && "HOVolumes::get_equivalent_particles input vector size mismatch");

    // prepare the data arrays for the element packet
    std::vector<float> x;
    std::vector<Int> idx;
    std::vector<float> val;
    size_t thisn = 0;

    // get this array so we can reference it more easily
    const std::array<Vector<S>,Dimensions>& ptx = soln_p.get_pos();

    // loop over volume elements which have appreciable strength change
    for (size_t i=0; i<soln_p.get_n(); ++i) if (std::fabs(_circ[i]) > 1.e-7) {

      // convert input _circ and underlying geometry into particles
      // just one per element for starters
      x.push_back(ptx[0][i]);
      x.push_back(ptx[1][i]);
      val.push_back(_circ[i]);
      thisn++;

      // use size of element and particle radius to ensure overlap
    }

    return ElementPacket<S>(x, idx, val, thisn, (uint8_t)0);
  }
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
  Vector<S>             solnarea;   // area of each solution node
  Vector<S>             gtop_wgt;   // grid-to-particle weight (per solution node)
  Vector<S>             ptog_wgt;   // particle-to-grid weight (per solution node)

  // boundary elements (only used when this is an euler/hybrid Collection)
  Surfaces<S>             wall_s;   // wall boundary surface (from msh file)
  Surfaces<S>             open_s;   // open boundary surface (from msh file)
  Points<S>               open_p;   // solution nodes at open boundary (from HO Solver)
  Points<S>               soln_p;   // solution nodes (from HO Solver)

private:

};

