/*
 * HOVolumes.h - Specialized class for high-order volumes
 *
 * (c)2020-1 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
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
#include "Volumes.h"
#include "Surfaces.h"
#include "Points.h"
#include "Reflect.h"

// versions of the HO solver
#ifdef HOFORTRAN
#include "hofortran_interface.h"
#elif HOCXX
#include "HO_2D.hpp"
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
  HOVolumes(const ElementPacket<float>& _velems,
            const ElementPacket<float>& _welems,
            const ElementPacket<float>& _oelems,
            const elem_t _e,
            const move_t _m,
            std::shared_ptr<Body> _bp)
    : Volumes<S>(_velems, _e, _m, _bp),
      wall_s(Surfaces<S>(_welems, _e, _m, _bp)),
      open_s(Surfaces<S>(_oelems, _e, _m, _bp)),
      open_p(Points<S>(ElementPacket<float>((uint8_t)0), _e, _m, _bp, 0.0)),
      soln_p(Points<S>(ElementPacket<float>((uint8_t)0), _e, _m, _bp, 0.0)),
      inlet_s(Surfaces<S>(ElementPacket<float>((uint8_t)1), _e, _m, _bp)),
      outlet_s(Surfaces<S>(ElementPacket<float>((uint8_t)1), _e, _m, _bp))
      //reduced_p(Points<S>(ElementPacket<S>((uint8_t)0), _e, _m, _bp, 0.0))
#ifdef HOCXX
      ,solver()
#endif
  { std::cout << "  in HOVolumes()" << std::endl; }

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

#ifdef HOCXX
  HO_2D& get_ho_solver() { return solver; }
#endif

  const bool have_inlet() const  { return (inlet_s.get_npanels()  > 0) ? true : false; }
  const bool have_outlet() const { return (outlet_s.get_npanels() > 0) ? true : false; }

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
    const std::vector<Int>& tempgeom = this->get_idx();
    std::vector<uint32_t> returnidx(tempgeom.begin(), tempgeom.end());
    return returnidx;
  }
  std::vector<uint32_t> get_wall_idx() {
    const std::vector<Int>& tempgeom = wall_s.get_idx();
    std::vector<uint32_t> returnidx(tempgeom.begin(), tempgeom.end());
    return returnidx;
  }
  std::vector<uint32_t> get_open_idx() {
    const std::vector<Int>& tempgeom = open_s.get_idx();
    std::vector<uint32_t> returnidx(tempgeom.begin(), tempgeom.end());
    return returnidx;
  }
  std::vector<uint32_t> get_inlet_idx() {
    const std::vector<Int>& tempgeom = inlet_s.get_idx();
    std::vector<uint32_t> returnidx(tempgeom.begin(), tempgeom.end());
    return returnidx;
  }
  std::vector<uint32_t> get_outlet_idx() {
    const std::vector<Int>& tempgeom = outlet_s.get_idx();
    std::vector<uint32_t> returnidx(tempgeom.begin(), tempgeom.end());
    return returnidx;
  }

  // return a new vector of the inlet/outlet velocities
  std::vector<double> get_inlet_vel() {
    const Vector<S>& normvel = inlet_s.get_norm_bcs();
    const std::vector<double> retvels(normvel.begin(), normvel.end());
    return retvels;
  }
  std::vector<double> get_outlet_vel() {
    const Vector<S>& normvel = outlet_s.get_norm_bcs();
    const std::vector<double> retvels(normvel.begin(), normvel.end());
    return retvels;
  }

  // retrieve solution points from solver, set them here
  void set_soln_pts(std::vector<double> _in) {
    assert(_in.size() > 0 && "ERROR (set_soln_pts): received zero solution points from solver");

    // convert input vector to local coords
    std::vector<float> x(_in.begin(), _in.end());

    // and dummy val array with zeros
    std::vector<float> val(x.size()/2);

    // ready to pass a packet to the ctor
    soln_p.add_new(ElementPacket<float>(x, std::vector<Int>(), val, x.size()/2, 0), 0.0);
  }

  // retrieve solution points on open boundary from solver, set them here
  void set_open_pts(std::vector<double> _in) {
    assert(_in.size() > 0 && "ERROR (set_open_pts): received zero open bc points from solver");

    // convert input vector to local coords
    std::vector<float> x(_in.begin(), _in.end());

    // and dummy val array with zeros
    std::vector<float> val(x.size()/2);

    // ready to pass a packet to the ctor
    open_p.add_new(ElementPacket<float>(x, std::vector<Int>(), val, x.size()/2, 0), 0.0);
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

  // append a Surfaces to the object for the inlet
  void add_inlet(const ElementPacket<float>& _in) {

    // ensure that this packet really is Surfaces
    assert(_in.idx.size() != 0 && "Input ElementPacket is empty");
    assert(_in.ndim == 1 && "Input ElementPacket is not Surfaces");

    inlet_s = Surfaces<S>(_in, this->E, this->M, this->B);
  }

  // append a Surfaces to the object for the outlet
  void add_outlet(const ElementPacket<float>& _in) {

    // ensure that this packet really is Surfaces
    assert(_in.idx.size() != 0 && "Input ElementPacket is empty");
    assert(_in.ndim == 1 && "Input ElementPacket is not Surfaces");

    outlet_s = Surfaces<S>(_in, this->E, this->M, this->B);
  }

  // need to ensure that inlet flows equal outlet flows, I believe?
  void conserve_iolet_volume() {
    std::cout << "  in conserve_iolet_volume with " << inlet_s.get_npanels() << " and " << outlet_s.get_npanels() << " panels" << std::endl;

    // only do this if there are inlet AND outlet surfaces
    if (inlet_s.get_npanels() == 0 or outlet_s.get_npanels() == 0) return;

    const float inrate = inlet_s.get_total_inflow();
    const float outrate = outlet_s.get_total_outflow();

    std::cout << "  total HOVolume inflow, outflow " << inrate << " " << outrate << std::endl;

    if (inrate > std::numeric_limits<float>::epsilon()) {
      if (outrate > std::numeric_limits<float>::epsilon()) {
        // there is finite inflow and outflow, ensure their magnitudes are matched

        std::cout << "    scaling outflows by " << inrate/outrate << std::endl;

        // correct the outflow to match
        outlet_s.scale_outflow(inrate/outrate);

      } else {
        // there is inflow, but zero outflow - set outflow to accept all inflow
        // note: this is different than what we do in Simulation::conserve_iolet_volume() !

        std::cout << "    outflow was 0 - setting to match inflows" << std::endl;

        const float outarea = outlet_s.get_total_arclength();
        outlet_s.set_outflow(inrate/outarea);
      }
    }
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
#ifdef HOFORTRAN
  void set_soln_areas() {
#elif HOCXX
  void set_soln_areas(HO_2D& solver) {
#endif

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
#elif HOCXX
    solver.getsolnareas_d((int32_t)hoarea.size(), hoarea.data());
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

    std::cout << "In set_overlap_weights" << std::endl;

    // get this array so we can reference it more easily
    //const std::array<Vector<S>,Dimensions>& ptx = soln_p.get_pos();

    // size the arrays properly
    gtop_wgt.resize(soln_p.get_n());
    ptog_wgt.resize(soln_p.get_n());

    // get distances from test points to outer boundary
    Vector<S> dist_to_open = get_nearest_distances(open_s, soln_p);
    Vector<S> dist_to_wall = get_nearest_distances(wall_s, soln_p);

    // now, based on the location of the solution nodes, set the weights
    for (size_t i=0; i<soln_p.get_n(); ++i) {

      // location normalized by range, values are always [0..1]
      const S normrad = dist_to_wall[i] / (dist_to_wall[i] + dist_to_open[i]);

      // grid-to-particle weight - typically peak in the middle and fall off near boundaries
      if (std::abs(normrad-gtop_center) < gtop_width) {
        gtop_wgt[i] = 0.5 + 0.5*std::cos((1.0/gtop_width)*M_PI*(normrad-gtop_center));
      } else {
        gtop_wgt[i] = 0.0;
      }

      // particle-to-grid weight - typically peak at the open BC
      if (std::abs(normrad-ptog_center) < ptog_width) {
        ptog_wgt[i] = 0.5 + 0.5*std::cos((1.0/ptog_width)*M_PI*(normrad-ptog_center));
      } else {
        ptog_wgt[i] = 0.0;
      }

      //std::cout << "soln pt " << i << " at " << ptx[0][i] << " " << ptx[1][i] << " has wgts " << gtop_wgt[i] << " " << ptog_wgt[i] << std::endl;
      //std::cout << std::fixed << std::setprecision(3) << "soln pt " << i << " at " << ptx[0][i] << " " << ptx[1][i] << " has dists " << dist_to_open[i] << " " << dist_to_wall[i] << " and wgts " << gtop_wgt[i] << " " << ptog_wgt[i] << std::endl;
    }

    std::cout << "    set_overlap_weights on " << soln_p.get_n() << " solution nodes" << std::endl;
  }

  //
  // return a particle version of the elements - one particle per *solution* node
  // new particles have given vdelta core size
  // only generate new particles with strength greater than a threshold
  //
  ElementPacket<float> get_equivalent_particles(const Vector<S>& _circ, const S _vd, const S _thresh) {

    assert(_circ.size() == soln_p.get_n() && "HOVolumes::get_equivalent_particles input vector size mismatch");

    // prepare the data arrays for the element packet (there's an "idx" in Volumes)
    std::vector<float> _x;
    std::vector<Int> _idx;
    std::vector<float> _val;
    size_t thisn = 0;

    // get this array so we can reference it more easily
    const std::array<Vector<S>,Dimensions>& ptx = soln_p.get_pos();

    // loop over volume elements which have appreciable strength change
    for (size_t i=0; i<soln_p.get_n(); ++i) if (std::fabs(_circ[i]) > _thresh) {

      // convert input _circ and underlying geometry into particles
      // just one per element for starters
      _x.push_back(ptx[0][i]);
      _x.push_back(ptx[1][i]);
      _val.push_back(_circ[i]);
      thisn++;

      // use size of element and particle radius to ensure overlap
    }

    return ElementPacket<float>(_x, _idx, _val, thisn, (uint8_t)0);
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
  //Points<S>            reduced_p;   // reduced set of solution nodes

  // optional boundary elements (one inlet and one outlet for now)
  Surfaces<S>            inlet_s;	// inlet boundary surface (from msh file)
  Surfaces<S>           outlet_s;	// outlet boundary surface (from msh file)

  // the HO Solver
#ifdef HOFORTRAN
  // none needed
#elif HOCXX
  HO_2D solver;
#endif

private:

};

