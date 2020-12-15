/*
 * Hybrid.h - coordinate with an external Eulerian flow solver to compute near-body flow
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Collection.h"
#include "Convection.h"
#include "BEM.h"
#include "Merge.h"

// versions of the HO solver
#ifdef HOFORTRAN
#include "hofortran_interface.h"
#else
#include "dummysolver.h"
#endif

#include <iostream>
#include <vector>
//#include <numeric>	// for transform_reduce (C++17)


//
// Contain and process the hybrid solver
//
// templatized on 'S'torage, 'A'ccumulator (calculation) types, and element 'I'ndex types
//
template <class S, class A, class I>
class Hybrid {
public:
  Hybrid()
    : active(false),
      initialized(false),
      elementOrder(1),
      timeOrder(1),
      numSubsteps(100),
      preconditioner("none"),
      solverType("fgmres")
#ifndef HOFORTRAN
      ,solver()
#endif
      //vrm(),
      //h_nu(0.1)
    {}

  const bool is_active() const { return active; }
  void set_active(const bool _do_hybrid) { active = _do_hybrid; }
  void activate() { active = true; }
  void deactivate() { active = false; }

  void init( std::vector<HOVolumes<S>>&);
  void reset();
  void first_step( const double,
             const std::array<double,Dimensions>&,
             std::vector<Collection>&,
             std::vector<Collection>&,
             BEM<S,I>&,
             Convection<S,A,I>&,
             std::vector<HOVolumes<S>>&);
  void step( const double,
             const double,
             const float,
             const std::array<double,Dimensions>&,
             std::vector<Collection>&,
             std::vector<Collection>&,
             BEM<S,I>&,
             Convection<S,A,I>&,
             std::vector<HOVolumes<S>>&,
             const float);

  // read/write parameters
  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

  void draw_advanced();

private:
  // are we even using the hybrid scheme?
  bool active;
  bool initialized;

  // parameters from json for the solver
  int elementOrder;
  int timeOrder;
  int numSubsteps;
  std::string preconditioner;
  std::string solverType;

  // the HO Solver
#ifndef HOFORTRAN
  DummySolver::Solver solver;
#endif

  // local copies of particle data
  //Particles<S> temp;

  // might need vrm-like solver to redistribute strengths
  //VRM<S,double,2> redistributor;

  // execution environment for velocity summations (not BEM)
  //ExecEnv hyb_env;
};


//
// Initialize external high-order (HO) solver
//
template <class S, class A, class I>
void Hybrid<S,A,I>::init(std::vector<HOVolumes<S>>& _euler) {
  std::cout << "Inside Hybrid::init with " << _euler.size() << " volumes" << std::endl;
  assert(_euler.size() == 1 && "ERROR Hybrid::init does not support more than 1 hybrid volume");

#ifdef HOFORTRAN
  // set in the Fortran solver
  (void) set_defaults();
  (void) enable_hybrid();
  (void) set_elemorder((int32_t)elementOrder);
#else
  // set in the dummy C++ solver
  solver.set_elemorder_d((int32_t)elementOrder);
#endif

  for (auto &coll : _euler) {
    // transform to current position
    coll.move(0.0, 0.0);

#ifdef HOFORTRAN
    // and, set the mesh in the Fortran solver
    {
    // make temporary vectors to convert data types
    std::vector<double> nodes_as_dble = coll.get_node_pos();

    std::vector<uint32_t> elemidxu = coll.get_elem_idx();
    std::vector<int32_t> elemidx(elemidxu.begin(), elemidxu.end());

    std::vector<uint32_t> wallidxu = coll.get_wall_idx();
    std::vector<int32_t> wallidx(wallidxu.begin(), wallidxu.end());

    std::vector<uint32_t> openidxu = coll.get_open_idx();
    std::vector<int32_t> openidx(openidxu.begin(), openidxu.end());

    // now we can call this - HACK - need to find order of geom mesh
    (void) load_mesh_arrays_d((int32_t)coll.get_geom_elem_order(),
                        (int32_t)nodes_as_dble.size(), nodes_as_dble.data(),
                        (int32_t)elemidx.size(), elemidx.data(),
                        (int32_t)wallidx.size(), wallidx.data(),
                        (int32_t)openidx.size(), openidx.data());
    }
#else
    // call the external solver with the current geometry
    // this will calculate the Jacobian and other cell-specific properties
    (void) solver.init_d(coll.get_node_pos(), coll.get_elem_idx(),
                         coll.get_wall_idx(), coll.get_open_idx());
#endif

#ifdef HOFORTRAN
    // ask fort solver for the open BC solution nodes, and the full internal solution nodes
    {
      // again, since fortran is dumb, we need extra steps
      int32_t solnptlen = getsolnptlen();
      std::cout << "There are " << (solnptlen/2) << " solution nodes" << std::endl;
      std::vector<double> solnpts(solnptlen);
      (void) getsolnpts_d(solnptlen, solnpts.data());
      std::cout << "  First soln point is at " << solnpts[0] << " " << solnpts[1] << std::endl;
      std::cout << "  Secnd soln point is at " << solnpts[2] << " " << solnpts[3] << std::endl;
      coll.set_soln_pts(solnpts);

      // repeat for open boundary nodes
      solnptlen = getopenptlen();
      std::cout << "There are " << (solnptlen/2) << " open boundary solution nodes" << std::endl;
      solnpts.resize(solnptlen);
      (void) getopenpts_d(solnptlen, solnpts.data());
      std::cout << "  First open point is at " << solnpts[0] << " " << solnpts[1] << std::endl;
      std::cout << "  Secnd open point is at " << solnpts[2] << " " << solnpts[3] << std::endl;
      coll.set_open_pts(solnpts);
    }
#else
    // and the dummy for the open BC solution nodes, and the full internal solution nodes
    coll.set_open_pts(solver.getopenpts_d());
    coll.set_soln_pts(solver.getsolnpts_d());
#endif

  }

  initialized = true;
}


//
// Simulation is reset - do we need to clear the Euler grid initialization?
//
template <class S, class A, class I>
void Hybrid<S,A,I>::reset() {
  initialized = false;

#ifdef HOFORTRAN
  // clear out memory of the grid
  (void) clean_up();
#endif
}


//
// Send first set of velocities to solver
//
template <class S, class A, class I>
void Hybrid<S,A,I>::first_step(const double                   _time,
                         const std::array<double,Dimensions>& _fs,
                         std::vector<Collection>&             _vort,
                         std::vector<Collection>&             _bdry,
                         BEM<S,I>&                            _bem,
                         Convection<S,A,I>&                   _conv,
                         std::vector<HOVolumes<S>>&           _euler) {

  if (not active) return;
  if (not initialized) init(_euler);

  std::cout << "Inside Hybrid::first_step at t=" << _time << std::endl;

  //
  // solve for velocity at each open-boundary solution node
  //

  // make a list of euler boundary regions
  std::vector<Collection> euler_bdrys;
  for (auto &coll : _euler) {
    // transform to current position
    coll.move(_time, 0.0);

    // isolate open/outer boundaries
    euler_bdrys.emplace_back(coll.get_bc_nodes(_time));
  }

  // get vels and vorts on each euler region - and force it
  _conv.find_vels(_fs, _vort, _bdry, euler_bdrys, velonly, true);

  for (auto &coll : euler_bdrys) {
    // convert to transferable packet
    std::array<Vector<S>,Dimensions> openvels = std::visit([=](auto& elem) { return elem.get_vel(); }, coll);

    // convert to a std::vector<double>
    std::vector<double> packedvels(Dimensions*openvels[0].size());
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<openvels[d].size(); ++i) {
        packedvels[Dimensions*i+d] = openvels[d][i];
      }
    }

    // transfer BC packet to solver
#ifdef HOFORTRAN
    (void) setopenvels_d((int32_t)packedvels.size(), packedvels.data());
#else
    (void) solver.setopenvels_d(packedvels);
#endif
  }

  // get vels and vorts on each euler region - and force it
  _conv.find_vort(_vort, _bdry, euler_bdrys);

  for (auto &coll : euler_bdrys) {
    // now prepare the open boundary vorticity values
    Vector<S> volvort = std::visit([=](auto& elem) { return elem.get_vort(); }, coll);
    std::vector<double> vorts(volvort.begin(), volvort.end());

    // transfer BC packet to solver
#ifdef HOFORTRAN
    (void) setopenvort_d((int32_t)vorts.size(), vorts.data());
#else
    // nothing here
#endif
  }

  //
  // now do the same for the vorticity at each solution node
  //

  // make a list of euler volume regions
  std::vector<Collection> euler_vols;
  for (auto &coll : _euler) {
    // transform to current position
    coll.move(_time, 0.0);

    // isolate open/outer boundaries
    euler_vols.emplace_back(coll.get_vol_nodes(_time));
  }

  // get vorts on each euler region - and force it
  _conv.find_vort(_vort, _bdry, euler_vols);

  for (auto &coll : euler_vols) {
    // convert to transferable packet
    Vector<S> volvort = std::visit([=](auto& elem) { return elem.get_vort(); }, coll);

    // convert to a std::vector<double>
    std::vector<double> vorts(volvort.begin(), volvort.end());

    // transfer BC packet to solver
#ifdef HOFORTRAN
    (void) setsolnvort_d((int32_t)vorts.size(), vorts.data());
#else
    (void) solver.setsolnvort_d(vorts);
#endif
  }
}

//
// Forward integration step
//
template <class S, class A, class I>
void Hybrid<S,A,I>::step(const double                         _time,
                         const double                         _dt,
                         const float                          _re,
                         const std::array<double,Dimensions>& _fs,
                         std::vector<Collection>&             _vort,
                         std::vector<Collection>&             _bdry,
                         BEM<S,I>&                            _bem,
                         Convection<S,A,I>&                   _conv,
                         std::vector<HOVolumes<S>>&           _euler,
                         const float                          _vd) {

  if (not active) return;
  if (not initialized) init(_euler);

  std::cout << "Inside Hybrid::step at t=" << _time << " and dt=" << _dt << std::endl;

  //
  // part A - prepare BCs for Euler solver
  //

  // update the BEM solution
  solve_bem<S,A,I>(_time, _fs, _vort, _bdry, _bem);

  // make a list of euler boundary regions
  std::vector<Collection> euler_bdrys;
  for (auto &coll : _euler) {
    // transform to current position
    coll.move(_time, 0.0);

    // isolate open/outer boundaries
    euler_bdrys.emplace_back(coll.get_bc_nodes(_time));
  }

  // get vels on each euler region - and force it
  _conv.find_vels(_fs, _vort, _bdry, euler_bdrys, velonly, true);

  for (auto &coll : euler_bdrys) {
    // convert to transferable packet
    std::array<Vector<S>,Dimensions> openvels = std::visit([=](auto& elem) { return elem.get_vel(); }, coll);

    // convert to a std::vector<double>
    std::vector<double> packedvels(Dimensions*openvels[0].size());
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<openvels[d].size(); ++i) {
        packedvels[Dimensions*i+d] = openvels[d][i];
      }
    }

    // transfer BC packet to solver
#ifdef HOFORTRAN
    (void) setopenvels_d((int32_t)packedvels.size(), packedvels.data());
#else
    (void) solver.setopenvels_d(packedvels);
#endif
  }

  // get vorts on each euler region
  _conv.find_vort(_vort, _bdry, euler_bdrys);

  for (auto &coll : euler_bdrys) {
    // now prepare the open boundary vorticity values
    Vector<S> volvort = std::visit([=](auto& elem) { return elem.get_vort(); }, coll);
    std::vector<double> vorts(volvort.begin(), volvort.end());

    // transfer BC packet to solver
#ifdef HOFORTRAN
    (void) setopenvort_d((int32_t)vorts.size(), vorts.data());
#else
    // nothing here
#endif
  }


  //
  // part B - call Euler solver
  //

  // call solver - solves all Euler volumes at once?
#ifdef HOFORTRAN
  (void) solveto_d((double)_dt, (int32_t)numSubsteps, (int32_t)timeOrder, (double)_re);
#else
  (void) solver.solveto_d((double)_time, (int32_t)numSubsteps, (int32_t)timeOrder, (double)_re);
#endif

  // pull results from external solver (assume just one for now)
  //for (auto &coll : _euler) {
  //  std::vector<double> allvorts = solver.getallvorts_d();
  //  assert(allvorts.size() == coll.get_vol_nodes(_time).get_n() && "ERROR (Hybrid::step) vorticity from solver is not the right size");

    // assign to the HOVolume
  //  coll.set_soln_vort(allvorts);
  //}

  // convert vorticity results to drawable and writable Volume elements (OpenGL stuff) - later


  //
  // part C - update particle strengths accordionly
  //

  std::cout << "Inside Hybrid::step updating particle strengths" << std::endl;

  // here's one way to do it:
  // identify all free vortex particles inside of euler regions and remove them
  //   (add up how much circulation we remove)
  //   then we'd need to re-run BEM

  // another way:
  // don't remove any existing particles, just add new ones over the old ones
  //   and let merge take care of the extra density - this is what we do here:

  for (auto &coll : _euler) {

    Points<S>& solnpts = coll.get_vol_nodes(_time);
    const size_t thisn = solnpts.get_n();

    // pull results from external solver (assume just one for now)
    std::vector<double> eulvort;
#ifdef HOFORTRAN
    {
      // again, since fortran is dumb, we need extra steps
      int32_t solnptlen = getsolnptlen() / 2;
      std::cout << "There are " << solnptlen << " solution nodes" << std::endl;
      eulvort.resize(solnptlen);
      (void) getallvorts_d(solnptlen, eulvort.data());
      //std::cout << "  vort 2014 from solver " << eulvort[2013] << std::endl;
      //std::cout << "  vorts from solver " << eulvort[0] << " " << eulvort[1] << " " << eulvort[2] << std::endl;
      //std::cout << "               more " << eulvort[3] << " " << eulvort[4] << " " << eulvort[5] << std::endl;
    }
#else
    eulvort = solver.getallvorts_d();
#endif
    assert(eulvort.size() == thisn && "ERROR (Hybrid::step) vorticity from solver is not the right size");

    if (true) {
      std::cout << "  vorticity from eulerian " << std::endl;
      const auto& locs = solnpts.get_pos();
      for (size_t i=0; i<thisn; ++i) {
        if (std::abs(locs[1][i]/locs[0][i]-0.4868) < 0.01 && locs[0][i] > 0.0) {
          const S dist = std::sqrt(locs[0][i]*locs[0][i]+locs[1][i]*locs[1][i]);
          std::cout << "  " << i << " " << dist << " " << eulvort[i] << std::endl;
        }
      }
    }

    // find the Lagrangian-computed vorticity on all solution nodes (make a vector of one collection)
    std::vector<Collection> euler_vols;
    euler_vols.emplace_back(solnpts);	// this is a COPY
    //_conv.find_vels(_fs, _vort, _bdry, euler_vols, velandvort, true);
    _conv.find_vort(_vort, _bdry, euler_vols);
    Points<float>& solvedpts = std::get<Points<S>>(euler_vols[0]);
    Vector<S>& lagvort = solvedpts.get_vort();
    assert(lagvort.size() == thisn && "ERROR (Hybrid::step) vorticity from particle sim is not the right size");

    if (true) {
      std::cout << "  vorticity from lagrangian " << std::endl;
      const auto& locs = solvedpts.get_pos();
      for (size_t i=0; i<thisn; ++i) {
        if (std::abs(locs[1][i]/locs[0][i]-0.4868) < 0.01 && locs[0][i] > 0.0) {
          const S dist = std::sqrt(locs[0][i]*locs[0][i]+locs[1][i]*locs[1][i]);
          std::cout << "  " << i << " " << dist << " " << lagvort[i] << std::endl;
        }
      }
    }


    // find the initial vorticity error/deficit
    // subtract the Lagrangian-computed vort from the actual Eulerian vort on those nodes
    // now we have the amount of vorticity we need to re-add to the Lagrangian side
    Vector<S> circ;
    circ.resize(thisn);
    // computes eulvort - lagvort = circ
    std::transform(eulvort.begin(), eulvort.end(),
                   lagvort.begin(),
                   circ.begin(),
                   std::minus<S>());

    // scale by masked cell area
    // uses a mask for the solution nodes (elements here!) to indicate which we will consider,
    // and which are too close to the wall (and thus too thin) to require correction
    // use one full vdelta for this (input _vd)
    (void) coll.set_soln_areas();
    const Vector<S>& area = coll.get_soln_area();
    assert(area.size() == thisn && "ERROR (Hybrid::step) volume area vector is not the right size");

    std::transform(circ.begin(), circ.end(),
                   area.begin(),
                   circ.begin(),
                   std::multiplies<S>());

    // and multiply again by the weight assigned to the grid-to-particle operation
    (void) coll.set_overlap_weights(0.7, 0.2, 1.0, 0.2);
    const Vector<S>& gtop = coll.get_gtop_wgt();

    std::transform(circ.begin(), circ.end(),
                   gtop.begin(),
                   circ.begin(),
                   std::multiplies<S>());

    //size_t num_printed = 0;
    //std::cout << "    indx  area * ( eulvort - lagvort ) = circ" << std::endl;
    //for (size_t i=0; i<thisn and num_printed<5; ++i) {
    //  if (area[i] > 1.e-6) {
    //    std::cout << "    " << i << "  " << area[i] << " * ( " << eulvort[i] << " - " << lagvort[i] << " ) = " << circ[i] << std::endl;
    //    num_printed++;
    //  }
    //}
    //std::cout << "    " << 2013 << "  " << area[2013] << " * ( " << eulvort[2013] << " - " << lagvort[2013] << " ) = " << circ[2013] << std::endl;
    //std::cout << "    " << 1975 << "  " << area[1975] << " * ( " << eulvort[1975] << " - " << lagvort[1975] << " ) = " << circ[1975] << std::endl;
    //std::cout << "    " << 1937 << "  " << area[1937] << " * ( " << eulvort[1937] << " - " << lagvort[1937] << " ) = " << circ[1937] << std::endl;
    if (false) {
      std::cout << "  feedback calculation " << std::endl;
      std::cout << "    indx  rad   area * gtop * ( eulvort - lagvort ) = circ" << std::endl;
      const auto& locs = solvedpts.get_pos();
      for (size_t i=0; i<thisn; ++i) {
        if (std::abs(locs[1][i]/locs[0][i]-0.4868) < 0.01 && locs[0][i] > 0.0) {
          const S dist = std::sqrt(locs[0][i]*locs[0][i]+locs[1][i]*locs[1][i]);
          std::cout << "    " << i << "  " << dist << " " << area[i] << " * " << gtop[i] << " ( " << eulvort[i] << " - " << lagvort[i] << " ) = " << circ[i] << std::endl;
        }
      }
    }


    // find a scaling factor for the error
    double totalcircmag = 0.0;
    for (size_t i=0; i<thisn; ++i) { totalcircmag += std::fabs(eulvort[i]*area[i]*gtop[i]); }
    // if totalcircmag is too small, this may never converge...
    if (totalcircmag < 1.e-8) totalcircmag = 1.0;

    // measure this error
    double maxerror = 0.001;
    double thiserror = 0.0;
    for (size_t i=0; i<thisn; ++i) { thiserror += std::fabs(circ[i]); }
    thiserror /= totalcircmag;
    std::cout << "  initial error " << thiserror << std::endl;

    // iterate toward a solution
    int iter = 0;
    int maxiter = 20;
    while (thiserror>maxerror and iter<maxiter) {

      // resolve it via VRM
      //   for each sub-node of each HO quad, run a VRM onto the existing set of Lagrangian
      //   particles adding where necessary to satisfy at least 0th and 1st moments, if not 2nd
      // resolve it via naive insertion
      //   just create a new particle at the centroid of each element, let merge deal

      // generate particles to fill the deficit
      //   note that cells may be a lot larger than particles!
      ElementPacket<S> newparts = coll.get_equivalent_particles(circ, (S)_vd);

      // make them a new collection in _vort (so that we can delete them later?)
      // or add them to the first collection?
      Points<float>& thevorts = std::get<Points<S>>(_vort[0]);
      thevorts.add_new(newparts, _vd);

      // merge here
      // HACK - hard-coding overlap and merge thresh!
      merge_operation<S>(_vort, 1.5, 0.5, false);

      // find the new Lagrangian vorticity on the solution nodes
      // and the new vorticity error/deficit
      _conv.find_vort(_vort, _bdry, euler_vols);

      // convert the new vort to a new circ
      std::transform(eulvort.begin(), eulvort.end(),
                     lagvort.begin(),
                     circ.begin(),
                     std::minus<S>());

      std::transform(circ.begin(), circ.end(),
                     area.begin(),
                     circ.begin(),
                     std::multiplies<S>());

      std::transform(circ.begin(), circ.end(),
                     gtop.begin(),
                     circ.begin(),
                     std::multiplies<S>());

      //num_printed = 0;
      //std::cout << "    indx  area * ( eulvort - lagvort ) = circ" << std::endl;
      //for (size_t i=0; i<thisn and num_printed<5; ++i) {
      //  if (area[i] > 1.e-6) {
      //    std::cout << "    " << i << "  " << area[i] << " * ( " << eulvort[i] << " - " << lagvort[i] << " ) = " << circ[i] << std::endl;
      //    num_printed++;
      //  }
      //}
      //std::cout << "    " << 2013 << "  " << area[2013] << " * ( " << eulvort[2013] << " - " << lagvort[2013] << " ) = " << circ[2013] << std::endl;
      //std::cout << "    " << 1975 << "  " << area[1975] << " * ( " << eulvort[1975] << " - " << lagvort[1975] << " ) = " << circ[1975] << std::endl;
      //std::cout << "    " << 1937 << "  " << area[1937] << " * ( " << eulvort[1937] << " - " << lagvort[1937] << " ) = " << circ[1937] << std::endl;
    if (true) {
      std::cout << "  feedback calculation " << std::endl;
      std::cout << "    indx  rad   area * gtop * ( eulvort - lagvort ) = circ" << std::endl;
      const auto& locs = solvedpts.get_pos();
      for (size_t i=0; i<thisn; ++i) {
        if (std::abs(locs[1][i]/locs[0][i]-0.4868) < 0.01 && locs[0][i] > 0.0) {
          const S dist = std::sqrt(locs[0][i]*locs[0][i]+locs[1][i]*locs[1][i]);
          std::cout << "    " << i << "  " << dist << " " << area[i] << " * " << gtop[i] << " * ( " << eulvort[i] << " - " << lagvort[i] << " ) = " << circ[i] << std::endl;
        }
      }
    }

      // measure this error
      thiserror = 0.0;
      for (size_t i=0; i<thisn; ++i) { thiserror += std::fabs(circ[i]); }
      thiserror /= totalcircmag;
      std::cout << "  iter " << (iter+1) << " has error " << thiserror << std::endl;

      iter++;
    }

    //
    // part D - update HO vorticity from particle solution
    //
    // Daeninck calls this "Euler adjustment region"
    // we use a set of weights to adjust the difference

    std::cout << "Inside Hybrid::step updating particle strengths" << std::endl;

    // convert the new vort to a new circ
    std::transform(eulvort.begin(), eulvort.end(),
                   lagvort.begin(),
                   circ.begin(),
                   std::minus<S>());

    // now circulation contains the vorticity difference between the euler and lagrangian solutions

    // scale it by the particle-to-grid weights
    const Vector<S>& ptog = coll.get_ptog_wgt();
    std::transform(circ.begin(), circ.end(),
                   ptog.begin(),
                   circ.begin(),
                   std::multiplies<S>());

    // and add it to the HO vorticities
    std::transform(eulvort.begin(), eulvort.end(),
                   circ.begin(),
                   eulvort.begin(),
                   std::minus<S>());

    if (true) {
      std::cout << "  feedback calculation " << std::endl;
      std::cout << "    indx  rad   ptog   neweulvort   lagvort" << std::endl;
      const auto& locs = solvedpts.get_pos();
      for (size_t i=0; i<thisn; ++i) {
        if (std::abs(locs[1][i]/locs[0][i]-0.4868) < 0.01 && locs[0][i] > 0.0) {
          const S dist = std::sqrt(locs[0][i]*locs[0][i]+locs[1][i]*locs[1][i]);
          std::cout << "    " << i << "  " << dist << "  " << ptog[i] << "  " << eulvort[i] << "  " << lagvort[i] << std::endl;
        }
      }
    }


    // finally, replace the vorticity in the HO solver with these new values
#ifdef HOFORTRAN
    (void) setsolnvort_d((int32_t)eulvort.size(), eulvort.data());
#else
    (void) solver.setsolnvort_d(eulvort);
#endif

  }

  // done!!!
}

//
// read/write parameters to json
//

// read "simparams" json object
template <class S, class A, class I>
void Hybrid<S,A,I>::from_json(const nlohmann::json simj) {
  if (simj.find("hybrid") != simj.end()) {
    nlohmann::json j = simj["hybrid"];

    active = j.value("enabled", false);
    elementOrder = j.value("elementOrder", 1);
    timeOrder = j.value("timeOrder", 1);
    numSubsteps = j.value("numSubsteps", 100);
    preconditioner = j.value("preconditioner", "none");
    solverType = j.value("solverType", "fgmres");
  }
}

// create and write a json object for all diffusion parameters
template <class S, class A, class I>
void Hybrid<S,A,I>::add_to_json(nlohmann::json& simj) const {
  nlohmann::json j;
  j["enabled"] = active;
  j["elementOrder"] = elementOrder; //1-5
  j["timeOrder"] = timeOrder; //1,2,4
  j["numSubsteps"] = numSubsteps; //1-1000
  j["preconditioner"] = preconditioner;
  j["solverType"] = solverType;

  simj["hybrid"] = j;
}

#ifdef USE_IMGUI
//
// draw advanced options parts of the GUI
//
template <class S, class A, class I>
void Hybrid<S,A,I>::draw_advanced() {

  ImGui::Separator();
  ImGui::Spacing();
  ImGui::Text("Hybrid/Grid settings");

  ImGui::Checkbox("Enabled", &active);

  if (active) {

  ImGui::SliderInt("Element Order", &elementOrder, 1, 5);
  ImGui::SliderInt("Time Order", &timeOrder, 1, 4);

  /*
  const int numTimeOrders = 3;
  static int timeI = 0;
  const char* timeOrders[] = {"1", "2", "4"};
  ImGui::Combo("Select Time Order", &timeI, timeOrders, numTimeOrders);
  switch (timeI) {
    case 0: timeOrder = 1; break;
    case 1: timeOrder = 2; break;
    case 2: timeOrder = 4; break;
  }
  */

  if (ImGui::InputInt("Number of Substeps", &numSubsteps)) {
    if (numSubsteps < 1) { numSubsteps = 1; }
    else if (numSubsteps > 1000) { numSubsteps = 1000; }
  }
  
  const int numPreconditioners = 1;
  static int preconI = 0;
  const char* preconditioners[] = {"none"};
  ImGui::Combo("Select Preconditioner", &preconI, preconditioners, numPreconditioners);
  preconditioner = preconditioners[preconI];
  
  const int numSolvers = 1;
  static int solverI = 0;
  const char* solvers[] = {"fgmres"};
  ImGui::Combo("Select Solver", &solverI, solvers, numSolvers);
  solverType = solvers[solverI];
  }
}
#endif
