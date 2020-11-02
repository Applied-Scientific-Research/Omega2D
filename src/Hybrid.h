/*
 * Hybrid.h - coordinate with an external Eulerian flow solver to compute near-body flow
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"

#include <iostream>
#include <vector>

/*
// these could be in headers exposed by the external solver
extern "C" float external_euler_init_f_(int*, const float*, const float*, const float*, const float*,
                                        int*, const float*, const float*, float*, float*);
extern "C" float external_euler_init_d_(int*, const double*, const double*, const double*, const double*,
                                        int*, const double*, const double*, double*, double*);
*/


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
      initialized(false)//,
      //vrm(),
      //h_nu(0.1)
    {}

  const bool is_active() const { return active; }
  void set_active(const bool _do_hybrid) { active = _do_hybrid; }
  void activate() { active = true; }
  void deactivate() { active = false; }

  void init( std::vector<Collection>&);
  void reset();
  void step( const double,
             const double,
             const std::array<double,Dimensions>&,
             std::vector<Collection>&,
             std::vector<Collection>&,
             std::vector<Collection>&,
             std::vector<Collection>&);

private:
  // are we even using the hybrid scheme?
  bool active;
  bool initialized;

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
void Hybrid<S,A,I>::init(std::vector<Collection>& _euler) {
  std::cout << "Inside Hybrid::init with " << _euler.size() << " volumes" << std::endl;

  // call the external solver with the current geometry
  // this will calculate the Jacobian and other cell-specific properties
  //flops = external_euler_init_f_(&ns, sx[0].data(), sx[1].data(),    ss.data(),    sr.data(),
  //                               &nt, tx[0].data(), tx[1].data(), tu[0].data(), tu[1].data());

  initialized = true;
}


//
// Simulation is reset - do we need to clear the Euler grid initialization?
//
template <class S, class A, class I>
void Hybrid<S,A,I>::reset() {
  initialized = false;
}


//
// Forward integration step
//
template <class S, class A, class I>
void Hybrid<S,A,I>::step(const double                         _time,
                         const double                         _dt,
                         const std::array<double,Dimensions>& _fs,
                         std::vector<Collection>&             _vort,
                         std::vector<Collection>&             _bdry,
                         std::vector<Collection>&             _fldpt,
                         std::vector<Collection>&             _euler) {

  if (not active) return;

  if (not initialized) init(_euler);

  std::cout << "Inside Hybrid::step at t=" << _time << " and dt=" << _dt << std::endl;

  // part A - prepare BCs for Euler solver

  // transform all elements to their appropriate locations
  // isolate open/outer boundaries
  //Points<> euler_bdry = _euler.get_bc_nodes(_time);
  // find vels there
  //Convection::find_vels<S,A,I>(_fs, _vort, _bdry, euler_bdry);
  // convert to transferable packet
  //find_vels<S,A,I>(_time, _fs, _vort, _bdry, _bem);  ???

  // part B - call Euler solver

  // transfer BC packet and call solver
  //flops = external_euler_solve_f_(&ns, sx[0].data(), sx[1].data(),    ss.data(),    sr.data(),
  //                                &nt, tx[0].data(), tx[1].data(), tu[0].data(), tu[1].data());

  // pull results from external solver
  //external_euler_vorts_f_(stuff);

  // convert vorticity results to drawable and writable Volume elements

  // part C - update particle strengths accordionly

  // "update" strengths on vortex particles within the boundary

  // here's one way to do it:
  // identify all free vortex particles inside of euler regions and remove them
  //   (add up how much circulation we remove)
  // re-run BEM and compute vorticity on all HO volume nodes
  // subtract the Lagrangian-computed vort from the actual Eulerian vort on those nodes
  // now we have the amount of vorticity we need to re-add to the Lagrangian side
  // for each sub-node of each HO quad, run a VRM onto the existing set of Lagrangian
  //   particles adding where necessary to satisfy at least 0th and 1st moments, if not 2nd

}

