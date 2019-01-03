/*
 * Convection.h - a class for forward integration of elements and their strengths
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Collection.h"
#include "CollectionHelper.h"
#include "NewInfluence.h"
#include "Boundaries.h"
#include "Influence.h"
#include "Vorticity.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <variant>


//
// One step of convection of elements, allowing for boundary conditions
//
// templatized on 'S'torage, 'A'ccumulator types, and 'I'ndex types
//
template <class S, class A, class I>
class Convection {
public:
  Convection() {}
  void advect_1st(const double,
                  const std::array<float,2>&,
                  Vorticity<S,I>&,
                  Boundaries<S,I>&);
  void advect_1st(const double,
                  const std::array<float,2>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&);
  void advect_2nd(const double,
                  const std::array<float,2>&,
                  Vorticity<S,I>&,
                  Boundaries<S,I>&);

private:
  // local copies of particle data
  //Particles<S> temp;
};


//
// first-order Euler forward integration
//
template <class S, class A, class I>
void Convection<S,A,I>::advect_1st(const double _dt,
                                   const std::array<float,2>& _fs,
                                   Vorticity<S,I>& _vort,
                                   Boundaries<S,I>& _bdry) {
  //std::cout << "  inside advect_1st with dt=" << _dt << std::endl;

  // part A - unknowns

  // when bodies are passed in, we must solve for their strengths first
  // eventually, this will do it all:
  //_bdry.find_strengths(_fs, _vort);
  // but for now, try this
  if (_bdry.exists()) {
    _bdry.reset_vels();
    add_influence<S,A,I>(_vort, _bdry);
    _bdry.scale_and_add_freestream(_fs);
    _bdry.find_strengths();
  }

  // part B - knowns

  // zero out the velocity vector - only do this once, here
  _vort.reset_vels();

  // add the influence of the body
  if (_bdry.exists()) add_influence<S,A,I>(_bdry, _vort);

  // add velocities from the other particles
  add_influence<S,A,I>(_vort, _vort);

  // no change in strength or radius for now

  // scale and add freestream
  _vort.scale_and_add_freestream(_fs);

  // advect in place
  _vort.step_in_place((S)_dt);

  //std::cout << "After 1st order convection, particles are:" << std::endl;
  //if (n>0) std::cout << "  part 0 with str " << x[2] << " is at " << x[0] << " " << x[1] << std::endl;
}


template <class S, class A, class I>
void Convection<S,A,I>::advect_1st(const double _dt,
                                   const std::array<float,Dimensions>& _fs,
                                   std::vector<Collection>& _vort,
                                   std::vector<Collection>& _bdry,
                                   std::vector<Collection>& _fldpt) {

  std::cout << "  inside advect_1st with dt=" << _dt << std::endl;

  // part A - unknowns

  //solve_bem(_fs, _vort, _bdry);

  // part B - knowns

  //find_vels(_fs, _vort, _bdry, _fldpt);

  // need this for dispatching velocity influence calls, template param is accumulator type
  // should the solution_t be an argument to the constructor?
  InfluenceVisitor<A> visitor;

  // find the influence on every vorticity element
  for (auto &targ : _vort) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl << std::flush;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src : _vort) {
      // how do I specify the solver?
      std::visit(visitor, src, targ);
    }
    // accumulate from boundaries
    for (auto &src : _bdry) {
      std::visit(visitor, src, targ);
    }
    // add freestream and divide by 2pi
    std::visit([=](auto& elem) { elem.finalize_vels(_fs); }, targ);
  }

  // part C - convection here

  std::cout << std::endl << "Convection step" << std::endl;

  // move every movable element
  for (auto &coll : _vort) {
    std::visit([=](auto& elem) { elem.move(_dt); }, coll);
  }
  for (auto &coll : _bdry) {
    std::visit([=](auto& elem) { elem.move(_dt); }, coll);
  }
  for (auto &coll : _fldpt) {
    std::visit([=](auto& elem) { elem.move(_dt); }, coll);
  }
}


//
// second-order RK2 forward integration
//
template <class S, class A, class I>
void Convection<S,A,I>::advect_2nd(const double _dt,
                                   const std::array<float,2>& _fs,
                                   Vorticity<S,I>& _vort,
                                   Boundaries<S,I>& _bdry) {
  //std::cout << "  inside advect_2nd with dt=" << dt << std::endl;

  // take the first Euler step

  // perform the first BEM
  if (_bdry.exists()) {
    _bdry.reset_vels();
    // move boundaries to time t
    add_influence<S,A,I>(_vort, _bdry);
    _bdry.scale_and_add_freestream(_fs);
    _bdry.find_strengths();
  }

  // find the derivatives
  _vort.reset_vels();
  if (_bdry.exists()) add_influence<S,A,I>(_bdry, _vort);
  add_influence<S,A,I>(_vort, _vort);
  _vort.scale_and_add_freestream(_fs);

  // advect into an intermediate system
  Vorticity<S,I> interm_vort = _vort;	// uses copy constructor, not assignment operator
  interm_vort.step_in_place((S)_dt);

  // now _vort has its original positions and the velocities evaluated there
  // and interm_vort has the positions at t+dt

  // begin the 2nd step

  // perform the second BEM
  if (_bdry.exists()) {
    _bdry.reset_vels();
    // move boundaries to t+dt
    add_influence<S,A,I>(interm_vort, _bdry);
    _bdry.scale_and_add_freestream(_fs);
    _bdry.find_strengths();
  }

  // find the derivatives
  interm_vort.reset_vels();
  if (_bdry.exists()) add_influence<S,A,I>(_bdry, interm_vort);
  add_influence<S,A,I>(interm_vort, interm_vort);
  interm_vort.scale_and_add_freestream(_fs);

  // _vort still has its original positions and the velocities evaluated there
  // but interm_vort now has the velocities at t+dt

  // advect using the combination of both velocities
  _vort.step_in_place((S)_dt, interm_vort);


  //std::cout << "After 1st order convection, particles are:" << std::endl;
  //for (size_t i=0; i<4*n; i+=4) {
  //  std::cout << "  " << i/4 << "   " << u[i] << " " << u[i+1] << "   " << x[i] << " " << x[i+1] << std::endl;
  //}
  //if (n>0) std::cout << "  part 0 with str " << x[2] << " is at " << x[0] << " " << x[1] << std::endl;

  //if (n>0) std::cout << "  part 0 with str " << x[2] << " is at " << x[0] << " " << x[1] << std::endl;
}

