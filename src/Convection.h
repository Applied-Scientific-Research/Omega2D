/*
 * Convection.h - a class for forward integration of elements and their strengths
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Collection.h"
#include "CollectionHelper.h"
#include "Influence.h"
#include "BEM.h"
#include "BEMHelper.h"

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
  void find_vels( const std::array<double,Dimensions>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&);
  void advect_1st(const double,
                  const double,
                  const std::array<double,Dimensions>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  BEM<S,I>&);
  void advect_2nd(const double,
                  const double,
                  const std::array<double,Dimensions>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  BEM<S,I>&);

private:
  // local copies of particle data
  //Particles<S> temp;
};


//
// helper function to find velocities at a given state, assuming BEM is solved
//
template <class S, class A, class I>
void Convection<S,A,I>::find_vels(const std::array<double,Dimensions>& _fs,
                                  std::vector<Collection>&             _vort,
                                  std::vector<Collection>&             _bdry,
                                  std::vector<Collection>&             _targets) {

  //if (_targets.size() > 0) std::cout << std::endl << "Solving for velocities" << std::endl;
  //if (_targets.size() > 0) std::cout << std::endl;

  // need this for dispatching velocity influence calls, template param is accumulator type
  // should the solution_t be an argument to the constructor?
  InfluenceVisitor<A> visitor;

  // add vortex and source strengths to account for rotating bodies
  for (auto &src : _bdry) {
    std::visit([=](auto& elem) { elem.add_solved_rot_strengths(1.0); }, src);
  }

  // find the influence on every field point/tracer element
  for (auto &targ : _targets) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl;

    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);

    // accumulate from vorticity
    for (auto &src : _vort) {
      std::visit(visitor, src, targ);
    }

    // accumulate from boundaries
    for (auto &src : _bdry) {
      // call the Influence routine for these collections
      std::visit(visitor, src, targ);
    }

    // add freestream and divide by 2pi
    std::visit([=](auto& elem) { elem.finalize_vels(_fs); }, targ);
  }

  // remove vortex and source strengths due to rotation
  for (auto &src : _bdry) {
    std::visit([=](auto& elem) { elem.add_solved_rot_strengths(-1.0); }, src);
  }
}

//
// first-order Euler forward integration
//
template <class S, class A, class I>
void Convection<S,A,I>::advect_1st(const double _time,
                                   const double _dt,
                                   const std::array<double,Dimensions>& _fs,
                                   std::vector<Collection>&             _vort,
                                   std::vector<Collection>&             _bdry,
                                   std::vector<Collection>&             _fldpt,
                                   BEM<S,I>&                            _bem) {

  std::cout << "Inside advect_1st with dt=" << _dt << std::endl;

  // part A - unknowns

  solve_bem<S,A,I>(_time, _fs, _vort, _bdry, _bem);

  // part B - knowns

  find_vels(_fs, _vort, _bdry, _vort);
  find_vels(_fs, _vort, _bdry, _fldpt);

  // part C - convection here

  std::cout << std::endl << "Convection step" << std::endl;

  // move every movable element
  for (auto &coll : _vort) {
    std::visit([=](auto& elem) { elem.move(_time, _dt); }, coll);
  }
  for (auto &coll : _bdry) {
    std::visit([=](auto& elem) { elem.move(_time, _dt); }, coll);
  }
  for (auto &coll : _fldpt) {
    std::visit([=](auto& elem) { elem.move(_time, _dt); }, coll);
  }
}


//
// second-order RK2 forward integration
//
template <class S, class A, class I>
void Convection<S,A,I>::advect_2nd(const double _time,
                                   const double _dt,
                                   const std::array<double,Dimensions>& _fs,
                                   std::vector<Collection>&             _vort,
                                   std::vector<Collection>&             _bdry,
                                   std::vector<Collection>&             _fldpt,
                                   BEM<S,I>&                            _bem) {

  std::cout << "Inside advect_2nd with dt=" << _dt << std::endl;

  // take the first Euler step ---------

  // perform the first BEM
  solve_bem<S,A,I>(_time, _fs, _vort, _bdry, _bem);

  // find the derivatives
  find_vels(_fs, _vort, _bdry, _vort);
  find_vels(_fs, _vort, _bdry, _fldpt);

  // advect into an intermediate system
  std::vector<Collection> interim_vort = _vort;
  for (auto &coll : interim_vort) {
    std::visit([=](auto& elem) { elem.move(_time, _dt); }, coll);
  }
  // now _vort has its original positions and the velocities evaluated there
  // and interm_vort has the positions at t+dt

  // do the same for fldpt
  std::vector<Collection> interim_fldpt = _fldpt;
  for (auto &coll : interim_fldpt) {
    std::visit([=](auto& elem) { elem.move(_time, _dt); }, coll);
  }

  // begin the 2nd step ---------

  // perform the second BEM
  //solve_bem<S,A,I>(_time + _dt, _fs, interim_vort, interim_bdry, _bem);
  solve_bem<S,A,I>(_time + _dt, _fs, interim_vort, _bdry, _bem);

  // find the derivatives
  //find_vels(_fs, interim_vort, interim_bdry, interim_fldpt);
  find_vels(_fs, interim_vort, _bdry, interim_vort);
  find_vels(_fs, interim_vort, _bdry, interim_fldpt);

  // _vort still has its original positions and the velocities evaluated there
  // but interm_vort now has the velocities at t+dt

  // advect using the combination of both velocities
  auto v1p = _vort.begin();
  auto v2p = interim_vort.begin();
  for (size_t i = 0; i < _vort.size(); ++i) {
    Collection& c1 = *v1p;
    Collection& c2 = *v2p;
    // switch based on what type is actually held in the std::variant
    if (std::holds_alternative<Points<float>>(c1) and std::holds_alternative<Points<float>>(c2)) {
      Points<float>& p1 = std::get<Points<float>>(c1);
      Points<float>& p2 = std::get<Points<float>>(c2);
      p1.move(_time, _dt, 0.5, p1, 0.5, p2);
    }
    ++v1p;
    ++v2p;
  }

  v1p = _fldpt.begin();
  v2p = interim_fldpt.begin();
  for (size_t i = 0; i < _fldpt.size(); ++i) {
    Collection& c1 = *v1p;
    Collection& c2 = *v2p;
    // switch based on what type is actually held in the std::variant
    if (std::holds_alternative<Points<float>>(c1) and std::holds_alternative<Points<float>>(c2)) {
      Points<float>& p1 = std::get<Points<float>>(c1);
      Points<float>& p2 = std::get<Points<float>>(c2);
      p1.move(_time, _dt, 0.5, p1, 0.5, p2);
    }
    ++v1p;
    ++v2p;
  }
}

