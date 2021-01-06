/*
 * Convection.h - a class for forward integration of elements and their strengths
 *
 * (c)2017-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Collection.h"
#include "CollectionHelper.h"
#include "Influence.h"
#include "BEM.h"
#include "BEMHelper.h"
#include "InfluenceVort.h"
#include "Reflect.h"
#include "GuiHelper.h"

#include <json/json.hpp>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <variant>
#include <cassert>


//
// One step of convection of elements, allowing for boundary conditions
//
// templatized on 'S'torage, 'A'ccumulator types, and 'I'ndex types
//
template <class S, class A, class I>
class Convection {
public:
  Convection()
    : convection_order(2),
      conv_env()
    {}

  void find_vort( std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&);
  void find_vels( const std::array<double,Dimensions>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  const results_t _results = velonly,
                  const bool _force = false);
  void advect(    const double,
                  const double,
                  const std::array<double,Dimensions>&,
                  const S,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  BEM<S,I>&);
  void advect_1st(const double,
                  const double,
                  const std::array<double,Dimensions>&,
                  const S,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  BEM<S,I>&);
  void advect_2nd(const double,
                  const double,
                  const std::array<double,Dimensions>&,
                  const S,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  std::vector<Collection>&,
                  BEM<S,I>&);

#ifdef USE_IMGUI
  void draw_advanced();
#endif

  // read/write parameters
  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

private:
  // local copies of particle data
  //Particles<S> temp;

  // integrator order
  int32_t convection_order;

  // execution environment for velocity summations (not BEM)
  ExecEnv conv_env;
};


//
// helper function to find vorticity at a given state
//
template <class S, class A, class I>
void Convection<S,A,I>::find_vort(std::vector<Collection>&             _vort,
                                  std::vector<Collection>&             _bdry,
                                  std::vector<Collection>&             _targets) {

  //if (_targets.size() > 0) std::cout << std::endl << "Solving for velocities" << std::endl;
  //if (_targets.size() > 0) std::cout << std::endl;

  // find the influence on every target collection that is Points
  for (auto &targ : _targets) if (std::holds_alternative<Points<S>>(targ)) {
    Points<S>& ptarg = std::get<Points<S>>(targ);

    std::cout << "  Solving vorticity on" << to_string(targ) << std::endl;

    // zero vels and vorticity
    ptarg.zero_vels();

    // accumulate from vorticity, but only from Points
    for (auto &src : _vort) if (std::holds_alternative<Points<S>>(src)) {
      Points<S>& psrc = std::get<Points<S>>(src);
      points_affect_points_vorticity<float,A>(psrc, ptarg, conv_env);
    }

    // accumulate from reactive, but only from Points
    for (auto &src : _bdry) if (std::holds_alternative<Points<S>>(src)) {
      Points<S>& psrc = std::get<Points<S>>(src);
      points_affect_points_vorticity<float,A>(psrc, ptarg, conv_env);
    }

    // finalize vels and vorticity by dividing by constant
    ptarg.finalize_vels(std::array<double,Dimensions>({0.0,0.0}));
  }
}


//
// helper function to find velocities at a given state, assuming BEM is solved
//
template <class S, class A, class I>
void Convection<S,A,I>::find_vels(const std::array<double,Dimensions>& _fs,
                                  std::vector<Collection>&             _vort,
                                  std::vector<Collection>&             _bdry,
                                  std::vector<Collection>&             _targets,
                                  const results_t                      _results,
                                  const bool                           _force) {

  //if (_targets.size() > 0) std::cout << std::endl << "Solving for velocities" << std::endl;
  //if (_targets.size() > 0) std::cout << std::endl;

  // need this for dispatching velocity influence calls, template param is accumulator type
  // member variable is passed-in execution environment
  InfluenceVisitor<A> visitor = {ResultsType(_results), conv_env};

  // add vortex and source strengths to account for rotating bodies
  for (auto &src : _bdry) {
    std::visit([=](auto& elem) { elem.add_solved_rot_strengths(1.0); }, src);
  }

  // find the influence on every field point/tracer element
  for (auto &targ : _targets) {

    // only do this for Lagrangian items, or if told to (as when exporting to vtk)
    const move_t tmt = std::visit([=](auto& elem) { return elem.get_movet(); }, targ);
    if (not (_force or tmt == lagrangian)) continue;

    std::cout << "  Solving" << ResultsType(_results).to_string() << " on" << to_string(targ) << std::endl;

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

    // add freestream and divide by constant
    std::visit([=](auto& elem) { elem.finalize_vels(_fs); }, targ);
  }

  // remove vortex and source strengths due to rotation
  for (auto &src : _bdry) {
    std::visit([=](auto& elem) { elem.add_solved_rot_strengths(-1.0); }, src);
  }
}

//
// general forward integrator
//
template <class S, class A, class I>
void Convection<S,A,I>::advect(const double                         _time,
                               const double                         _dt,
                               const std::array<double,Dimensions>& _fs,
                               const S                              _ips,
                               std::vector<Collection>&             _vort,
                               std::vector<Collection>&             _bdry,
                               std::vector<Collection>&             _fldpt,
                               BEM<S,I>&                            _bem) {

  assert(convection_order > 0 and convection_order < 3 && "Convection integrator orders over 2 unsupported");

  if (convection_order == 1) advect_1st(_time, _dt, _fs, _ips, _vort, _bdry, _fldpt, _bem);
  else advect_2nd(_time, _dt, _fs, _ips, _vort, _bdry, _fldpt, _bem);
}


//
// first-order Euler forward integration
//
template <class S, class A, class I>
void Convection<S,A,I>::advect_1st(const double                         _time,
                                   const double                         _dt,
                                   const std::array<double,Dimensions>& _fs,
                                   const S                              _ips,
                                   std::vector<Collection>&             _vort,
                                   std::vector<Collection>&             _bdry,
                                   std::vector<Collection>&             _fldpt,
                                   BEM<S,I>&                            _bem) {

  std::cout << "Inside Convection::advect_1st with dt=" << _dt << std::endl;

  // part A - unknowns

  // push away particles inside or too close to the body
  assert(M_PI != 0); // Can't divide by 0
  clear_inner_layer<S>(1, _bdry, _vort, 0.5/std::sqrt(2.0*M_PI), _ips);
  // and solve the bem
  solve_bem<S,A,I>(_time, _fs, _vort, _bdry, _bem);

  // part B - knowns

  find_vels(_fs, _vort, _bdry, _vort);
  find_vels(_fs, _vort, _bdry, _fldpt);

  // part C - convection here

  //std::cout << std::endl << "Convection step" << std::endl;

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
void Convection<S,A,I>::advect_2nd(const double                         _time,
                                   const double                         _dt,
                                   const std::array<double,Dimensions>& _fs,
                                   const S                              _ips,
                                   std::vector<Collection>&             _vort,
                                   std::vector<Collection>&             _bdry,
                                   std::vector<Collection>&             _fldpt,
                                   BEM<S,I>&                            _bem) {

  std::cout << "Inside Convection::advect_2nd with dt=" << _dt << std::endl;

  // take the first Euler step ---------

  // push away particles inside or too close to the body
  assert(M_PI != 0); // Can't divide by 0
  clear_inner_layer<S>(1, _bdry, _vort, 0.5/std::sqrt(2.0*M_PI), _ips);
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

  // push away particles inside or too close to the body
  clear_inner_layer<S>(1, _bdry, interim_vort, 0.5/std::sqrt(2.0*M_PI), _ips);
  // perform the second BEM
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


#ifdef USE_IMGUI
//
// draw advanced options parts of the GUI
//
template <class S, class A, class I>
void Convection<S,A,I>::draw_advanced() {

  //ImGui::Separator();
  ImGui::Spacing();
  ImGui::Text("Convection settings");

  static bool use_internal_solver = conv_env.is_internal();
#ifdef EXTERNAL_VEL_SOLVE
  ImGui::Checkbox("Use internal velocity solver", &use_internal_solver);
  conv_env.set_internal(use_internal_solver);
  ImGui::SameLine();
  ShowHelpMarker("Use the internal method to calculate velocities. Uncheck to use an external solver.");
#endif

  if (use_internal_solver) {
#ifdef USE_VC
  #ifdef USE_OGL_COMPUTE
    static int acc_item = 2;
    const char* acc_items[] = { "x86 (CPU)", "Vc SIMD (CPU)", "OpenGL (GPU)" };
    ImGui::PushItemWidth(240);
    ImGui::Combo("Select instructions", &acc_item, acc_items, 3);
    ImGui::PopItemWidth();
    switch(acc_item) {
        case 0: conv_env.set_instrs(cpu_x86); break;
        case 1: conv_env.set_instrs(cpu_vc); break;
        case 2: conv_env.set_instrs(gpu_opengl); break;
    } // end switch
  #else
    static int acc_item = 1;
    const char* acc_items[] = { "x86 (CPU)", "Vc SIMD (CPU)" };
    ImGui::PushItemWidth(240);
    ImGui::Combo("Select instructions", &acc_item, acc_items, 2);
    ImGui::PopItemWidth();
    switch(acc_item) {
        case 0: conv_env.set_instrs(cpu_x86); break;
        case 1: conv_env.set_instrs(cpu_vc); break;
    } // end switch
  #endif
#else
  #ifdef USE_OGL_COMPUTE
    static int acc_item = 1;
    const char* acc_items[] = { "x86 (CPU)", "OpenGL (GPU)" };
    ImGui::PushItemWidth(240);
    ImGui::Combo("Select instructions", &acc_item, acc_items, 2);
    ImGui::PopItemWidth();
    switch(acc_item) {
        case 0: conv_env.set_instrs(cpu_x86); break;
        case 1: conv_env.set_instrs(gpu_opengl); break;
    } // end switch
  #else
    // none!
    ImGui::Text("Instructions are x86 (CPU)");
  #endif
#endif

    // now, depending on which was selected, allow different summation algorithms
    //static int algo_item = 0;
    const accel_t accel_selected = conv_env.get_instrs();
    if (accel_selected == cpu_x86) {
      ImGui::Text("Algorithm is direct, O(N^2)");
      conv_env.set_summation(direct);
      //const char* algo_items[] = { "direct, O(N^2)", "treecode, O(NlogN)" };
      //ImGui::PushItemWidth(240);
      //ImGui::Combo("Select algorithm", &algo_item, algo_items, 2);
      //ImGui::PopItemWidth();
      //switch(algo_item) {
      //  case 0: conv_env.set_summation(direct); break;
      //  case 1: conv_env.set_summation(barneshut); break;
      //} // end switch
    } else {
      ImGui::Text("Algorithm is direct, O(N^2)");
      conv_env.set_summation(direct);
    }
  }
}
#endif

//
// read/write parameters to json
//

// read "simparams" json object
template <class S, class A, class I>
void Convection<S,A,I>::from_json(const nlohmann::json j) {

  if (j.find("timeOrder") != j.end()) {
    convection_order = j["timeOrder"];
    std::cout << "  setting forward integrator order= " << convection_order << std::endl;
  }

  // may eventually set method by name?
}

// create and write a json object for all diffusion parameters
template <class S, class A, class I>
void Convection<S,A,I>::add_to_json(nlohmann::json& j) const {
  j["timeOrder"] = convection_order;
}

