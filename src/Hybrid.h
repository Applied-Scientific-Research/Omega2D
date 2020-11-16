/*
 * Hybrid.h - coordinate with an external Eulerian flow solver to compute near-body flow
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Collection.h"
#include "BEM.h"

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
      initialized(false),
      elementOrder(1),
      timeOrder(1),
      numSubsteps(100),
      preconditioner("none"),
      solverType("fgmres")//,
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
             BEM<S,I>&,
             std::vector<Collection>&);

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
  uint8_t timeOrder;
  int numSubsteps;
  std::string preconditioner;
  std::string solverType;

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
                         BEM<S,I>&                            _bem,
                         std::vector<Collection>&             _euler) {

  if (not active) return;

  if (not initialized) init(_euler);

  std::cout << "Inside Hybrid::step at t=" << _time << " and dt=" << _dt << std::endl;

  // part A - prepare BCs for Euler solver

  // transform all elements to their appropriate locations
  for (auto &coll : _euler) {
    std::visit([=](auto& elem) { elem.move(_time, _dt); }, coll);
  }
  // update the BEM solution
  solve_bem<S,A,I>(_time, _fs, _vort, _bdry, _bem);
  // get vels on each euler region
  for (auto &coll : _euler) {
    // isolate open/outer boundaries
    //Points<A> euler_bdry = std::visit([=](auto& elem) { elem.get_bc_nodes(_time); }, coll);
    //Points<A> euler_bdry = coll.get_bc_nodes(_time);

    // find vels there
    //Convection::find_vels<A,A,I>(_fs, _vort, _bdry, euler_bdry);

    // convert to transferable packet
    //std::array<Vector<A>,Dimensions> openvels = euler_bdry.get_vel();

    // transfer BC packet to solver
    //(void) hosolver_setopenvels_d_(coll.get_n(), openvels[0], openvels[1]);
  }

  // part B - call Euler solver

  // call solver
  //(void) hosolver_solveto_d_(_time);

  // pull results from external solver
  //std::vector<A> allvorts;
  //(void) hosolver_getallvorts_d_(np, allvorts);

  // convert vorticity results to drawable and writable Volume elements

  // part C - update particle strengths accordionly

  // "update" strengths on vortex particles within the boundary

  // here's one way to do it:
  // identify all free vortex particles inside of euler regions and remove them
  //   (add up how much circulation we remove) - or not?

  // re-run BEM and compute vorticity on all HO volume nodes - bem already run
  //Points<S> euler_vol = coll.get_vol_nodes(_time);
  //Convection::find_vels<S,A,I>(_fs, _vort, _bdry, euler_vol, velandvort);

  // subtract the Lagrangian-computed vort from the actual Eulerian vort on those nodes
  // now we have the amount of vorticity we need to re-add to the Lagrangian side
  // for each sub-node of each HO quad, run a VRM onto the existing set of Lagrangian
  //   particles adding where necessary to satisfy at least 0th and 1st moments, if not 2nd

  // simple way: just create a new particle at the centroid of each element, let merge deal

}

//
// read/write parameters to json
//

// read "simparams" json object
template <class S, class A, class I>
void Hybrid<S,A,I>::from_json(const nlohmann::json j) {
  active = j.value("enabled", false);
  elementOrder = j.value("elementOrder", 1);
  timeOrder = j.value("timeOrder", 1);
  numSubsteps = j.value("numSubsteps", 100);
  preconditioner = j.value("preconditioner", "none");
  solverType = j.value("solverType", "fgmres");
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
  if (ImGui::InputInt("Element Order", &elementOrder)) {
    if (elementOrder < 1) { elementOrder = 1; }
    else if (elementOrder > 5) { elementOrder = 5; }
  }

  const int numTimeOrders = 3;
  static int timeI = 0;
  const char* timeOrders[] = {"1", "2", "4"};
  ImGui::Combo("Select Time Order", &timeI, timeOrders, numTimeOrders);
  switch (timeI) {
    case 0: timeOrder = 1; break;
    case 1: timeOrder = 2; break;
    case 2: timeOrder = 4; break;
  }

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
#endif
