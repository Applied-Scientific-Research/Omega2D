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
#include "dummysolver.h"

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
      solverType("fgmres"),
      solver()
      //vrm(),
      //h_nu(0.1)
    {}

  const bool is_active() const { return active; }
  void set_active(const bool _do_hybrid) { active = _do_hybrid; }
  void activate() { active = true; }
  void deactivate() { active = false; }

  void init( std::vector<HOVolumes<S>>&);
  void reset();
  void step( const double,
             const double,
             const std::array<double,Dimensions>&,
             std::vector<Collection>&,
             std::vector<Collection>&,
             BEM<S,I>&,
             std::vector<HOVolumes<S>>&);

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
  // uint32_t
  int numSubsteps;
  std::string preconditioner;
  std::string solverType;

  // the HO Solver
  DummySolver::Solver solver;

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

  solver.set_re_d_(100.0);
  solver.set_elemorder_d_(elementOrder);
  solver.set_timeorder_d_(timeOrder);
  solver.set_numsteps_d_(numSubsteps);

  for (auto &coll : _euler) {
    // transform to current position
    coll.move(0.0, 0.0);

    // call the external solver with the current geometry
    // this will calculate the Jacobian and other cell-specific properties
    //(void) solver.init_d_(coll.get_pos(), coll.get_idx(),
    //                        coll.get_wall_idx(), coll.get_open_idx());

    // and ask it for the open BC solution nodes, and the full internal solution nodes
    //coll.set_open_pts(solver.getopenpts_d_());
    //coll.set_soln_pts(solver.getsolnpts_d_());
  }

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
                         std::vector<HOVolumes<S>>&           _euler) {

  if (not active) return;

  if (not initialized) init(_euler);

  std::cout << "Inside Hybrid::step at t=" << _time << " and dt=" << _dt << std::endl;

  // part A - prepare BCs for Euler solver

  // transform all elements to their appropriate locations (do we need to do this?)
  //for (auto &coll : _vort) {
  //}

  // update the BEM solution
  solve_bem<S,A,I>(_time, _fs, _vort, _bdry, _bem);
  // get vels on each euler region
  for (auto &coll : _euler) {
    // transform to current position
    coll.move(_time, _dt);

    // isolate open/outer boundaries
    //Points<A> euler_bdry = std::visit([=](auto& elem) { elem.get_bc_nodes(_time); }, coll);
    Points<S> euler_bdry = coll.get_bc_nodes(_time);

    // find vels there
    //Convection<S,A,I>::find_vels(_fs, _vort, _bdry, euler_bdry);

    // convert to transferable packet
    //std::array<Vector<S>,Dimensions> openvels = euler_bdry.get_vel();

    // transfer BC packet to solver
    //(void) solver.setopenvels_d_(coll.get_n(), openvels[0], openvels[1]);
  }

  // part B - call Euler solver

  // call solver
  (void) solver.solveto_d_(_time);

  // pull results from external solver
  //std::vector<A> allvorts;
  //(void) solver.getallvorts_d_(np, allvorts);

  // convert vorticity results to drawable and writable Volume elements

  // part C - update particle strengths accordionly

  // "update" strengths on vortex particles within the boundary

  // here's one way to do it:
  // identify all free vortex particles inside of euler regions and remove them
  //   (add up how much circulation we remove) - or not?

  // re-run BEM and compute vorticity on all HO volume nodes - bem already run
  for (auto &coll : _euler) {
    Points<S> euler_vol = coll.get_vol_nodes(_time);
    //Convection<S,A,I>::find_vels(_fs, _vort, _bdry, euler_vol, velandvort);
  }

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
  ImGui::SliderInt("Element Order", &elementOrder, 1, 5);

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
