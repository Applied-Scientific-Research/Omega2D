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
#include "dummysolver.h"

#include <iostream>
#include <vector>
//#include <numeric>	// for transform_reduce (C++17)

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
  void first_step( const double,
             const std::array<double,Dimensions>&,
             std::vector<Collection>&,
             std::vector<Collection>&,
             BEM<S,I>&,
             Convection<S,A,I>&,
             std::vector<HOVolumes<S>>&);
  void step( const double,
             const double,
             const std::array<double,Dimensions>&,
             std::vector<Collection>&,
             std::vector<Collection>&,
             BEM<S,I>&,
             Convection<S,A,I>&,
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
  int timeOrder;
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
  solver.set_elemorder_d_((uint8_t)elementOrder);
  solver.set_timeorder_d_((uint8_t)timeOrder);
  solver.set_numsteps_d_((uint32_t)numSubsteps);

  for (auto &coll : _euler) {
    // transform to current position
    coll.move(0.0, 0.0);

    // call the external solver with the current geometry
    // this will calculate the Jacobian and other cell-specific properties
    (void) solver.init_d_(coll.get_node_pos(), coll.get_elem_idx(),
                          coll.get_wall_idx(), coll.get_open_idx());

    // and ask it for the open BC solution nodes, and the full internal solution nodes
    coll.set_open_pts(solver.getopenpts_d_());
    coll.set_soln_pts(solver.getsolnpts_d_());
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
    (void) solver.setopenvels_d_(packedvels);
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

  // get vels and vorts on each euler region - and force it
  _conv.find_vels(_fs, _vort, _bdry, euler_vols, velandvort, true);

  for (auto &coll : euler_vols) {
    // convert to transferable packet
    Vector<S> volvort = std::visit([=](auto& elem) { return elem.get_vort(); }, coll);

    // convert to a std::vector<double>
    std::vector<double> vorts(volvort.begin(), volvort.end());

    // transfer BC packet to solver
    (void) solver.setsolnvort_d_(vorts);
  }
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
                         Convection<S,A,I>&                   _conv,
                         std::vector<HOVolumes<S>>&           _euler) {

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
    (void) solver.setopenvels_d_(packedvels);
  }


  //
  // part B - call Euler solver
  //

  // call solver - solves all Euler volumes at once?
  (void) solver.solveto_d_(_time);

  // pull results from external solver (assume just one for now)
  //for (auto &coll : _euler) {
  //  std::vector<double> allvorts = solver.getallvorts_d_();
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

  for (auto &coll : _euler) {

    Points<S>& solnpts = coll.get_vol_nodes(_time);
    const size_t thisn = solnpts.get_n();

    // pull results from external solver (assume just one for now)
    std::vector<double> eulvort = solver.getallvorts_d_();
    assert(eulvort.size() == thisn && "ERROR (Hybrid::step) vorticity from solver is not the right size");

    // find the Lagrangian-computed vorticity on all solution nodes (make a vector of one collection)
    std::vector<Collection> euler_vols;
    euler_vols.emplace_back(solnpts);	// this is a COPY
    _conv.find_vels(_fs, _vort, _bdry, euler_vols, velandvort, true);
    //Vector<S>& lagvort = solnpts.get_vort();	// can't do this, that object doesn't have the results!
    Points<float>& solvedpts = std::get<Points<S>>(euler_vols[0]);
    Vector<S>& lagvort = solvedpts.get_vort();
    assert(lagvort.size() == thisn && "ERROR (Hybrid::step) vorticity from particle sim is not the right size");


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
    // HACK - when solution points are no longer 1:1 with volume elements, this will need work
    const Vector<S>& area = coll.get_maskarea();
    assert(area.size() == thisn && "ERROR (Hybrid::step) volume area vector is not the right size");
    std::transform(circ.begin(), circ.end(),
                   area.begin(),
                   circ.begin(),
                   std::multiplies<S>());

    for (size_t i=0; i<std::min(thisn,(size_t)60); ++i) {
      std::cout << "    " << i << "  " << area[i] << " * ( " << eulvort[i] << " - " << lagvort[i] << " ) = " << circ[i] << std::endl;
    }

    // measure this error
    double maxerror = 0.1;
    double thiserror = 0.0;
    // trying to get fancy - failed!
    //double thiserror = std::accumulate(circ.begin(), circ.end(), 0.0);
    //double thiserror = std::transform_reduce(circ.begin(), circ.end(), 0.0, std::plus<>{}, //std::fabs);
    //                                         static_cast<double (*)(double)>(std::fabs));
    for (size_t i=0; i<thisn; ++i) { thiserror += std::fabs(circ[i]); }
    std::cout << "  initial error " << thiserror << std::endl;

    // iterate toward a solution
    int iter = 0;
    int maxiter = 10;
    while (thiserror>maxerror and iter<maxiter) {

      // resolve it via VRM
      //   for each sub-node of each HO quad, run a VRM onto the existing set of Lagrangian
      //   particles adding where necessary to satisfy at least 0th and 1st moments, if not 2nd
      // resolve it via naive insertion
      //   just create a new particle at the centroid of each element, let merge deal

      // generate particles to fill the deficit
      //   note that cells may be a lot larger than particles!

      // merge here? seems smart

      // find the new Lagrangian vorticity on the solution nodes
      // and the new vorticity error/deficit

      // measure this error
      thiserror = 0.0;
      std::cout << "  iter " << (iter+1) << " has error " << thiserror << std::endl;

      iter++;
    }
  }

  // then merge all particles that we created

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
