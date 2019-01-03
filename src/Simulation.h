/*
 * Simulation.h - a class to control a 2d vortex particle sim
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Body.h"
#include "Collection.h"
#include "Boundaries.h"
#include "BoundaryFeature.h"
#include "Convection.h"
#include "Diffusion.h"
#include "Vorticity.h"

#include <string>
#include <vector>
#include <future>
#include <chrono>

template <class T>
bool is_future_ready(std::future<T> const& f) {
    if (!f.valid()) return false;
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

//
// A set of particles, can be sources or targets
//
class Simulation {
public:
  Simulation();

  // imgui needs access to memory locations
  float* addr_re();
  float* addr_dt();
  float* addr_fs();

  // get the derived parameters
  float get_ips();
  float get_hnu();
  float get_vdelta();
  float get_time();

  // get runtime status
  size_t get_npanels();
  size_t get_nparts();

  // inviscid case needs this
  void set_re_for_ips(float);

  // receive and add a set of particles
  void add_particles(std::vector<float>);
  void add_tracers(std::vector<float>);
  void add_boundary(bdryType, std::vector<float>);

  // act on stuff
  //void set_amr(const bool);
  void set_diffuse(const bool);
  void reset();
  void async_step();
  void step();
  void init_bcs();
  bool is_initialized();
  void set_initialized();
  std::string check_simulation(const size_t, const size_t);
  bool test_for_new_results();

  // graphics pass-through calls
  void initGL(std::vector<float>&, float*, float*);
  void updateGL();
  void drawGL(std::vector<float>&, float*, float*);

private:
  // primary simulation params
  float re;
  float dt;
  float fs[2];

  // Object to contain all Lagrangian elements
  Vorticity<float,uint16_t> vort;
  std::vector<Collection> vort2;	// active elements

  // Object to contain all Reactive elements
  //   inside is the vector of bodies and inlets and liftinglines/kuttapoints
  //   and the Panels list of all unknowns discretized representations
  Boundaries<float,uint16_t> bdry;
  std::vector<Collection> bdry2;	// reactive-active elements like BEM surfaces

  // Object with all of the non-reactive, non-active (inert) points
  std::vector<Collection> fldpt;	// tracers and field points

  // Diffusion will resolve exchange of strength among particles and between panels and particles
  // Note that NNLS needs doubles for its compute type or else it will fail
  //   but velocity evaluations using Vc need S=A=float
  Diffusion<float,double,uint16_t> diff;

  // Convection class takes bodies, panels, and vector of Particles objects
  //   and performs 1st, 2nd, etc. RK forward integration
  //   inside here are non-drawing Particles objects used as temporaries in the multi-step methods
  //   also copies of the panels, which will be recreated for each step, and solutions to unknowns
  // Note that with Vc, the storage and accumulator classes have to be the same
#ifdef USE_VC
  Convection<float,float,uint16_t> conv;
#else
  Convection<float,double,uint16_t> conv;
#endif

  // state
  double time;
  bool sim_is_initialized;
  bool step_has_started;
  bool step_is_finished;
  std::future<void> stepfuture;  // this future needs to be listed after the big four: diff, conv, ...
};

