/*
 * Simulation.h - a class to control a 2D vortex particle sim
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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
#include "Body.h"
#include "Collection.h"
#include "BEM.h"
#include "Convection.h"
#include "Diffusion.h"
#if defined(HOFORTRAN) || defined(HOCXX)
#include "HOVolumes.h"
#include "Hybrid.h"
#endif
#include "ElementPacket.h"
#include "StatusFile.h"

#ifdef USE_GL
#include "RenderParams.h"
#endif

#include <string>
#include <vector>
#include <future>

#ifdef USE_VC
#define STORE float
#define ACCUM float
#else
#define STORE float
#define ACCUM double
#endif


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
  float get_re() const;
  float get_dt() const;
  float get_ips() const;
  float get_hnu() const;
  float get_vdelta() const;
  float get_time() const;
  float get_end_time() const;
  bool using_end_time() const;
  size_t get_nstep() const;
  size_t get_max_steps() const;
  bool using_max_steps() const;
  float get_output_dt() const;
  std::string get_description();
  bool autostart();
  bool quitonstop();

  // setters
  void set_description(const std::string);
  void set_end_time(const double);
  void unset_end_time();
  void set_max_steps(const size_t);
  void unset_max_steps();
  void set_output_dt(const double);
  void set_auto_start(const bool);
  void set_quit_on_stop(const bool);

  // access status file
  void set_status_file_name(const std::string);
  std::string get_status_file_name();

  // get runtime status
  size_t get_npanels();
  size_t get_nparts();
  size_t get_nfldpts();

  // inviscid case needs this
  void set_re_for_ips(float);

  // receive and add a set of elements
  void add_elements(const ElementPacket<float>, const elem_t, const move_t, std::shared_ptr<Body>);
  void file_elements(std::vector<Collection>&, const ElementPacket<float>, const elem_t, const move_t, std::shared_ptr<Body>);
  void add_hybrid(const std::vector<ElementPacket<float>>, std::shared_ptr<Body>);

  // access body list
  void add_body(std::shared_ptr<Body>);
  std::shared_ptr<Body> get_last_body();
  std::shared_ptr<Body> get_pointer_to_body(const std::string);
  std::vector<std::shared_ptr<Body>>::iterator bodies_begin() { return bodies.begin(); }
  std::vector<std::shared_ptr<Body>>::iterator bodies_end() { return bodies.end(); }

  // get/set vrm/amr triggers
  void set_amr(const bool);
  void set_diffuse(const bool);
  const bool get_amr() const { return diff.get_amr(); };
  const bool get_diffuse() const { return diff.get_diffuse(); };

  // act on stuff
  void reset();
  void clear_bodies();
  void async_first_step();
  void first_step();
  void async_step();
  void step();
  void dump_stats_to_status();
  std::array<float,Dimensions> calculate_simple_forces();
  bool is_initialized();
  void set_initialized();
  std::string check_initialization();
  std::string check_simulation();
  bool do_any_bodies_move();
  bool any_nonzero_bcs();
  void conserve_iolet_volume();
  bool test_for_new_results();
  std::vector<std::string> write_vtk(const int _index = -1,
                                     const bool _do_bdry = true,
                                     const bool _do_flow = true,
                                     const bool _do_measure = true);
  bool test_vs_stop();
  bool test_vs_stop_async();

  // read to and write from a json object
  void from_json(const nlohmann::json);
  nlohmann::json to_json() const;
  void flow_from_json(const nlohmann::json);
  nlohmann::json flow_to_json() const;
  void runtime_from_json(const nlohmann::json);
  nlohmann::json runtime_to_json() const;

#ifdef USE_GL
  // graphics pass-through calls
  void updateGL();
  void drawGL(std::vector<float>&, RenderParams&);
#endif

#ifdef USE_IMGUI
  void draw_advanced();
#endif

private:
  // primary simulation params
  float re;
  float dt;
  float fs[Dimensions];

  // List of independent rigid bodies (and one for ground)
  std::vector< std::shared_ptr<Body> > bodies;

  // Object to contain all Lagrangian elements
  std::vector<Collection> vort;		// active elements

  // Object to contain all Reactive elements
  //   inside is the vector of bodies and inlets and liftinglines/kuttapoints
  //   and the Surfaces list of all unknowns discretized representations
  std::vector<Collection> bdry;		// reactive-active elements like BEM surfaces

  // Object with all of the non-reactive, non-active (inert) points
  std::vector<Collection> fldpt;	// tracers and field points

#if defined(HOFORTRAN) || defined(HOCXX)
  // Object with all of the near-body, Eulerian solution regions
  std::vector<HOVolumes<STORE>> euler;	// hybrid/euler regions
#endif

  // The need to solve for the unknown strengths of reactive elements inside both the
  //   diffusion and convection steps necessitates a BEM object here
  BEM<STORE,Int> bem;

  // Diffusion will resolve exchange of strength among particles and between panels and particles
  // Note that NNLS needs doubles for its compute type or else it will fail
  //   but velocity evaluations using Vc need S=A=float
  Diffusion<STORE,ACCUM,Int> diff;

  // Convection class takes bodies, panels, and vector of Particles objects
  //   and performs 1st, 2nd, etc. RK forward integration
  //   inside here are non-drawing Particles objects used as temporaries in the multi-step methods
  //   also copies of the panels, which will be recreated for each step, and solutions to unknowns
  // Note that with Vc, the storage and accumulator classes have to be the same
  Convection<STORE,ACCUM,Int> conv;

#if defined(HOFORTRAN) || defined(HOCXX)
  // Hybrid class controls data exchange to and from an external Euler solver
  //   to be used for near-body regions, and routines to affect existing vortex
  //   particles from grid-based vorticity
  Hybrid<STORE,ACCUM,Int> hybr;
#endif

  // status file
  StatusFile sf;

  // state
  std::string description;
  double time;
  double output_dt;
  double end_time;
  bool use_end_time;
  float overlap_ratio;
  float core_size_ratio;
  size_t nstep;
  bool use_max_steps;
  size_t max_steps;
  bool auto_start;
  bool quit_on_stop;
  bool sim_is_initialized;
  bool step_has_started;
  bool step_is_finished;
  std::future<void> stepfuture;  // this future needs to be listed after the big four: diff, conv, ...
};

