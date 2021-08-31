/*
 * Simulation.cpp - a class to control a 2D vortex particle sim
 *
 * (c)2017-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#include "Simulation.h"
#include "Reflect.h"
#include "BEMHelper.h"
#include "GuiHelper.h"
#include "FutureHelper.h"

#include <cassert>
#include <cmath>
#include <cfenv> // Catch fp exceptions
#include <limits>
#include <variant>

#ifdef _WIN32
#pragma STDC FENV_ACCESS ON // For fp exceptions
#endif


// constructor
Simulation::Simulation()
  : re(100.0),
    dt(0.01),
    fs{0.0,0.0},
    bodies(),
    vort(),
    bdry(),
    fldpt(),
#if defined(HOFORTRAN) || defined(HOCXX)
    euler(),
#endif
    bem(),
    diff(),
    conv(),
#if defined(HOFORTRAN) || defined(HOCXX)
    hybr(),
#endif
    sf(),
    description(),
    time(0.0),
    output_dt(0.0),
    end_time(100.0),
    use_end_time(false),
    overlap_ratio(2.0),
    core_size_ratio(std::sqrt(6)),
    nstep(0),
    use_max_steps(false),
    max_steps(100),
    auto_start(false),
    quit_on_stop(false),
    sim_is_initialized(false),
    step_has_started(false),
    step_is_finished(false),
    stepfuture()
  {}

// addresses for use in imgui
float* Simulation::addr_re() { return &re; }
float* Simulation::addr_dt() { return &dt; }
float* Simulation::addr_fs() { return fs; }

// getters
float Simulation::get_re() const { return re; }
float Simulation::get_dt() const { return dt; }
float Simulation::get_hnu() const { return std::sqrt(dt/re); }
float Simulation::get_ips() const { return core_size_ratio * get_hnu(); }
float Simulation::get_vdelta() const { return overlap_ratio * get_ips(); }
float Simulation::get_time() const { return (float)time; }
float Simulation::get_end_time() const { return (float)end_time; }
bool Simulation::using_end_time() const { return use_end_time; }
size_t Simulation::get_nstep() const { return nstep; }
size_t Simulation::get_max_steps() const { return max_steps; }
bool Simulation::using_max_steps() const { return use_max_steps; }
float Simulation::get_output_dt() const { return (float)output_dt; }
std::string Simulation::get_description() { return description; }
bool Simulation::autostart() { return auto_start; }
bool Simulation::quitonstop() { return quit_on_stop; }

// setters
void Simulation::set_description(const std::string desc) { description = desc; }
void Simulation::set_end_time(const double net) { end_time = net; use_end_time = true; }
void Simulation::unset_end_time() { use_end_time = false; }
void Simulation::set_max_steps(const size_t nms) { max_steps = nms; use_max_steps = true; }
void Simulation::unset_max_steps() { use_max_steps = false; }
void Simulation::set_output_dt(const double nodt) { output_dt = nodt; }
void Simulation::set_auto_start(const bool autos) { auto_start = autos; }
void Simulation::set_quit_on_stop(const bool qos) { quit_on_stop = qos; }

// access status file
void Simulation::set_status_file_name(const std::string _fn) { sf.set_filename(_fn); }
std::string Simulation::get_status_file_name() { return sf.get_filename(); }

// status
size_t Simulation::get_npanels() {
  size_t n = 0;
  for (auto &coll: bdry) {
    //std::visit([&n](auto& elem) { n += elem.get_npanels(); }, coll);
    // only proceed if the last collection is Surfaces
    if (std::holds_alternative<Surfaces<STORE>>(coll)) {
      Surfaces<STORE>& surf = std::get<Surfaces<STORE>>(coll);
      n += surf.get_npanels();
    }
  }
  return n;
}

size_t Simulation::get_nparts() {
  size_t n = 0;
  for (auto &coll: vort) {
    std::visit([&n](auto& elem) { n += elem.get_n(); }, coll);
  }
  return n;
}

size_t Simulation::get_nfldpts() {
  size_t n = 0;
  for (auto &coll : fldpt) {
    std::visit([&n](auto& elem) { n += elem.get_n(); }, coll);
  }
  return n;
}

// like a setter
void Simulation::set_re_for_ips(const float _ips) {
  re = std::pow(core_size_ratio, 2) * dt / pow(_ips, 2);
  diff.set_diffuse(false);
}

void Simulation::set_diffuse(const bool _do_diffuse) {
  diff.set_diffuse(_do_diffuse);
}

void Simulation::set_amr(const bool _do_amr) {
  diff.set_amr(_do_amr);
}


//
// json read/write
//

// read "simparams" json object
void
Simulation::from_json(const nlohmann::json j) {

  if (j.find("nominalDt") != j.end()) {
    dt = j["nominalDt"];
    std::cout << "  setting dt= " << dt << std::endl;
  }

  if (j.find("outputDt") != j.end()) {
    output_dt = j["outputDt"];
    std::cout << "  setting output dt= " << output_dt << std::endl;
  }

  //if (j.find("nominalDx") != j.end()) {
  //  dx = j["nominalDx"];
  //  std::cout << "  setting dx= " << dx << std::endl;
  //}

  if (j.find("maxSteps") != j.end()) {
    use_max_steps = true;
    max_steps = j["maxSteps"];
    std::cout << "  setting max_steps= " << max_steps << std::endl;
  } else {
    use_max_steps = false;
  }

  if (j.find("endTime") != j.end()) {
    use_end_time = true;
    end_time = j["endTime"];
    std::cout << "  setting end_time= " << end_time << std::endl;
  } else {
    use_end_time = false;
  }

  if (j.find("overlapRatio") != j.end()) {
    overlap_ratio = j["overlapRatio"];
    std::cout << "  setting overlap ratio= " << overlap_ratio << std::endl;
  }

  if (j.find("coreSizeRatioSqrd") != j.end()) {
    core_size_ratio = std::sqrt((float)j["coreSizeRatioSqrd"]);
    std::cout << "  setting core size ratio (nominal separation over h_nu) = " << core_size_ratio << std::endl;
  }

  // Convection will find and set "timeOrder"
  conv.from_json(j);

  // Diffusion will find and set "viscous", "VRM" and "AMR" parameters
  diff.from_json(j);

#if defined(HOFORTRAN) || defined(HOCXX)
  // set hybrid Eulerian-Lagrangian solution parameters
  hybr.from_json(j);
#endif
}

// create and write a json object for "simparams"
nlohmann::json
Simulation::to_json() const {
  nlohmann::json j;

  j["nominalDt"] = dt;
  j["outputDt"] = output_dt;
  if (using_max_steps()) j["maxSteps"] = get_max_steps();
  if (using_end_time()) j["endTime"] = get_end_time();
  j["overlapRatio"] = overlap_ratio;
  j["coreSizeRatioSqrd"] = std::pow(core_size_ratio,2);

  // Convection will write "timeOrder"
  conv.add_to_json(j);

  // Diffusion will write "viscous", "VRM" and "AMR" parameters
  diff.add_to_json(j);

#if defined(HOFORTRAN) || defined(HOCXX)
  // Hybrid will create a "hybrid" section
  hybr.add_to_json(j);
#endif

  return j;
}

// set "flowparams" json object
void
Simulation::flow_from_json(const nlohmann::json j) {

  if (j.find("Re") != j.end()) {
    re = j["Re"];
    std::cout << "  setting re= " << re << std::endl;
  }
  if (j.find("Uinf") != j.end()) {
    // eventually support an expression for Uinf instead of just a single float
    std::vector<float> new_fs = {0.0, 0.0, 0.0};
    new_fs.resize(Dimensions);
    if (j["Uinf"].is_array()) {
      new_fs = j["Uinf"].get<std::vector<float>>();
    } else if (j["Uinf"].is_number()) {
      new_fs[0] = j["Uinf"].get<float>();
    }
    for (size_t i=0; i<Dimensions; ++i) fs[i] = new_fs[i];
    std::cout << "  setting freestream to " << fs[0] << " " << fs[1] << std::endl;
  }
}

// create and write a json object for "flowparams"
nlohmann::json
Simulation::flow_to_json() const {
  nlohmann::json j;

  j["Re"] = re;
  j["Uinf"] = {fs[0], fs[1]};

  return j;
}

// set "runtime" json object
void
Simulation::runtime_from_json(const nlohmann::json j) {

  sf.from_json(j);

  if (j.find("autoStart") != j.end()) {
    bool autostart = j["autoStart"];
    set_auto_start(autostart);
    std::cout << "  autostart? " << autostart << std::endl;
  }
  if (j.find("quitOnStop") != j.end()) {
    bool qos = j["quitOnStop"];
    set_quit_on_stop(qos);
    std::cout << "  quit on stop? " << qos << std::endl;
  }
}

// create and write a json object for "runtime"
nlohmann::json
Simulation::runtime_to_json() const {
  nlohmann::json j;

  sf.add_to_json(j);

  j["autoStart"] = auto_start;
  j["quitOnStop"] = quit_on_stop;

  return j;
}


#ifdef USE_IMGUI
//
// ImGui code
//
void Simulation::draw_advanced() {

  // set the execution environment in Convection.h
  conv.draw_advanced();

  // set the diffusion parameters in Diffusion.h
  diff.draw_advanced();
  
#if defined(HOFORTRAN) || defined(HOCXX)
  // set the hybrid parameters in Hybrid.h
  hybr.draw_advanced();
#endif
}
#endif


#ifdef USE_GL
//
// OpenGL-specific code
//

void Simulation::updateGL() {
  for (auto &coll : vort) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
  for (auto &coll : bdry) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
  for (auto &coll : fldpt) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
}

void Simulation::drawGL(std::vector<float>& _projmat,
                        RenderParams&       _rparams) {

  if (step_is_finished) {
    _rparams.tracer_size = get_ips() * _rparams.tracer_scale;
    for (auto &coll : vort) {
      std::visit([&](auto& elem) { elem.drawGL(_projmat, _rparams, get_vdelta()); }, coll);
    }
    for (auto &coll : bdry) {
      std::visit([&](auto& elem) { elem.drawGL(_projmat, _rparams, get_vdelta()); }, coll);
    }
    for (auto &coll : fldpt) {
      std::visit([&](auto& elem) { elem.drawGL(_projmat, _rparams, get_vdelta()); }, coll);
    }
  }
}
#endif

bool Simulation::is_initialized() { return sim_is_initialized; }

void Simulation::set_initialized() { sim_is_initialized = true; }

void Simulation::reset() {

  // must wait for step() to complete, if it's still working
  if (stepfuture.valid()) {
    stepfuture.wait();
    stepfuture.get();
  }

  // now reset everything else
  time = 0.0;
  nstep = 0;
  vort.clear();
  bdry.clear();
  fldpt.clear();
  bem.reset();
#if defined(HOFORTRAN) || defined(HOCXX)
  hybr.reset(euler);
  euler.clear();
#endif
  sf.reset_sim();
  sim_is_initialized = false;
  step_has_started = false;
  step_is_finished = false;
}

void Simulation::clear_bodies() {
  bodies.clear();
}

// Write a set of vtu files for the particles and panels
std::vector<std::string> Simulation::write_vtk(const int _index,
                                               const bool _do_bdry,
                                               const bool _do_flow,
                                               const bool _do_measure) {

  std::cout << "Inside Simulation::write_vtk at t=" << time << std::endl;

  // solve the BEM (before any VTK or status file output)
  //std::cout << "Updating element vels" << std::endl;
  std::array<double,2> thisfs = {fs[0], fs[1]};
  //clear_inner_layer<STORE>(1, bdry, vort, 1.0/std::sqrt(2.0*M_PI), get_ips());
  solve_bem<STORE,ACCUM,Int>(time, thisfs, vort, bdry, bem);

  // special - only here do we cacluate the vorticity as well as velocity, true means force
  //if (_do_flow)    conv.find_vels(thisfs, vort, bdry, vort, velandshear, true);
  if (_do_flow)    conv.find_vels(thisfs, vort, bdry, vort, velandvort, true);
  if (_do_measure) conv.find_vels(thisfs, vort, bdry, fldpt, velandvort, true);
  if (_do_bdry)    conv.find_vels(thisfs, vort, bdry, bdry, velandvort, true);

  // may eventually want to avoid clobbering by maintaining an internal count of the
  //   number of simulations run from this execution of the GUI
  std::vector<std::string> files;

  // if a positive number was passed in, use that instead of the current step
  size_t stepnum = 0;
  if (_index < 0) {
    stepnum = nstep;
  } else {
    stepnum = (size_t)_index;
  }

  // ask Vtk to write files for each collection
  if (_do_flow) { 
    size_t idx = 0;
    for (auto &coll : vort) {
      std::visit([&](auto &&elem) { files.emplace_back(elem.write_vtk(idx++, stepnum, time)); }, coll);
    }
  }
  if (_do_measure) {
    size_t idx = 0;
    for (auto &coll : fldpt) {
      std::visit([&](auto &&elem) { files.emplace_back(elem.write_vtk(idx++, stepnum, time)); }, coll);
    }
  }
  if (_do_bdry) {
    size_t idx = 0;
    for (auto &coll : bdry) {
      std::visit([&](auto &&elem) { files.emplace_back(elem.write_vtk(idx++, stepnum, time)); }, coll);
    }
  }

  if (false) {
    // use the hack-y way to find vorticity and write it
    conv.find_vort(vort, bdry, fldpt);
    size_t idx = 0;
    for (auto &coll : fldpt) {
      std::visit([&](auto &&elem) { files.emplace_back(elem.write_vtk(idx++, stepnum+10, time)); }, coll);
    }
  }

#if defined(HOFORTRAN) || defined(HOCXX)
  if (hybr.is_active()) {
    // there is only one hybrid volume allowed now, so no counting needed
    hybr.trigger_write(stepnum, euler);
    files.emplace_back("and a HO grid file");
  }
#endif

  return files;
}

//
// Check all aspects of the initialization for conditions that prevent a run from starting
//
std::string Simulation::check_initialization() {
  std::string retstr;

  // are any flow features particle generators? - HOW DO WE DO THIS?
  const bool are_generators = false;

  // Check for no bodies and no particles
  if (get_npanels() == 0 and get_nparts() == 0 and not are_generators) {
    retstr.append("No flow features and no bodies. Add one or both, reset, and run.\n");
  }

  // Check for a body and no particles
  if (get_npanels() > 0 and get_nparts() == 0) {

    const bool zero_freestream = (fs[0]*fs[0]+fs[1]*fs[1] < std::numeric_limits<float>::epsilon());
    const bool no_body_movement = not do_any_bodies_move();
    const bool all_zero_bcs = not any_nonzero_bcs();

    // AND no freestream
    if (zero_freestream and no_body_movement and all_zero_bcs) {
      retstr.append("No flow features, zero freestream speed, no movement, and no driven boundaries - try adding one of these.\n");
      return retstr;
    }

    // AND no viscosity
    if (not diff.get_diffuse()) {
      retstr.append("You have a solid body, but no diffusion. It will not shed vorticity. Turn on viscosity or add a flow feature, reset, and run.\n");
    }
  }

  // Check for very large BEM problem
  if (get_npanels() > 21000) {
    retstr.append("Boundary features have too many panels, program will run out of memory. Reduce Reynolds number or increase time step or both.\n");
  }

  // adjust outlet speeds to conserve volume
  (void)conserve_iolet_volume();

  return retstr;
}

//
// Check dynamic aspects of the simulation for conditions that should stop the run
//
std::string Simulation::check_simulation() {
  std::string retstr;

  // Are there any dynamic problems in 2D that could blow a run?

  return retstr;
}

//
// check all bodies for movement
//
bool Simulation::do_any_bodies_move() {
  bool some_move = false;
  for (size_t i=0; i<bodies.size(); ++i) {
    auto thisvel = bodies[i]->get_vel(time);
    auto nextvel = bodies[i]->get_vel(time+dt);
    auto thisrot = bodies[i]->get_rotvel(time);
    auto nextrot = bodies[i]->get_rotvel(time+dt);
    if (std::abs(thisvel[0]) + std::abs(thisvel[1]) + std::abs(thisrot) +
        std::abs(nextvel[0]) + std::abs(nextvel[1]) + std::abs(nextrot) >
        std::numeric_limits<float>::epsilon()) {
      some_move = true;
    }
  }
  return some_move;
}

//
// check all bodies for non-zero BCs
//
bool Simulation::any_nonzero_bcs() {
  bool all_are_zero = true;
  // loop over bdry Collections
  for (auto &src : bdry) {
    float max_bc = std::visit([=](auto& elem) { return elem.get_max_bc_value(); }, src);
    if (std::abs(max_bc) > std::numeric_limits<float>::epsilon()) all_are_zero = false;
  }
  return not all_are_zero;
}

//
// check all inlets and outlets to ensure volume conservation
//
void Simulation::conserve_iolet_volume() {

  // current inlet and outlet volume flow rates (global summation, includes the Surfaces
  //   from HOVolumes loaded in via meshes but used on the Lagrangian side)
  float inrate = 0.0;
  float outrate = 0.0;
  int num_tested = 0;

  // loop over bdry Collections and find total inlet and outlet rates
  for (auto &src : bdry) {
    if (std::holds_alternative<Surfaces<STORE>>(src)) {
      Surfaces<STORE>& surf = std::get<Surfaces<STORE>>(src);
      if (surf.get_elemt() == reactive) {
        inrate += surf.get_total_inflow();
        outrate += surf.get_total_outflow();
        num_tested++;
      }
    }
  }
  if (num_tested == 0) return;

  std::cout << "  total lagrangian sim inflow, outflow " << inrate << " " << outrate << std::endl;

  if (inrate > std::numeric_limits<float>::epsilon()) {
    if (outrate > std::numeric_limits<float>::epsilon()) {
      // there is finite inflow and outflow, ensure their magnitudes are matched

      std::cout << "    scaling outflows by " << inrate/outrate << std::endl;

      // correct the outflow to match
      for (auto &src : bdry) {
        if (std::holds_alternative<Surfaces<STORE>>(src)) {
          Surfaces<STORE>& surf = std::get<Surfaces<STORE>>(src);
          if (surf.get_elemt() == reactive) {
            surf.scale_outflow(inrate/outrate);
          }
        }
      }

    } else {
      // there is inflow, but zero outflow - just zero the inflow and continue

      std::cout << "    no outflows defined - zeroing all inflows" << std::endl;

      for (auto &src : bdry) {
        if (std::holds_alternative<Surfaces<STORE>>(src)) {
          Surfaces<STORE>& surf = std::get<Surfaces<STORE>>(src);
          if (surf.get_elemt() == reactive) {
            surf.zero_inflow();
          }
        }
      }
    }
  }

  return;
}

//
// query and get() the future if possible
//
bool Simulation::test_for_new_results() {

  if (not step_has_started) {
    // if we haven't made an async call yet
    return true;

  } else if (is_future_ready(stepfuture)) {
    // if we did, and it's ready to give us the results
    stepfuture.get();

#ifdef USE_GL
    // tell flow objects to update their values to the GPU
    updateGL();
#endif

    // set flag indicating that at least one step has been solved
    step_is_finished = true;
    step_has_started = false;

    return true;
  }

  // async call is not finished, do not try calling it again
  return false;
}

//
// call this from a real-time GUI - will only start the first step if no other steps are being worked on
//
void Simulation::async_first_step() {
  step_has_started = true;
  stepfuture = std::async(std::launch::async, [this](){first_step();});
}

//
// initialize the system so we can start drawing things
//
void Simulation::first_step() {
  std::cout << std::endl << "Taking step " << nstep << " at t=" << time << std::endl;

  // we wind up using this a lot
  std::array<double,2> thisfs = {fs[0], fs[1]};

  // this is the first step, just solve BEM and return - it's time=0

  // update BEM and find vels on any particles but DO NOT ADVECT
  conv.advect_1st(time, 0.0, thisfs, get_ips(), vort, bdry, fldpt, bem);

#if defined(HOFORTRAN) || defined(HOCXX)
  // call HO grid solver, but only to send first velocity results and initialize vorticity
  hybr.first_step(time, thisfs, vort, bdry, bem, conv, euler);
#endif

  // and write status file
  dump_stats_to_status();
}

//
// call this from a real-time GUI - will only run a new step if the last one is done
//
void Simulation::async_step() {
  step_has_started = true;
  stepfuture = std::async(std::launch::async, [this](){step();});
}

//
// here's the vortex method: convection and diffusion with operator splitting
//
void Simulation::step() {
  // Catch Division by 0
  // unsigned int current_word = 0;
  // _controlfp_s(&current_word, _EM_UNDERFLOW | _EM_OVERFLOW | _EM_INEXACT, _MCW_EM);

  std::cout << std::endl << "Taking step " << nstep << " at t=" << time << " with n=" << get_nparts() << std::endl;

  const bool use_2nd_order_operator_splitting = true;

  // we wind up using this a lot
  std::array<double,2> thisfs = {fs[0], fs[1]};

#ifdef PLUGIN_AVRM
  // calculate shear rates for adaptive vrm, and force it
  if (diff.get_amr()) conv.find_vels(thisfs, vort, bdry, vort, velandshear, true);
#endif

  if (use_2nd_order_operator_splitting) {
    // operator splitting requires one half-step diffuse (use coefficients from previous step, if available)
    diff.step(time, 0.5*dt, re, overlap_ratio, get_vdelta(), thisfs, vort, bdry, bem);
  } else {
    // for simplicity's sake, just run one full diffusion step here
    diff.step(time, dt, re, overlap_ratio, get_vdelta(), thisfs, vort, bdry, bem);
  }

  // advect with no diffusion (must update BEM strengths)
  conv.advect(time, dt, thisfs, get_ips(), vort, bdry, fldpt, bem);

  if (use_2nd_order_operator_splitting) {
#ifdef PLUGIN_AVRM
    // calculate shear rates for adaptive vrm, and force it
    if (diff.get_amr()) conv.find_vels(thisfs, vort, bdry, vort, velandshear, true);
#endif

    // operator splitting requires another half-step diffuse (must compute new coefficients)
    diff.step(time+dt, 0.5*dt, re, overlap_ratio, get_vdelta(), thisfs, vort, bdry, bem);
  }

#if defined(HOFORTRAN) || defined(HOCXX)
  // call HO grid solver to recalculate vorticity at the end of this time step
  hybr.step(time, dt, re, thisfs, vort, bdry, bem, conv, euler, overlap_ratio, get_vdelta());
#endif

  // update time
  time += (double)dt;

  // push field points out of objects every few steps
  if (nstep%5 == 0) clear_inner_layer<STORE>(1, bdry, fldpt, (STORE)0.0, (STORE)(0.5*get_ips()));

  // only increment step here!
  nstep++;

  // and write status file
  dump_stats_to_status();
}

//
// close out the step with some work and output to the status file
//
void Simulation::dump_stats_to_status() {
  if (sf.is_active()) {
    // the basics
    sf.append_value("time",(float)time);
    sf.append_value("Nv",(int)get_nparts());

    // more advanced info

    std::array<double,2> thisfs = {fs[0], fs[1]};

    // push away particles inside or too close to the body
    //clear_inner_layer<STORE>(1, bdry, vort, 1.0/std::sqrt(2.0*M_PI), get_ips());
    // solve the BEM (before any VTK or status file output)
    solve_bem<STORE,ACCUM,Int>(time, thisfs, vort, bdry, bem);

    // but do we really need to do these?
    //conv.find_vels(thisfs, vort, bdry, vort);
    //conv.find_vels(thisfs, vort, bdry, fldpt);
    //conv.find_vels(thisfs, vort, bdry, bdry);

    // add up the total circulation
    float tot_circ = 0.0;
    for (auto &src : vort) {
      tot_circ += std::visit([=](auto& elem) { return elem.get_total_circ(time); }, src);
    }
    // then add up the circulation in bodies
    for (auto &src : bdry) {
      tot_circ += std::visit([=](auto& elem) { return elem.get_total_circ(time); }, src);
      tot_circ += std::visit([=](auto& elem) { return elem.get_body_circ(time); }, src);
    }
    sf.append_value("circ",tot_circ);

    // now forces
    std::array<float,Dimensions> impulse = calculate_simple_forces();
    sf.append_value("fx",impulse[0]);
    sf.append_value("fy",impulse[1]);
    if (Dimensions > 2) sf.append_value("fz",impulse[3]);

    // write here
    sf.write_line();
  }
}

// Use impulse method to calculate total forces
std::array<float,Dimensions>
Simulation::calculate_simple_forces() {

  static double last_time = 0.0;
  static std::array<float,Dimensions> last_impulse = {0.0};
  std::array<float,Dimensions> this_impulse = {0.0};

  // calculate impulse from particles
  for (auto &src : vort) {
    std::array<STORE,Dimensions> this_imp = std::visit([=](auto& elem) { return elem.get_total_impulse(); }, src);
    for (size_t i=0; i<Dimensions; ++i) this_impulse[i] += this_imp[i];
  }
  // then add up the impulse from bodies - DO WE NEED TO RE-SOLVE BEM FIRST?
  for (auto &src : bdry) {
    std::array<STORE,Dimensions> this_imp = std::visit([=](auto& elem) { return elem.get_total_impulse(); }, src);
    for (size_t i=0; i<Dimensions; ++i) this_impulse[i] += this_imp[i];
  }

  // reset the "last" values if time is zero
  if (time < 0.1*dt) {
    last_time = -dt;
    last_impulse = this_impulse;
  }

  // find the time derivative of the impulses
  std::array<float,Dimensions> forces;
  for (size_t i=0; i<Dimensions; ++i) forces[i] = (this_impulse[i] - last_impulse[i]) / (time - last_time);

  // save the last condition
  last_time = time;
  for (size_t i=0; i<Dimensions; ++i) last_impulse[i] = this_impulse[i];

  return forces;
}

// Add elements - any kind
void Simulation::add_elements(const ElementPacket<float> _elems,
                              const elem_t _et, const move_t _mt,
                              std::shared_ptr<Body> _bptr) {

  // skip out early if nothing's here
  if (_elems.nelem == 0) return;

  // now split on which Collection will receive this
  if (_et == active) {
    // it's active vorticity, add to vort
    file_elements(vort, _elems, active, _mt, _bptr);
    // in that routine, we will look for a match for move type, body pointer, and points/surfs/vols
  } else if (_et == reactive) {
    file_elements(bdry, _elems, reactive, _mt, _bptr);
  } else {
    file_elements(fldpt, _elems, inert, _mt, _bptr);
  }
}

// File the new elements into the correct collection
void Simulation::file_elements(std::vector<Collection>& _collvec,
                               const ElementPacket<float> _elems,
                               const elem_t _et, const move_t _mt,
                               std::shared_ptr<Body> _bptr) {

  // search the collections list for a match (same movement type, Body, elem dims)
  size_t imatch = 0;
  bool no_match = true;
  for (size_t i=0; i<_collvec.size(); ++i) {
    // assume match
    bool this_match = true;

    // check movement type
    const move_t tmt = std::visit([=](auto& elem) { return elem.get_movet(); }, _collvec[i]);
    if (_mt != tmt) {
      this_match = false;

    } else if (tmt == bodybound) {
      // check body pointer
      std::shared_ptr<Body> tbp = std::visit([=](auto& elem) { return elem.get_body_ptr(); }, _collvec[i]);
      if (_bptr != tbp) this_match = false;
    }

    // check Collections element dimension
    auto& coll = _collvec[i];
    if (std::holds_alternative<Points<STORE>>(coll) and _elems.ndim != 0) {
      this_match = false;
    } else if (std::holds_alternative<Surfaces<STORE>>(coll) and _elems.ndim != 1) {
      this_match = false;
    } else if (std::holds_alternative<Volumes<STORE>>(coll) and _elems.ndim != 2) {
      this_match = false;
    }

    if (this_match) {
      imatch = i;
      no_match = false;
    }
  }

  // if no match, or no collections exist
  if (no_match) {
    // make a new collection according to element dimension
    if (_elems.ndim == 0) {
      _collvec.push_back(Points<STORE>(_elems, _et, _mt, _bptr, get_vdelta()));
    } else if (_elems.ndim == 1) {
      _collvec.push_back(Surfaces<STORE>(_elems, _et, _mt, _bptr));
    } else if (_elems.ndim == 2) {
      _collvec.push_back(Volumes<STORE>(_elems, _et, _mt, _bptr));
    }

  } else {
    // found a match - get the Collection that matched
    auto& coll = _collvec[imatch];

    // proceed to add the correct object type
    if (_elems.ndim == 0) {
      Points<STORE>& pts = std::get<Points<STORE>>(coll);
      pts.add_new(_elems, get_vdelta());
    } else if (_elems.ndim == 1) {
      Surfaces<STORE>& surf = std::get<Surfaces<STORE>>(coll);
      surf.add_new(_elems);
    } else if (_elems.ndim == 2) {
      Volumes<STORE>& vols = std::get<Volumes<STORE>>(coll);
      vols.add_new(_elems);
    }
  }
}

// Add elements - cells for hybrid calculation
void Simulation::add_hybrid(const std::vector<ElementPacket<float>> _elems,
                            std::shared_ptr<Body> _bptr) {

  // skip out early if nothing's here
  if (_elems.size() == 0) return;

#if defined(HOFORTRAN) || defined(HOCXX)
  // or if hybrid isn't turned on
  if (not hybr.is_active()) return;

  std::cout << "In Simulation::add_hybrid" << std::endl;
  std::cout << "  incoming vector has " << _elems.size() << " ElementPackets" << std::endl;

  // make sure we've got the right data
  assert((_elems.size() == 3 ||  _elems.size() == 5) && "Bad number of ElementPackets in add_hybrid");

  //_elems[0] is the volume elements - always add unique Collection to euler
  //_elems[1] is the wall boundaries
  //_elems[2] is the open boundaries
  euler.emplace_back(HOVolumes<STORE>(_elems[0], _elems[1], _elems[2], hybrid, fixed, _bptr));
  std::cout << "  euler now has " << euler.size() << " HOVolumes" << std::endl;

  // alternate way to assign the wall and open bc elements
  //euler.back().add_wall(_elems[1]);
  //euler.back().add_open(_elems[2]);

  // generate inlets and outlets, if available
  if (_elems.size() > 3) {
    euler.back().add_inlet(_elems[3]);
    euler.back().add_outlet(_elems[4]);
  }

  // tell the HOVolume to run its own conserve_iolet_volume routine to set/scale outlet rates
  euler.back().conserve_iolet_volume();
#endif
}

// add a new Body with the given name
void Simulation::add_body(std::shared_ptr<Body> _body) {
  bodies.emplace_back(_body);
  std::cout << "  added new body (" << _body->get_name() << "), now have " << bodies.size() << std::endl;
}

// return a Body pointer to the last Body in the array
std::shared_ptr<Body> Simulation::get_last_body() {
  std::shared_ptr<Body> bp;

  if (bodies.size() == 0) {
    std::cout << "  no last body found, creating (ground)" << std::endl;
    bp = std::make_shared<Body>();
    bp->set_name("ground");
    add_body(bp);
  } else {
    bp = bodies.back();
    std::cout << "  returning last body (" << bp->get_name() << ")" << std::endl;
  }

  return bp;
}

// return a Body pointer to the body matching the given name
std::shared_ptr<Body> Simulation::get_pointer_to_body(const std::string _name) {
  std::shared_ptr<Body> bp;

  for (auto &bptr : bodies) {
    if (_name.compare(bptr->get_name()) == 0) {
      std::cout << "  found body matching name (" << _name << ")" << std::endl;
      bp = bptr;
    }
  }

  // or ground if none match
  if (not bp) {
    std::cout << "  no body matching (" << _name << ") found, creating (ground)" << std::endl;
    bp = std::make_shared<Body>();
    bp->set_name("ground");
    add_body(bp);
  }

  return bp;
}

// check vs. step and time to see if simulation should pause/stop
bool Simulation::test_vs_stop() {
  bool should_stop = false;
  if (using_max_steps() and get_max_steps() == nstep) {
    std::cout << "Stopping at step " << get_max_steps() << std::endl;
    should_stop = true;
  }
  if (using_end_time() and get_end_time() <= time+0.5*dt){
    std::cout << "Stopping at time " << get_end_time() << std::endl;
    should_stop = true;
  }
  return should_stop;
}

// this is different because we have to trigger when last step is still running
bool Simulation::test_vs_stop_async() {
  bool should_stop = false;
  static bool already_reported = false;

  if (using_max_steps() and get_max_steps() == nstep+1) {
    if (not already_reported) {
      std::cout << std::endl << "Stopping at step " << get_max_steps() << std::endl;
      already_reported = true;
    }
    should_stop = true;
  }

  if (using_end_time() and get_end_time() >= time+0.5*dt
                       and get_end_time() <= time+1.5*dt) {
    if (not already_reported) {
      std::cout << std::endl << "Stopping at time " << get_end_time() << std::endl;
      already_reported = true;
    }
    should_stop = true;
  }

  // reset the toggle
  if (not should_stop) already_reported = false;

  return should_stop;
}

