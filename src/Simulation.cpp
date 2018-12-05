/*
 * Simulation.cpp - a class to control a 2d vortex particle sim
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Simulation.h"

#include <cassert>
#include <cmath>

// constructor
Simulation::Simulation()
  : re(100.0),
    dt(0.01),
    fs{0.0,0.0},
    vort(),
    bdry(),
    diff(),
    conv(),
    time(0.0),
    sim_is_initialized(false),
    step_has_started(false),
    step_is_finished(false)
  {}

// addresses for use in imgui
float* Simulation::addr_re() { return &re; }
float* Simulation::addr_dt() { return &dt; }
float* Simulation::addr_fs() { return fs; }

// getters
float Simulation::get_hnu() { return std::sqrt(dt/re); }
float Simulation::get_ips() { return diff.get_nom_sep_scaled() * get_hnu(); }
float Simulation::get_vdelta() { return diff.get_particle_overlap() * get_ips(); }
float Simulation::get_time() { return (float)time; }

// status
size_t Simulation::get_npanels() { return bdry.get_npanels(); }
size_t Simulation::get_nparts() { return vort.get_n(); }

// like a setter
void Simulation::set_re_for_ips(const float _ips) {
  re = std::pow(diff.get_nom_sep_scaled(), 2) * dt / pow(_ips, 2);
  diff.set_diffuse(false);
}

void Simulation::set_diffuse(const bool _do_diffuse) {
  diff.set_diffuse(_do_diffuse);
}

//void Simulation::set_amr(const bool _do_amr) {
//  diff.set_amr(_do_amr);
//  diff.set_diffuse(true);
//}

void Simulation::initGL(std::vector<float>& _projmat,
                        float*              _poscolor,
                        float*              _negcolor) {
  bdry.initGL(_projmat, _poscolor, _negcolor);
  vort.initGL(_projmat, _poscolor, _negcolor);
}
void Simulation::updateGL() {
  bdry.updateGL();
  vort.updateGL();
}
void Simulation::drawGL(std::vector<float>& _projmat,
                        float*              _poscolor,
                        float*              _negcolor) {
  if (step_is_finished) {
    bdry.drawGL(_projmat, _poscolor, _negcolor);
    vort.drawGL(_projmat, _poscolor, _negcolor);
  }
}

//
// main must indicate that panels should be made
//   because initGL and updateGL need to send data soon
//
void Simulation::init_bcs() {
  // are panels even made? do this first
  bdry.make_panels(get_ips());
}

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
  vort.reset();
  bdry.reset();
  sim_is_initialized = false;
  step_has_started = false;
  step_is_finished = false;
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

    // tell flow objects to update their values to the GPU
    updateGL();

    // set flag indicating that at least one step has been solved
    step_is_finished = true;
    step_has_started = false;

    return true;
  }

  // async call is not finished, do not try calling it again
  return false;
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
  std::cout << "taking step at t=" << time << " with n=" << vort.get_n() << std::endl;

  // we wind up using this a lot
  std::array<float,2> thisfs = reinterpret_cast<std::array<float,2>&>(fs);

  // are panels even made? do this first
  bdry.make_panels(get_ips());

  // for simplicity's sake, just run one full diffusion step here
  diff.step(dt, re, get_vdelta(), thisfs, vort, bdry);

  // operator splitting requires one half-step diffuse (use coefficients from previous step, if available)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);

  // advect with no diffusion (must update BEM strengths)
  //conv.advect_1st(dt, thisfs, vort, bdry);
  conv.advect_2nd(dt, thisfs, vort, bdry);

  // operator splitting requires another half-step diffuse (must compute new coefficients)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);

  // update strength for coloring purposes (eventually should be taken care of automatically)
  vort.update_max_str();

  // update dt and return
  time += (double)dt;
}

// set up the particles
void Simulation::add_particles(std::vector<float> _xysr) {

  if (_xysr.size() == 0) return;

  // make sure we're getting full particles
  assert(_xysr.size() % 4 == 0);

  // add the vdelta to each particle and pass it on
  for (size_t i=3; i<_xysr.size(); i+=4) {
    _xysr[i] = get_vdelta();
  }

  // also add to vorticity
  vort.add_new(_xysr);
}

// add geometry
void Simulation::add_boundary(bdryType _type, std::vector<float> _params) {
  if (_type == circle) {
    assert(_params.size() == 3);
    bdry.add(Circle<float>(_params[0], _params[1], _params[2]));
  }
}

