/*
 * Simulation.cpp - a class to control a 2d vortex particle sim
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Simulation.h"
#include "Points.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <variant>

// constructor
Simulation::Simulation()
  : re(100.0),
    dt(0.01),
    fs{0.0,0.0},
    vort(),
    vort2(),
    bdry(),
    bdry2(),
    fldpt(),
    bem(),
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
size_t Simulation::get_npanels() { //return bdry.get_npanels(); }
  size_t n = 0;
  for (auto &coll: bdry2) {
    //std::visit([&n](auto& elem) { n += elem.get_npanels(); }, coll);
    // only proceed if the last collection is Surfaces
    if (std::holds_alternative<Surfaces<float>>(coll)) {
      Surfaces<float>& surf = std::get<Surfaces<float>>(coll);
      n += surf.get_npanels();
    }
  }
  return n;
}

size_t Simulation::get_nparts() {
  size_t n = 0;
  for (auto &coll: vort2) {
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

#ifdef USE_GL
void Simulation::initGL(std::vector<float>& _projmat,
                        float*              _poscolor,
                        float*              _negcolor,
                        float*              _defcolor) {
  //bdry.initGL(_projmat, _poscolor, _negcolor);
  //vort.initGL(_projmat, _poscolor, _negcolor);
  for (auto &coll : vort2) {
    std::visit([=, &_projmat](auto& elem) { elem.initGL(_projmat, _poscolor, _negcolor, _defcolor); }, coll);
  }
  for (auto &coll : bdry2) {
    std::visit([=, &_projmat](auto& elem) { elem.initGL(_projmat, _poscolor, _negcolor, _defcolor); }, coll);
  }
  for (auto &coll : fldpt) {
    std::visit([=, &_projmat](auto& elem) { elem.initGL(_projmat, _poscolor, _negcolor, _defcolor); }, coll);
  }
}

void Simulation::updateGL() {
  bdry.updateGL();
  //vort.updateGL();
  for (auto &coll : vort2) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
  for (auto &coll : bdry2) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
  for (auto &coll : fldpt) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
}

void Simulation::drawGL(std::vector<float>& _projmat,
                        float*              _poscolor,
                        float*              _negcolor,
                        float*              _defcolor,
                        float               _tracersize) {
  if (step_is_finished) {
    //bdry.drawGL(_projmat, _poscolor, _negcolor);
    //vort.drawGL(_projmat, _poscolor, _negcolor);
    const float tracersz = get_ips() * _tracersize;
    for (auto &coll : vort2) {
      std::visit([=, &_projmat](auto& elem) { elem.drawGL(_projmat, _poscolor, _negcolor, _defcolor, tracersz); }, coll);
    }
    for (auto &coll : bdry2) {
      std::visit([=, &_projmat](auto& elem) { elem.drawGL(_projmat, _poscolor, _negcolor, _defcolor, tracersz); }, coll);
    }
    for (auto &coll : fldpt) {
      std::visit([=, &_projmat](auto& elem) { elem.drawGL(_projmat, _poscolor, _negcolor, _defcolor, tracersz); }, coll);
    }
  }
}
#endif

//
// main must indicate that panels should be made
//   because initGL and updateGL need to send data soon
//
void Simulation::init_bcs() {
  // are panels even made? do this first
  bdry.make_panels(get_ips());
  //bdry2.make_panels(get_ips());
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
  vort2.clear();
  bdry.reset();
  bdry2.clear();
  fldpt.clear();
  bem.reset();
  sim_is_initialized = false;
  step_has_started = false;
  step_is_finished = false;
}

//
// Check all aspects of the simulation for conditions that should stop the run
//
std::string Simulation::check_simulation(const size_t _nff, const size_t _nbf) {
  std::string retstr;

  // Check for no bodies and no particles
  if (_nbf == 0 and vort.get_n() == 0) {
    retstr.append("No flow features and no bodies. Add one or both, reset, and run.\n");
  }

  // Check for a body and no particles
  if (_nbf > 0 and vort.get_n() == 0) {

    // AND no freestream
    if (fs[0]*fs[0]+fs[1]*fs[1] < std::numeric_limits<float>::epsilon()) {
      retstr.append("No flow features and zero freestream speed - try adding one or both.\n");
      return retstr;
    }

    // AND no viscosity
    if (not diff.get_diffuse()) {
      retstr.append("You have a solid body, but no diffusion. It will not shed vorticity. Turn on viscosity or add a flow feature, reset, and run.\n");
    }
  }

  // Check for conditions that lead to loss of accuracy
  // like vorticity-based Courant number
  //static bool ignore_warning = false;
  //float max_elong = 0.0;
  //for (auto &coll: vort) {
    //std::visit([&](auto& elem) { max_elong = std::max(max_elong, elem.get_max_elong(); }, coll);
  //}
  //if (max_elong > 2.0) retstr.append("Elongation threshold exceeded! Reset and reduce the time step size.\n");

  return retstr;
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
  std::cout << std::endl << "Taking step at t=" << time << " with n=" << vort.get_n() << std::endl;

  // we wind up using this a lot
  std::array<double,2> thisfs = {fs[0], fs[1]};

  // are panels even made? do this first
  //bdry.make_panels(get_ips());

  // for simplicity's sake, just run one full diffusion step here
  //diff.step(dt, re, get_vdelta(), thisfs, vort, bdry);
  //diff.step(dt, re, get_vdelta(), thisfs, vort2, bdry2, bem);

  // operator splitting requires one half-step diffuse (use coefficients from previous step, if available)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);

  // advect with no diffusion (must update BEM strengths)
  //std::cout << std::endl << "STARTING OLD ARCH" << std::endl;
  //conv.advect_1st(dt, thisfs, vort, bdry);
  //conv.advect_2nd(dt, thisfs, vort, bdry);

  // advect using new architecture
  //std::cout << std::endl << "STARTING NEW ARCH" << std::endl;
  conv.advect_1st(dt, thisfs, vort2, bdry2, fldpt, bem);
  //conv.advect_2nd(dt, thisfs, vort2, bdry2, fldpt, bem);

  // operator splitting requires another half-step diffuse (must compute new coefficients)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);

  // update strengths for coloring purposes (eventually should be taken care of automatically)
  //vort.update_max_str();

  // update dt and return
  time += (double)dt;
}

// add some vortex particles to both old and new arch
void Simulation::add_particles(std::vector<float> _xysr) {

  if (_xysr.size() == 0) return;

  // make sure we're getting full particles
  assert(_xysr.size() % 4 == 0);

  // add the vdelta to each particle and pass it on
  const float thisvd = get_vdelta();
  for (size_t i=3; i<_xysr.size(); i+=4) {
    _xysr[i] = thisvd;
  }

  // make a copy of the vector so that the new arch can process it, too
  std::vector<float> incopy = _xysr;

  // also add to vorticity
  vort.add_new(_xysr);

  // add to new archtecture's vorticity (vort2)
  // if no collections exist
  if (vort2.size() == 0) {
    // make a new collection
    vort2.push_back(Points<float>(incopy, active, lagrangian));      // vortons
  } else {
    // THIS MUST USE A VISITOR
    // HACK - add all particles to first collection
    //std::visit([&](auto& elem) { elem.add_new(incopy); }, vort2.back());
    auto& coll = vort2.back();
    // only proceed if the last collection is Points
    if (std::holds_alternative<Points<float>>(coll)) {
      Points<float>& pts = std::get<Points<float>>(coll);
      pts.add_new(incopy);
    }
  }
}

// add some tracer particles to new arch
void Simulation::add_tracers(std::vector<float> _xy) {

  if (_xy.size() == 0) return;

  // make sure we're getting full points
  assert(_xy.size() % 2 == 0);

  // add to new archtecture

  // if no collections exist
  if (fldpt.size() == 0) {
    // make a new collection
    fldpt.push_back(Points<float>(_xy, inert, lagrangian));      // vortons

  } else {
    // THIS MUST USE A VISITOR
    // HACK - add all particles to first collection
    //std::visit([&](auto& elem) { elem.add_new(_xy); }, fldpt.back());
    auto& coll = fldpt.back();
    // eventually we will want to check every collection for matching element and move types
    // only proceed if the collection is Points
    if (std::holds_alternative<Points<float>>(coll)) {
      Points<float>& pts = std::get<Points<float>>(coll);
      pts.add_new(_xy);
    }
  }
}

// add geometry
void Simulation::add_boundary(bdryType _type, std::vector<float> _params) {
  if (_type == circle) {
    assert(_params.size() == 3);
    bdry.add(Circle<float>(_params[0], _params[1], _params[2]));
  }
}

void Simulation::add_boundary(ElementPacket<float> _geom) {

  // if no collections exist
  if (bdry2.size() == 0) {
    // make a new collection - assume BEM panels
    bdry2.push_back(Surfaces<float>(_geom.x,
                                    _geom.idx,
                                    _geom.val,
                                    0, reactive, fixed));
  } else {

    auto& coll = bdry2.back();
    // only proceed if the last collection is Surfaces
    // eventually check each collection for a element and movement type match (i.e. reactive and fixed)
    if (std::holds_alternative<Surfaces<float>>(coll)) {
      Surfaces<float>& surf = std::get<Surfaces<float>>(coll);
      surf.add_new(_geom.x,
                   _geom.idx,
                   _geom.val);
    }
  }
}

