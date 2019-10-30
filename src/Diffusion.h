/*
 * Diffusion.h - a class for diffusion of strength from bodies to particles and among particles
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Core.h"
#include "Merge.h"
#include "Reflect.h"
#ifdef PLUGIN_AVRM
  #include "VRMadaptive.h"
#else
  #include "VRM.h"
#endif
#include "BEM.h"

#include "json/json.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>


//
// One step of convection of elements, allowing for boundary conditions
//
// templatized on 'S'torage, 'A'ccumulator, and 'I'ndex types
//
template <class S, class A, class I>
class Diffusion {
public:
  Diffusion()
    : vrm(),
      core_func(gaussian),
      is_inviscid(false),
      adaptive_radii(false),
      nom_sep_scaled(std::sqrt(8.0)),
      particle_overlap(1.5),
      shed_before_diffuse(false)
    {}

  const bool get_diffuse() const { return !is_inviscid; }
  void set_diffuse(const bool _do_diffuse) { is_inviscid = not _do_diffuse; }
  void set_amr(const bool _do_amr) { adaptive_radii = _do_amr; }
  S get_nom_sep_scaled() const { return nom_sep_scaled; }
  S get_nom_sep() { return nom_sep_scaled * vrm.get_hnu(); }
  S get_particle_overlap() const { return particle_overlap; }

  // take a full diffusion step
  void step(const double,
            const double,
            const S,
            const S,
            const std::array<double,Dimensions>&,
            std::vector<Collection>&,
            std::vector<Collection>&,
            BEM<S,I>& _bem);

  // read/write parameters
  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

private:
  // the VRM algorithm, template params are storage, solver, max moments
  VRM<S,double,2> vrm;

  // other necessary variables
  CoreType core_func;

  // toggle inviscid
  bool is_inviscid;

  // toggle adaptive particle sizes
  bool adaptive_radii;

  // nominal separation normalized by h_nu
  S nom_sep_scaled;

  // particle core size is nominal separation times this
  S particle_overlap;

  // method 1 (true) is to shed *at* the boundary, VRM those particles, then push out
  // method 2 (false) is to VRM, push out, *then* generate new particles at the correct distance
  bool shed_before_diffuse;
};


//
// take a diffusion step
//
template <class S, class A, class I>
void Diffusion<S,A,I>::step(const double                _time,
                            const double                _dt,
                            const S                     _re,
                            const S                     _vdelta,
                            const std::array<double,Dimensions>& _fs,
                            std::vector<Collection>&    _vort,
                            std::vector<Collection>&    _bdry,
                            BEM<S,I>&                   _bem) {

  if (is_inviscid) return;

  std::cout << "Inside Diffusion::step with dt=" << _dt << std::endl;

  // ensure that we have a current h_nu
  vrm.set_hnu(std::sqrt(_dt/_re));

#ifdef PLUGIN_AVRM
  // ensure that it knows to allow or disallow adaptive radii
  vrm.set_adaptive_radii(adaptive_radii);
#endif

  //
  // always re-run the BEM calculation before shedding
  //
  solve_bem<S,A,I>(_time, _fs, get_nom_sep(), _vort, _bdry, _bem);

  //
  // important for augmented BEM: reset the circulation counter
  //
  for (auto &coll : _bdry) {
    // run this step if the collection is Surfaces
    if (std::holds_alternative<Surfaces<S>>(coll)) {
      Surfaces<S>& surf = std::get<Surfaces<S>>(coll);
      // zeros reabsorbed circ
      surf.reset_augmentation_vars();
    }
  }

  //
  // generate particles at boundary surfaces
  //
  if (shed_before_diffuse) {
    for (auto &coll : _bdry) {

      // run this step if the collection is Surfaces
      if (std::holds_alternative<Surfaces<S>>(coll)) {
        Surfaces<S>& surf = std::get<Surfaces<S>>(coll);

        // generate particles just above the surface
        std::vector<S> new_pts = surf.represent_as_particles(0.01*(S)vrm.get_hnu(), _vdelta);

        // add those particles to the main particle list
        if (_vort.size() == 0) {
          // no collections yet? make a new collection
          _vort.push_back(Points<S>(new_pts, active, lagrangian, nullptr));      // vortons
        } else {
          // HACK - add all particles to first collection
          auto& coll = _vort.back();
          // only proceed if the last collection is Points
          if (std::holds_alternative<Points<S>>(coll)) {
            Points<S>& pts = std::get<Points<S>>(coll);
            pts.add_new(new_pts);
          }
        }
      }

      // Kutta points and lifting lines can generate points here
    }
  }

  //
  // diffuse strength among existing particles
  //
  // loop over active vorticity
  for (auto &coll : _vort) {

    // if no strength, skip
    if (std::visit([=](auto& elem) { return elem.is_inert(); }, coll)) continue;

    // run this step if the collection is Points
    if (std::holds_alternative<Points<S>>(coll)) {

      Points<S>& pts = std::get<Points<S>>(coll);
      std::cout << "    computing diffusion among " << pts.get_n() << " particles" << std::endl;

      // vectors are not passed as const, because they may be extended with new particles
      // this call also applies the changes, though we may want to save any changes into another
      //   vector of derivatives to be applied later
      vrm.diffuse_all(pts.get_pos(),
                      pts.get_str(),
                      pts.get_rad(),
                      core_func,
                      particle_overlap);

      // resize the rest of the arrays
      pts.resize(pts.get_rad().size());
    }
  }


  //
  // reflect interior particles to exterior because VRM only works in free space
  //
  (void) reflect_interior<S>(_bdry, _vort);


  //
  // merge any close particles to clean up potentially-dense areas
  // when this step runs after clear_inner_panp2, drag drops a little bit
  //
  (void) merge_operation<S>(_vort, particle_overlap, 0.3, adaptive_radii);


  //
  // clean up by removing the innermost layer - the one that will be represented by boundary strengths
  //
  if (shed_before_diffuse) {
    // use method which trims circulations under the threshold
    (void) clear_inner_layer<S>(0, _bdry, _vort, 0.0, _vdelta/particle_overlap);
  } else {
    // use method which simply pushes all still-active particles to be at or above a threshold distance
    // cutoff is a multiple of ips (these are the last two arguments)
    (void) clear_inner_layer<S>(1, _bdry, _vort, 1.0/std::sqrt(2.0*M_PI), _vdelta/particle_overlap);
  }


  //
  // generate particles above boundary surfaces
  //
  if (not shed_before_diffuse) {
    for (auto &coll : _bdry) {

      // run this step if the collection is Surfaces
      if (std::holds_alternative<Surfaces<S>>(coll)) {
        Surfaces<S>& surf = std::get<Surfaces<S>>(coll);

        // generate particles above the surface at the centroid of one step of
        //   diffusion from a flat plate
        std::vector<S> new_pts = surf.represent_as_particles(vrm.get_hnu()*std::sqrt(4.0/M_PI), _vdelta);

        // add those particles to the main particle list
        if (_vort.size() == 0) {
          // no collections yet? make a new collection
          _vort.push_back(Points<S>(new_pts, active, lagrangian, nullptr));      // vortons
        } else {
          // HACK - add all particles to first collection
          auto& coll = _vort.back();
          // only proceed if the last collection is Points
          if (std::holds_alternative<Points<S>>(coll)) {
            Points<S>& pts = std::get<Points<S>>(coll);
            pts.add_new(new_pts);
          }
        }
      }

      // Kutta points and lifting lines can generate points here
    }
  }

  //
  // merge again if clear did any work
  //
  if (_bdry.size() > 0) merge_operation<S>(_vort, particle_overlap, 0.3, adaptive_radii);


  // now is a fine time to reset the max active/particle strength
  for (auto &coll : _vort) {
    std::visit([=](auto& elem) { elem.update_max_str(); }, coll);
  }
}

//
// read/write parameters to json
//

// read "simparams" json object
template <class S, class A, class I>
void Diffusion<S,A,I>::from_json(const nlohmann::json j) {

  if (j.find("viscous") != j.end()) {
    std::string viscous = j["viscous"];
    if (viscous == "vrm") { set_diffuse(true); }
    else { set_diffuse(false); }
    std::cout << "  setting is_viscous= " << get_diffuse() << std::endl;
  }

  vrm.from_json(j);

#ifdef PLUGIN_AVRM
  // set adaptive-VRM-specific settings
  if (j.find("adaptiveSize") != j.end()) {
    // for now, just enable it, don't set parameters
    set_amr(true);
    std::cout << "  enabling amr" << std::endl;
  }
#endif
}

// create and write a json object for all diffusion parameters
template <class S, class A, class I>
void Diffusion<S,A,I>::add_to_json(nlohmann::json& j) const {
  //nlohmann::json j;

  j["viscous"] = get_diffuse() ? "vrm" : "none";
#ifdef PLUGIN_AVRM
  j["adaptiveSize"] = adaptive_radii;
#endif

  // eventually write other parameters
  //j["overlap"] = particle_overlap;
  //j["core"] = core_func;

  // VRM will write "VRM" and "AMR" parameters
  vrm.add_to_json(j);
}

