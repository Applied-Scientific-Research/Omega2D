/*
 * Diffusion.h - a class for diffusion of strength from bodies to particles and among particles
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Written by Mark J Stock <markjstock@gmail.com>
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
#include "PSE.h"
#include "RVM.h"
#include "CoreSpread.h"
#include "BEM.h"
#include "GuiHelper.h"

#include "json/json.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>


// eventually support multiple diffusion types, both inter-particle and boundary-to-fluid
enum PartDiffuseType {
  pd_none,
  pd_core,
  pd_rvm,
  pd_pse,
  pd_vrm,
  pd_avrm
  };
//enum BdryDiffuseType { none, single, vrm, hybrid };


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
      pse(),
      rvm(),
      coresp(),
      h_nu(0.1),
      core_func(gaussian),
      is_inviscid(false),
      pd_type(pd_vrm),
      adaptive_radii(false),
      nom_sep_scaled(std::sqrt(8.0)),
      particle_overlap(1.5),
      merge_thresh(0.4),
      shed_before_diffuse(true)
    {}

  void set_diffuse(const bool _do_diffuse) { is_inviscid = not _do_diffuse; }
  void set_amr(const bool _do_amr);
  const bool get_diffuse() const { return !is_inviscid; }
  const bool get_amr() const { return adaptive_radii; }
  S get_nom_sep_scaled() const { return nom_sep_scaled; }
  S get_nom_sep(const S _hnu) { return nom_sep_scaled * _hnu; }
  S get_nom_sep(const double _dt, const S _re) { return nom_sep_scaled * std::sqrt(_dt/_re); }
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

#ifdef USE_IMGUI
  void draw_advanced();
#endif

  // read/write parameters
  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

private:
  // the VRM algorithm, template params are storage, solver, max moments
  VRM<S,double,2> vrm;

  // the particle-strength-exchange method
  PSE<S,A> pse;

  // the random vortex method
  RVM<S> rvm;

  // the core spreading method
  CoreSpread<S> coresp;

  // other necessary variables
  S h_nu;
  CoreType core_func;

  // toggle inviscid
  bool is_inviscid;

  // diffusion method
  PartDiffuseType pd_type;

  // toggle adaptive particle sizes
  bool adaptive_radii;

  // nominal separation normalized by h_nu
  S nom_sep_scaled;

  // particle core size is nominal separation times this
  S particle_overlap;

  // merge aggressivity
  S merge_thresh;

  // method 1 (true) is to shed *at* the boundary, VRM those particles, then push out
  // method 2 (false) is to VRM, push out, *then* generate new particles at the correct distance
  bool shed_before_diffuse;
};

//
template <class S, class A, class I>
void Diffusion<S,A,I>::set_amr(const bool _do_amr) {
  adaptive_radii = _do_amr;
  if (_do_amr) set_diffuse(true);
}

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

  // don't let part diffusion type change during execution
  const PartDiffuseType curr_pd_type = pd_type;

  if (is_inviscid or curr_pd_type==pd_none) return;

  // some methods require different merging properties
  if (curr_pd_type==pd_rvm) merge_thresh = 0.0;
  if (curr_pd_type==pd_core) merge_thresh = 0.02;
  else merge_thresh = 0.2;

  std::cout << "Inside Diffusion::step with dt=" << _dt << std::endl;

  // ensure that we have a current h_nu
  assert((S)_re != 0); // Can't divide by 0
  h_nu = (S)std::sqrt(_dt/_re);

#ifdef PLUGIN_AVRM
  // ensure that it knows to allow or disallow adaptive radii
  if (curr_pd_type == pd_vrm) vrm.set_adaptive_radii(adaptive_radii);
#endif

  //
  // always re-run the BEM calculation before shedding
  //
  // first push away particles inside or too close to the body
  assert(M_PI != 0); // Can't divide by 0
  clear_inner_layer<S>(1, _bdry, _vort, 0.5/std::sqrt(2.0*M_PI), get_nom_sep(h_nu));
  solve_bem<S,A,I>(_time, _fs, _vort, _bdry, _bem);

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
        std::vector<S> new_pts = surf.represent_as_particles(0.01*(S)h_nu, _vdelta);

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

      if (curr_pd_type==pd_vrm) {
        // vectors are not passed as const, because they may be extended with new particles
        // this call also applies the changes, though we may want to save any changes into another
        //   vector of derivatives to be applied later
        vrm.diffuse_all(pts.get_pos(),
                        pts.get_str(),
                        pts.get_rad(),
                        h_nu, core_func,
                        particle_overlap);

        // resize the rest of the arrays
        pts.resize(pts.get_rad().size());

      } else if (curr_pd_type==pd_pse) {
        // vectors are not passed as const, because they may be extended with new particles
        // this call also applies the changes, though we may want to save any changes into another
        //   vector of derivatives to be applied later
        pse.diffuse_all(pts.get_pos(),
                        pts.get_str(),
                        pts.get_rad(),
                        h_nu, core_func,
                        particle_overlap);

        // resize the rest of the arrays
        pts.resize(pts.get_rad().size());

      } else if (curr_pd_type==pd_core) {
        // core-spreading only changes the radii, nothing else
        coresp.diffuse_all(pts.get_pos(),
                        pts.get_str(),
                        pts.get_rad(),
                        h_nu, core_func);

      } else if (curr_pd_type==pd_rvm) {
        // RVM only changes the positions, nothing else
        rvm.diffuse_all(pts.get_pos(),
                        pts.get_str(),
                        pts.get_rad(),
                        h_nu);
      }
    }
  }


  //
  // reflect interior particles to exterior because VRM only works in free space
  //
  (void) reflect_interior<S>(_bdry, _vort);


  //
  // merge any close particles to clean up potentially-dense areas
  //
  if (curr_pd_type != pd_rvm) merge_operation<S>(_vort, particle_overlap, merge_thresh, adaptive_radii);


  //
  // clean up by removing the innermost layer - the one that will be represented by boundary strengths
  //
  if (shed_before_diffuse) {
    // use method which trims circulations under the threshold
    //(void) clear_inner_layer<S>(0, _bdry, _vort, 0.0, _vdelta/particle_overlap); // THIS IS BAD
    (void) clear_inner_layer<S>(1, _bdry, _vort, 0.5/std::sqrt(2.0*M_PI), _vdelta/particle_overlap);
  } else {
    // use method which simply pushes all still-active particles to be at or above a threshold distance
    // cutoff is a multiple of ips (these are the last two arguments)
    (void) clear_inner_layer<S>(1, _bdry, _vort, 0.5/std::sqrt(2.0*M_PI), _vdelta/particle_overlap);
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
        std::vector<S> new_pts = surf.represent_as_particles(h_nu*std::sqrt(4.0/M_PI), _vdelta);

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
  if (_bdry.size() > 0 and curr_pd_type != pd_rvm) merge_operation<S>(_vort, particle_overlap, merge_thresh, adaptive_radii);


  // now is a fine time to reset the max active/particle strength
  for (auto &coll : _vort) {
    std::visit([=](auto& elem) { elem.update_max_str(); }, coll);
  }
}

#ifdef USE_IMGUI
//
// draw advanced options parts of the GUI
//
template <class S, class A, class I>
void Diffusion<S,A,I>::draw_advanced() {

  ImGui::Separator();
  ImGui::Spacing();
  ImGui::Text("Diffusion settings");


  // select diffusion type among available
  int diff_item = 3;		// default is VRM
  if (pd_type==pd_core) diff_item = 0;
  else if (pd_type==pd_rvm) diff_item = 1;
  else if (pd_type==pd_pse) diff_item = 2;
  else if (pd_type==pd_vrm) diff_item = 3;

  // draw the selector box, with the current one selected
  const char* diff_items[] = { "Core-spreading", "Random walk", "Particle strength exchange", "VRM solution" };
  ImGui::PushItemWidth(240);
  ImGui::Combo("Select diffusion method", &diff_item, diff_items, 4);
  ImGui::PopItemWidth();
  switch(diff_item) {
      case 0: pd_type = pd_core; break;
      case 1: pd_type = pd_rvm; break;
      case 2: pd_type = pd_pse; break;
      case 3: pd_type = pd_vrm; break;
  } // end switch


  // now, present options, depending on the diffusion type
  if (pd_type == pd_vrm) {

  bool relative_thresh = vrm.get_relative();
  ImGui::Checkbox("Thresholds are relative to strongest particle", &relative_thresh);
  ImGui::SameLine();
  ShowHelpMarker("If unchecked, the thresholds defined here are absolute and unscaled to the strongest particle.");
  vrm.set_relative(relative_thresh);

  ImGui::PushItemWidth(-270);
  float ignore_thresh = std::log10(vrm.get_ignore());
  ImGui::SliderFloat("Threshold to ignore", &ignore_thresh, -12, 0, "%.1f");
  ImGui::SameLine();
  ShowHelpMarker("During diffusion, ignore any particles with strength magnitude less than this power of ten threshold.");
  vrm.set_ignore(std::pow(10.f,ignore_thresh));
  ImGui::PopItemWidth();

#ifdef PLUGIN_SIMPLEX
  // bool toggle for NNLS vs. Simplex
  bool use_simplex = vrm.get_simplex();
  ImGui::Checkbox("VRM uses Simplex solver", &use_simplex);
  ImGui::SameLine();
  ShowHelpMarker("Use the proprietary Simplex solver for overdetermined systems. If unchecked, the Vorticity Redistribution Method uses a Non-Negative Least Squares solver from Eigen.");
  vrm.set_simplex(use_simplex);
#endif

#ifdef PLUGIN_AVRM
  // show the toggle for AMR
  bool use_amr = get_amr();
  ImGui::Checkbox("Allow adaptive resolution", &use_amr);
  ImGui::SameLine();
  ShowHelpMarker("Particle sizes will adapt as required to maintain resolution during the diffusion calculation. If unchecked, all particles will stay the same size.");
  set_amr(use_amr);

  if (use_amr) vrm.draw_advanced();
#endif

  } else if (pd_type == pd_pse) {

    bool relative_thresh = pse.get_relative();
    ImGui::Checkbox("Thresholds are relative to strongest particle", &relative_thresh);
    ImGui::SameLine();
    ShowHelpMarker("If unchecked, the thresholds defined here are absolute and unscaled to the strongest particle.");
    pse.set_relative(relative_thresh);

    ImGui::PushItemWidth(-270);
    float ignore_thresh = std::log10(pse.get_ignore());
    ImGui::SliderFloat("Threshold to ignore", &ignore_thresh, -12, 0, "%.1f");
    ImGui::SameLine();
    ShowHelpMarker("During diffusion, ignore any particles with strength magnitude less than this power of ten threshold.");
    pse.set_ignore(std::pow(10.f,ignore_thresh));
    ImGui::PopItemWidth();

    bool pse_uses_vols = pse.get_volumes();
    ImGui::Checkbox("Use volumes to weigh exchange", &pse_uses_vols);
    ImGui::SameLine();
    ShowHelpMarker("If unchecked, PSE will become inaccurate when particles are jumbled.");
    pse.set_volumes(pse_uses_vols);

  } else if (pd_type == pd_core) {
    // no special parameters

  } else if (pd_type == pd_rvm) {
    // no special parameters

  }
}
#endif

//
// read/write parameters to json
//

// read "simparams" json object
template <class S, class A, class I>
void Diffusion<S,A,I>::from_json(const nlohmann::json j) {

  if (j.find("viscous") != j.end()) {
    std::string viscous = j["viscous"];
    if (viscous == "vrm") {
      set_diffuse(true);
      pd_type = pd_vrm;
    } else if (viscous == "pse") {
      set_diffuse(true);
      pd_type = pd_pse;
    } else if (viscous == "random") {
      set_diffuse(true);
      pd_type = pd_rvm;
    } else if (viscous == "corespread") {
      set_diffuse(true);
      pd_type = pd_core;
    } else {
      // "none" or unsupported
      set_diffuse(false);
    }
  } else {
    set_diffuse(false);
  }
  std::cout << "  setting is_viscous= " << get_diffuse() << std::endl;

  // regardless, load some settings as they were
  vrm.from_json(j);
  pse.from_json(j);

#ifdef PLUGIN_AVRM
  // set adaptive-VRM-specific settings
  if (j.find("adaptiveSize") != j.end()) {
    if (j["adaptiveSize"]) {
      set_amr(true);
      std::cout << "  enabling amr" << std::endl;
    }
  }
#endif
}

// create and write a json object for all diffusion parameters
template <class S, class A, class I>
void Diffusion<S,A,I>::add_to_json(nlohmann::json& j) const {
  //nlohmann::json j;

  if (not get_diffuse()) {
    j["viscous"] = "none";
  } else if (pd_type == pd_core) {
    j["viscous"] = "corespread";
  } else if (pd_type == pd_rvm) {
    j["viscous"] = "random";
  } else if (pd_type == pd_pse) {
    j["viscous"] = "pse";
  } else if (pd_type == pd_vrm) {
    j["viscous"] = "vrm";
  }

#ifdef PLUGIN_AVRM
  j["adaptiveSize"] = adaptive_radii;
#endif

  // eventually write other parameters
  //j["overlap"] = particle_overlap;
  //j["core"] = core_func;

  // VRM always writes "VRM" and "AMR" parameters
  vrm.add_to_json(j);
  // PSE always writes its parameters too
  pse.add_to_json(j);
}

