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
#include "VRM.h"
#include "BEM.h"

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
      particle_overlap(1.5)
    {}

  const bool get_diffuse() { return !is_inviscid; }
  void set_diffuse(const bool _do_diffuse) { is_inviscid = not _do_diffuse; }
  //void set_amr(const bool _do_amr) { adaptive_radii = _do_amr; }
  S get_nom_sep_scaled() const { return nom_sep_scaled; }
  S get_particle_overlap() const { return particle_overlap; }

  void step(const double,
            const S,
            const S,
            const std::array<double,2>&,
            std::vector<Collection>& _vort,
            std::vector<Collection>& _bdry,
            BEM<S,I>& _bem);

private:
  // the VRM algorithm, template params are storage, compute, max moments
  // note that NNLS needs doubles for its compute type or else it will fail
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
};


//
// take a diffusion step
//
template <class S, class A, class I>
void Diffusion<S,A,I>::step(const double                _dt,
                            const S                     _re,
                            const S                     _vdelta,
                            const std::array<double,2>& _fs,
                            std::vector<Collection>&    _vort,
                            std::vector<Collection>&    _bdry,
                            BEM<S,I>&                   _bem) {

  if (is_inviscid) return;

  std::cout << "Inside Diffusion::step with dt=" << _dt << std::endl;

  // always re-run the BEM calculation before shedding
  solve_bem<S,A,I>(_fs, _vort, _bdry, _bem);

  //
  // generate particles at boundary surfaces
  //

  for (auto &coll : _bdry) {

    // run this step if the collection is Surfaces
    if (std::holds_alternative<Surfaces<S>>(coll)) {
      Surfaces<S>& surf = std::get<Surfaces<S>>(coll);

      // generate particles just above the surface
      std::vector<S> new_pts = surf.represent_as_particles(0.0001*(S)_dt, _vdelta);

      // add those particles to the main particle list
      if (_vort.size() == 0) {
        // no collections yet? make a new collection
        _vort.push_back(Points<S>(new_pts, active, lagrangian));      // vortons
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

  //
  // diffuse strength among existing particles
  //

  // ensure that we have a current h_nu
  vrm.set_hnu(std::sqrt(_dt/_re));

  // ensure that it knows to allow or disallow adaptive radii
  //vrm.set_adaptive_radii(adaptive_radii);

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
      pts.resize(pts.get_str().size());
    }
  }

  //
  // reflect interior particles to exterior because VRM only works in free space
  //

  //ReflectVisitor<S> rvisitor;
  // this should only function when _vort is Points and _bdry is Surfaces
  for (auto &targ : _vort) {
    if (std::holds_alternative<Points<S>>(targ)) {
      Points<S>& pts = std::get<Points<S>>(targ);

      for (auto &src : _bdry) {
        if (std::holds_alternative<Surfaces<S>>(src)) {
          Surfaces<S>& surf = std::get<Surfaces<S>>(src);

          // call the specific panels-affect-points routine
          (void) reflect_panp2<S>(surf, pts);
        }
      }
    }
  }

  //
  // merge any close particles to clean up potentially-dense areas
  //

  for (auto &coll : _vort) {

    // if inert, no need to merge (keep all tracer particles)
    if (std::visit([=](auto& elem) { return elem.is_inert(); }, coll)) continue;

    // run this step if the collection is Points
    if (std::holds_alternative<Points<S>>(coll)) {

      Points<S>& pts = std::get<Points<S>>(coll);

      //std::cout << "    merging among " << pts.get_n() << " particles" << std::endl;
      // last two arguments are: relative distance, allow variable core radii
      (void) merge_close_particles<S>(pts.get_pos(),
                                      pts.get_str(),
                                      pts.get_rad(),
                                      particle_overlap,
                                      0.3,
                                      adaptive_radii);

      // resize the rest of the arrays
      pts.resize(pts.get_str().size());
    }
  }

  //
  // clean up by removing the innermost layer - the one that will be represented by boundary strengths
  //

  // may need to do this multiple times to clear out concave zones!
  // this should only function when _vort is Points and _bdry is Surfaces
  for (auto &targ : _vort) {
    if (std::holds_alternative<Points<S>>(targ)) {
      Points<S>& pts = std::get<Points<S>>(targ);

      for (auto &src : _bdry) {
        if (std::holds_alternative<Surfaces<S>>(src)) {
          Surfaces<S>& surf = std::get<Surfaces<S>>(src);

          // call the specific panels-affect-points routine
          (void) clear_inner_panp2<S>(surf, pts, _vdelta/particle_overlap);
        }
      }
    }
  }

  // now is a fine time to reset the max active/particle strength
  for (auto &coll : _vort) {
    std::visit([=](auto& elem) { elem.update_max_str(); }, coll);
  }
}

