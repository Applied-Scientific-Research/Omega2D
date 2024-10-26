/*
 * BEMHelper.h - non-class to coordinate solving the BEM problem
 *
 * (c)2019-21 Applied Scientific Research, Inc.
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
#include "Collection.h"
#include "CollectionHelper.h"
#include "Influence.h"
#include "Coefficients.h"
#include "RHS.h"
#include "BEM.h"
#include "ResultsType.h"
#include "ExecEnv.h"

#include <cstdlib>
#include <iostream>
#include <array>
#include <vector>
#include <cassert>


//
// helper function to solve BEM equations on given state
//
template <class S, class A, class I>
void solve_bem(const double                         _time,
               const std::array<double,Dimensions>& _fs,
               std::vector<Collection>&             _vort,
               std::vector<Collection>&             _bdry,
               BEM<S,I>&                            _bem) {

  // no unknowns? no problem.
  if (_bdry.size() == 0) return;

  // save the simulation time from the last time we entered this function
  static double last_time = -99.9;

  // if this is the first time through after a reset, recalculate the row indices
  if (not _bem.is_A_current()) {
    I rowcnt = 0;
    // loop over boundary collections
    for (auto &targ : _bdry) {
      std::visit([=](auto& elem) { return elem.set_first_row(rowcnt); }, targ);
      rowcnt = std::visit([=](auto& elem) { return elem.get_next_row(); }, targ);
    }
  }

  // add vortex and source strengths to account for rotating bodies - unless we augment!
  //for (auto &src : _bdry) {
  //  std::visit([=](auto& elem) { elem.add_unit_rot_strengths(); }, src);
  //}

  // need this for dispatching velocity influence calls, template param is accumulator type,
  //   member variable is default execution environment
#ifdef USE_VC
  InfluenceVisitor<S,A> ivisitor = {ResultsType(velonly), ExecEnv(true, direct, cpu_vc)};
#else
  InfluenceVisitor<S,A> ivisitor = {ResultsType(velonly), ExecEnv(true, direct, cpu_x86)};
#endif
  RHSVisitor<S> rvisitor;

  //
  // update rhs first
  //

  std::cout << "  Solving for BEM RHS" << std::endl;

  // loop over boundary collections
  for (auto &targ : _bdry) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl;

    // transform the collection according to prescribed motion
    std::visit([=](auto& elem) { elem.transform(_time); }, targ);

    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);

    // accumulate from vorticity
    for (auto &src : _vort) {
      std::visit(ivisitor, src, targ);
    }

    // divide by factor and add freestream
    std::visit([=](auto& elem) { elem.finalize_vels(_fs); }, targ);

    // include the effects of the motion of the parent body
    std::visit([=](auto& elem) { elem.add_body_motion(-1.0, _time); }, targ);

    // find portion of RHS vector
    const size_t tstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, targ);
    const size_t tnum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, targ);

    // have to convert these velocities into BCs based on the target element and BC type!
    // and send it over to the BEM system
    std::vector<S> rhs = std::visit(rvisitor, targ);

    // optionally augment with an additional value
    if (std::visit([=](auto& elem) { return elem.is_augmented(); }, targ)) {
      assert(tnum-rhs.size()==1 && "Number of augmented rows is not 1");
      assert(std::holds_alternative<Surfaces<S>>(targ) && "Augmented boundary is not Surface!");

      // make this easier
      Surfaces<S>& surf = std::get<Surfaces<S>>(targ);

      // first, account for the circulation change due to the change in rotational speed
      const S self_circ = surf.get_last_body_circ();

      // and include the error from the previous solve
      const S last_error = surf.get_last_body_circ_error();

      // finally, find the amount of circulation re-absorbed by the body
      const S absorbed_circ = surf.get_reabsorbed();

      // combine and append
      std::cout << "    components of rhs: " << self_circ << " " << absorbed_circ << " " << last_error << std::endl;
      const S tot_circ = self_circ + absorbed_circ + last_error;
      rhs.push_back(tot_circ);
      std::cout << "    augmenting rhs with tot_circ= " << tot_circ << std::endl;
    }

    // finally, send it to the BEM
    _bem.set_rhs(tstart, tnum, rhs);
  }

  //
  // rhs is done, update A matrix now
  //

  // options from here are: rebuild every block, rebuild some blocks, rebuild no blocks
  bool rebuild_every_block = false;
  bool rebuild_some_blocks = false;

  if (_bem.is_A_current()) {
    // we've already built the A matrix once, so unmoving blocks are correct
    if (last_time == _time) {
      // even moving blocks didn't move again, so we don't have to recalculate them
    } else {
      // time is different, check every block pair for relative movement
      rebuild_some_blocks = true;
    }
  } else {
    // if this is the first call after a reset, we need to rebuild everything
    rebuild_every_block = true;
    std::cout << "  Solving for BEM matrix" << std::endl;
  }

  // actually make or remake the A matrix
  if (rebuild_every_block or rebuild_some_blocks) {

    auto start = std::chrono::steady_clock::now();

    // need this to inform bem that we need to re-init the solver
    _bem.panels_changed();

    // this is the dispatcher for Points/Surfaces on Points/Surfaces
    CoefficientVisitor<S> cvisitor;

    // loop over boundary collections
    for (auto &targ : _bdry) {
      //std::cout << "  Solving for influence coefficients on" << to_string(targ) << std::endl;

      // find portion of influence matrix
      const size_t tstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, targ);
      const size_t tnum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, targ);

      // assemble from all boundaries
      for (auto &src : _bdry) {

        // should we build/rebuild this block of the A matrix?
        bool rebuild_this_block = rebuild_every_block;
        if (rebuild_some_blocks) {
          // test for relative motion between these two blocks (translation only for now)
          std::shared_ptr<Body> tb = std::visit([=](auto& elem) { return elem.get_body_ptr(); }, targ);
          std::shared_ptr<Body> sb = std::visit([=](auto& elem) { return elem.get_body_ptr(); }, src);
          if (tb and sb) rebuild_this_block = tb->relative_motion_vs(sb, last_time, _time);
        }

        if (rebuild_this_block) {
          // find portion of influence matrix
          const size_t sstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, src);
          const size_t snum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, src);
          std::cout << "  Computing A matrix block [" << tstart << ":" << (tstart+tnum) << "] x [" << sstart << ":" << (sstart+snum) << "]" << std::endl;

          // for augmentation, find the induced velocity from the source on the target
          if (std::visit([=](auto& elem) { return elem.is_augmented(); }, src)) {
            std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
            std::visit([=](auto& elem) { elem.zero_strengths(); }, src);
            std::visit([=](auto& elem) { elem.add_unit_rot_strengths(); }, src);
            std::visit(ivisitor, src, targ);
            std::visit([=](auto& elem) { elem.zero_strengths(); }, src);
            std::visit([=](auto& elem) { elem.finalize_vels(std::array<double,Dimensions>({0.0,0.0})); }, targ);
          }

          // solve for the coefficients in this block
          Vector<S> coeffs = std::visit(cvisitor, src, targ);
          assert(coeffs.size() == tnum*snum && "Number of coefficients does not match predicted");
          // targets are rows, sources are cols
          _bem.set_block(tstart, tnum, sstart, snum, coeffs);
        }
      }
    }

    _bem.just_made_A();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    printf("    make A matrix:\t[%.4f] cpu seconds\n", (float)elapsed_seconds.count());
  }

  //
  // solve here
  //
  std::cout << "  Solving BEM for strengths" << std::endl;
  _bem.solve();
  //
  //
  //

  // copy strengths down to Points/Surfaces
  for (auto &targ : _bdry) {
    // find portion of RHS vector
    const size_t tstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, targ);
    const size_t tnum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, targ);

    // get that chunk
    Vector<S> new_s = _bem.get_str(tstart, tnum);

    // debug print
    if (false) {
      std::cout << "  Solution vector contains" << std::endl;
      for (size_t i=tstart; i<tstart+tnum; ++i) {
        std::cout << "    " << i << " \t" << new_s[i-tstart] << std::endl;
      }
    }

    // send it to the elements, including the augmented entry
    std::visit([=](auto& elem) { elem.set_str(tstart, new_s.size(), new_s);  }, targ);
  }

  // save the simulation time to compare to the next call
  last_time = _time;
}

