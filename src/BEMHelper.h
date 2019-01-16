/*
 * BEMHelper.h - non-class to coordinate solving the BEM problem
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Collection.h"
#include "CollectionHelper.h"
#include "NewInfluence.h"
#include "Coeffcients.h"
#include "RHS.h"
#include "BEM.h"

#include <cstdlib>
#include <iostream>
#include <array>
#include <vector>
#include <cassert>


//
// helper function to solve BEM equations on given state
//
template <class S, class A, class I>
void solve_bem(const std::array<double,Dimensions>& _fs,
               std::vector<Collection>&             _vort,
               std::vector<Collection>&             _bdry,
               BEM<S,I>&                            _bem) {

  // no unknowns? no problem.
  if (_bdry.size() == 0) return;

  // need this for dispatching velocity influence calls, template param is accumulator type
  InfluenceVisitor<A> ivisitor;
  RHSVisitor rvisitor;

  std::cout << "  Solving for BEM RHS" << std::endl;

  // loop over boundary collections
  for (auto &targ : _bdry) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src : _vort) {
      std::visit(ivisitor, src, targ);
    }
    // add freestream and divide by 2pi
    std::visit([=](auto& elem) { elem.finalize_vels(_fs); }, targ);

    // find portion of RHS vector
    const size_t tstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, targ);
    const size_t tnum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, targ);

    // have to convert these velocities into BCs based on the target element and BC type!
    // and send it over to the BEM system
    std::vector<S> rhs = std::visit(rvisitor, targ);
    //_bem.set_rhs(tstart, tnum, std::visit(rvisitor, targ) );
    _bem.set_rhs(tstart, tnum, rhs);
  }

  // check to see if we even have to make the A matrix

  if (not _bem.is_A_current()) {
    std::cout << "  Solving for BEM matrix" << std::endl;
    auto start = std::chrono::system_clock::now();

    // this is the dispatcher for Points/Surfaces on Points/Surfaces
    CoefficientVisitor cvisitor;

    // loop over boundary collections
    for (auto &targ : _bdry) {
      std::cout << "  Solving for influence coefficients on" << to_string(targ) << std::endl;

      // find portion of influence matrix
      const size_t tstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, targ);
      const size_t tnum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, targ);

      // assemble from all boundaries
      for (auto &src : _bdry) {

        // find portion of influence matrix
        const size_t sstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, src);
        const size_t snum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, src);

        // solve for the coefficients in this block
        Vector<S> coeffs = std::visit(cvisitor, src, targ);
        assert(coeffs.size() == tnum*snum);
        _bem.set_block(tstart, tnum, sstart, snum, coeffs);
      }
    }
    _bem.just_made_A();

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    printf("    make A matrix:\t[%.4f] cpu seconds\n", (float)elapsed_seconds.count());
  }

  std::cout << "  Solving BEM for strengths" << std::endl;
  _bem.solve();

  // copy strengths down to Points/Surfaces
  for (auto &targ : _bdry) {
    // find portion of RHS vector
    const size_t tstart = std::visit([=](auto& elem) { return elem.get_first_row(); }, targ);
    const size_t tnum = std::visit([=](auto& elem) { return elem.get_num_rows(); }, targ);

    // get that chunk
    Vector<S> new_s = _bem.get_str(tstart, tnum);
    // and send it to the elements
    std::visit([=](auto& elem) { elem.set_str(tstart, tnum, new_s);  }, targ);
  }
}

