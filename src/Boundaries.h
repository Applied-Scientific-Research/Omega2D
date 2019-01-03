/*
 * Boundaries.h - a class for all boundary conditions - elements with unknown strengths
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Body.h"
#include "Panels.h"

#define _USE_MATH_DEFINES

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

//------------------------------------------------------------------------
//
// A collection of rigid bodies and their associated BEM panels
//

template <class S, class I>
class Boundaries {
public:
  Boundaries() {}

  bool exists() const { return (panels.get_n() > 0);}
  void add(Circle<S> const &);
  void make_panels(S const);
  Panels<S,I> & get_panels() { return panels; }
  Panels<S,I> const & get_panels() const { return panels; }
  const size_t get_npanels() const { return panels.get_n();}
  void reset();
  void reset_vels();
  void scale_and_add_freestream(std::array<double,2> const &);
  void find_strengths();

  // graphics pass-through calls
  void initGL(std::vector<float>&, float*, float*);
  void updateGL();
  void drawGL(std::vector<float>&, float*, float*);

protected:

private:
  // a body represents a unified set of boundary conditions, as does an inlet, a lifting line, a Kutta condition, etc.
  //   all of which contain unknowns which need to be solved whenever they move or the particles move
  // there will rarely be more than a few dozen of these in any simulation
  // HACK - hire Nicholas again to make this a vector of smart pointers, or variants, or whatever
  std::vector<Circle<S>> bodies;

  // all panels and BEM problem
  Panels<S,I> panels;
};


//
// add a Body to the current list
//
template <class S, class I>
void Boundaries<S,I>::add(Circle<S> const& _b) {
  // add it to the list (yes, it copies, and is non-heap)
  bodies.emplace_back(_b);
}

// eventually this will be something like
//void Boundaries<S,I>::add(std::unique_ptr<Circle<S>>&& pointer) {
//    bodies.emplace_back( std::move(pointer) );
//}


//
// generate the physical panels
//
template <class S, class I>
void Boundaries<S,I>::make_panels(S const _nomsep) {
  if (not panels.are_panels_made()) {
    for (auto & b: bodies) {
      panels.add_panels_from_closed_body( b.discretize(_nomsep) );
    }
    panels.panels_are_made();
  }
}

//
// clear out all data for a restart
//
template <class S, class I>
void Boundaries<S,I>::reset() {
  bodies.clear();
  panels.reset();
}

//
// in preparation for velocity-finding, zero the velocities on each element
//
template <class S, class I>
void Boundaries<S,I>::reset_vels() {
  panels.reset_vels();
}

//
// after velocity-finding, scale by 1/2pi and add the freestream
//
template <class S, class I>
void Boundaries<S,I>::scale_and_add_freestream(const std::array<double,2>& _fs) {
  panels.scale_and_add_freestream(_fs);
}

//
// tell Panels to solve the BEM equation
//
template <class S, class I>
void Boundaries<S,I>::find_strengths() {
  panels.solve_bem();
}


//
// OpenGL calls, pass on to children
//
template <class S, class I>
void Boundaries<S,I>::initGL(std::vector<float>& _projmat,
                             float* _poscolor, float* _negcolor) {
  //std::cout << "inside Boundaries.initGL" << std::endl;
  panels.initGL(_projmat, _poscolor, _negcolor);
  //for (auto const& b: bodies) {
  //  b->initGL(_projmat, _poscolor, _negcolor);
  //}
  // or init panels?
}

template <class S, class I>
void Boundaries<S,I>::updateGL() {
  //for (auto const& b: bodies) {
  //  b->updateGL();
  //}
  // or update panels?
  panels.updateGL();
}

template <class S, class I>
void Boundaries<S,I>::drawGL(std::vector<float>& _projmat,
                             float* _poscolor, float* _negcolor) {
  //std::cout << "inside Boundaries.drawGL" << std::endl;
  panels.drawGL(_projmat, _poscolor, _negcolor);
  //for (auto const& b: bodies) {
  //  b->drawGL(_projmat, _poscolor, _negcolor);
  //}
  // or draw from panels?
}
