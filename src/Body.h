/*
 * Body.h - class for an independent solid boundary
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#define _USE_MATH_DEFINES // Required by MSVC to define M_PI,etc. in <cmath>
#include <cmath>
#include <vector>
#include <array>

//------------------------------------------------------------------------
//
// A single rigid body
//
// Draws and moves itself, so needs to know its panels? no.
//
// movement subtypes DynamicBody, KinematicBody, StaticBody
//
// how would we support inside-out bodies? say we're simulating flow inside a circle or box?
//

template <class S>
class Body {
public:
  Body();
  Body(const S, const S);
  virtual ~Body() {}

  // following are handled by the base class
  std::array<S,2> getPos();

  // following must be implemented by derived classes
  virtual std::vector<S> discretize(const S) = 0;

protected:
  // 2D position and velocity
  std::array<S,2> pos;
  std::array<S,2> vel;

  // angular position and velocity
  S apos;
  S avel;

private:
};

// delegating ctor
template <class S>
Body<S>::Body() :
  Body(0.0, 0.0)
  {}

// primary constructor
template <class S>
Body<S>::Body(const S _x, const S _y) :
  pos(std::array<S,2>({{_x, _y}})),
  vel(std::array<S,2>({{0.0, 0.0}})),
  apos(0.0),
  avel(0.0)
  {}

// getters/setters
template <class S>
std::array<S,2> Body<S>::getPos() {
  return pos;
}

