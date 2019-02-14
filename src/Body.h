/*
 * Body.h - class for an independent solid boundary
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"

#include <string>
#define _USE_MATH_DEFINES // Required by MSVC to define M_PI,etc. in <cmath>
#include <cmath>
#include <array>

using Vec = std::array<double,Dimensions>;

//------------------------------------------------------------------------
//
// A single rigid body
//
// Holds and reports only motions
//
// how would we support inside-out bodies? say we're simulating flow inside a circle or box?
//   those collections would be attached to the ground body - the default 0th body
//
// for use in 2D and 3D codes, should really use a custom Rotation object instead of a double
//

class Body {
public:
  Body();
  Body(const double, const double);
  Body(const std::string, const std::string);
  ~Body() = default;

  // following are handled by the base class
  Vec get_pos(const double);
  Vec get_vel(const double);
  double get_orient(const double);
  double get_rotvel(const double);

private:
  // string containing formulae to be parsed when needed
  std::string pos_func;
  std::string apos_func;

  // 2D position and velocity (initial, or constant)
  Vec pos;
  Vec vel;

  // angular position and velocity (initial, or constant)
  double apos;
  double avel;

  // enclosed volume (needed for total circulation of rotating body)
  double vol;
};

