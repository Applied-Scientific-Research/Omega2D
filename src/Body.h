/*
 * Body.h - class for an independent solid boundary
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"

#include <string>
#include <memory>
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

  // setters, as we may not construct the class at once
  void set_name(const std::string);
  void set_parent_name(const std::string);
  std::string get_name();

  // return positional and orientation data
  Vec get_pos(const double);
  Vec get_vel(const double);
  double get_orient(const double);
  double get_rotvel(const double);

private:
  // a name to refer to this body and echo when asked
  std::string name;

  // string containing formulae to be parsed when needed
  std::string pos_func;
  std::string apos_func;
  // why not std::variant<double, std::string> for these?

  // 2D position and velocity (initial, or constant)
  Vec pos;
  Vec vel;

  // angular position and velocity (initial, or constant)
  double apos;
  double avel;

  // enclosed volume (needed for total circulation of rotating body)
  double vol;

  // name of parent
  std::string parent;
  // pointer to parent
  std::shared_ptr<Body> pp;
};

