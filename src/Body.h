/*
 * Body.h - class for an independent solid boundary
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
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

#define TE_NAT_LOG
#include "tinyexpr/tinyexpr.h"
#include "json/json.hpp"

#include <string>
#include <memory>
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
  // destructor needs to call te_free on all te_expr pointers
  ~Body();// = default;

  // dump a json object for writing
  nlohmann::json to_json() const;

  // setters, as we may not construct the class at once
  void set_name(const std::string);
  void set_parent_name(const std::string);
  void set_pos(const size_t, const double);
  void set_pos(const size_t, const std::string);
  void set_rot(const double);
  void set_rot(const std::string);

  // return positional, orientation, or other data
  std::string get_name();
  Vec get_pos();
  Vec get_pos(const double);
  Vec get_vel();
  Vec get_vel(const double);

  // multiple ways to get the orientation
  double get_orient();
  double get_orient(const double);
  double get_rotvel();
  double get_rotvel(const double);

  // set and get the transform for a given time
  void transform(const double);

  // compare motion vs another Body
  bool relative_motion_vs(std::shared_ptr<Body>, const double, const double);

private:
  // a name to refer to this body and echo when asked
  std::string name;

  // string containing expression to be parsed when needed
  std::array<std::string,Dimensions> pos_expr;
  std::string apos_expr;
  // why not std::variant<double, std::string> for these?

  // needed by tinyexpr
  double this_time;
  std::vector<te_variable> func_vars;
  std::vector<te_expr*> pos_func;
  te_expr* apos_func;

  // 2D position and velocity (initial, or constant)
  Vec pos;
  Vec vel;

  // angular position and velocity (initial, or constant)
  double apos;
  double avel;

  // enclosed volume (needed for total circulation of rotating body)
  //double vol;

  // name of parent
  std::string parent;
  // pointer to parent
  std::shared_ptr<Body> pp;
};

