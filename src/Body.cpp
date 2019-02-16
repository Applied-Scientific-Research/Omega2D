/*
 * Body.cpp - class for an independent solid boundary
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Body.h"

#include <cassert>

//
// A single rigid body
//
// Note that Vec = std::array<double,Dimensions>
//

// delegating ctor
Body::Body() :
  Body(0.0, 0.0)
  {}

// primary constructor
Body::Body(const double _x, const double _y) :
  pos(Vec({{_x, _y}})),
  vel(Vec({{0.0, 0.0}})),
  apos(0.0),
  avel(0.0)
  {}

// getters/setters

void Body::set_name(const std::string _name) { name = _name; }
void Body::set_parent_name(const std::string _name) { parent = _name; }
std::string Body::get_name() { return name; }

void Body::set_pos(const size_t _i, const double _val) {
  assert(_i>=0 and _i<Dimensions);
  pos[_i] = _val;
}
void Body::set_pos(const size_t _i, const std::string _val) {
  assert(_i>=0 and _i<Dimensions);
  pos_func[_i] = _val;
}

Vec Body::get_pos(const double _time) {
  // for testing, return a sinusoid
  pos[0] = 0.5 * sin(2.0*_time);
  pos[1] = 0.5 * (1.0-cos(2.0*_time));
  return pos;
}

Vec Body::get_vel(const double _time) {
  // for testing, return a sinusoid
  vel[0] = cos(2.0*_time);
  vel[1] = 0.5 * 2.0*sin(2.0*_time);
  return vel;
}

double Body::get_orient(const double _time) {
  return apos;
}

double Body::get_rotvel(const double _time) {
  return avel;
}

