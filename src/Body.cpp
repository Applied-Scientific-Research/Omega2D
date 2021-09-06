/*
 * Body.cpp - class for an independent solid boundary
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
 */

#include "Body.h"

#include <cassert>
#include <iostream>

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
  this_time(0.0),
  apos_func(nullptr),
  pos(Vec({{_x, _y}})),
  vel(Vec({{0.0, 0.0}})),
  apos(0.0),
  avel(0.0)
{
  // time (t) is the only variable allowed in the equations
  func_vars.push_back({"t", &this_time});
  // ane make space for the compiled functions
  pos_func.resize(Dimensions);
}

Body::~Body() {
  // free the internal memory used by tinyexpr
  for (size_t i=0; i<pos_func.size(); ++i) te_free(pos_func[i]);
}


// create and write a json object to which to add geometries
nlohmann::json
Body::to_json() const {
  nlohmann::json j;

  if (not name.empty()) j["name"] = name;
  if (not parent.empty()) j["parent"] = parent;

  // translation has to be an array
  nlohmann::json jpos = nlohmann::json::array();
  for (size_t i=0; i<Dimensions; ++i) {
    if (pos_func[i]) {
      jpos.push_back(pos_expr[i]);
    } else {
      jpos.push_back(pos[i]);
    }
  }
  j["translation"] = jpos;

  if (apos_func) {
    j["rotation"] = apos_expr;
  } else {
    j["rotation"] = apos;
  }

  return j;
}


// getters/setters

void Body::set_name(const std::string _name) { name = _name; }
void Body::set_parent_name(const std::string _name) { parent = _name; }
std::string Body::get_name() { return name; }

void Body::set_pos(const size_t _i, const double _val) {
  assert(_i>=0 and _i<Dimensions && "Invalid index into array");
  pos[_i] = _val;
}

void Body::set_pos(const size_t _i, const std::string _val) {
  assert(_i>=0 and _i<Dimensions && "Invalid index into array");
  // store the expression locally
  pos_expr[_i] = _val;
  // compile it
  int ierr = 0;
  pos_func[_i] = te_compile(_val.c_str(), func_vars.data(), 1, &ierr);
  if (pos_func[_i]) {
    std::cout << "  read expression (" << pos_expr[_i] << ")" << std::endl;
    this_time = 0.0;
    std::cout << "  testing parsed expression, with t=0, value is " << te_eval(pos_func[_i]) << std::endl;
    this_time = 1.0;
    std::cout << "                                  t=1, value is " << te_eval(pos_func[_i]) << std::endl;
    //this_time = 2.0;
    //std::cout << "                                  t=2, value is " << te_eval(pos_func[_i]) << std::endl;
  } else {
    std::cout << "  Error parsing expression (" << _val << "), near character " << ierr << std::endl;
  }
}

void Body::set_rot(const double _val) {
  apos = _val;
}

void Body::set_rot(const std::string _val) {
  // store the expression locally
  apos_expr = _val;
  // compile it
  int ierr = 0;
  apos_func = te_compile(_val.c_str(), func_vars.data(), 1, &ierr);
  if (apos_func) {
    std::cout << "  read expression (" << apos_expr << ")" << std::endl;
    this_time = 0.0;
    std::cout << "  testing parsed expression, with t=0, value is " << te_eval(apos_func) << std::endl;
    this_time = 1.0;
    std::cout << "                                  t=1, value is " << te_eval(apos_func) << std::endl;
    //this_time = 2.0;
    //std::cout << "                                  t=2, value is " << te_eval(apos_func) << std::endl;
  } else {
    std::cout << "  Error parsing expression (" << _val << "), near character " << ierr << std::endl;
  }
}

void Body::transform(const double _time) {

  // evaluate all expressions at this time
  const double dt = 1.e-5;

  for (size_t i=0; i<Dimensions; ++i) {
    if (pos_func[i]) {
      this_time = _time;
      pos[i] = te_eval(pos_func[i]);

      this_time = _time + dt;
      const double pplus = te_eval(pos_func[i]);
      this_time = _time - dt;
      const double pminus = te_eval(pos_func[i]);
      vel[i] = (pplus - pminus) / (2.0*dt);
    }
  }

  if (apos_func) {
    this_time = _time;
    apos = te_eval(apos_func);
    apos = std::remainder(apos, 2.0*M_PI);

    this_time = _time + dt;
    const double pplus = te_eval(apos_func);
    this_time = _time - dt;
    const double pminus = te_eval(apos_func);
    // 2-point first derivative estimate
    avel = (pplus - pminus) / (2.0*dt);
  }
}

Vec Body::get_pos() {
  return pos;
}
Vec Body::get_pos(const double _time) {
  // for testing, return a loop
  //pos[0] = 0.5 * sin(2.0*_time);
  //pos[1] = 0.5 * (1.0-cos(2.0*_time));

  // for realsies, get the value or evaluate the expression
  this_time = _time;
  //std::cout << "  MOVING BODY (" << get_name() << ") at time " << _time << std::endl;
  for (size_t i=0; i<Dimensions; ++i) {
    if (pos_func[i]) {
      pos[i] = te_eval(pos_func[i]);
      //std::cout << "IDEAL POS " << (0.5 * (1.0-cos(2.0*_time))) << "  AND ACTUAL " << pos[i] << std::endl;
      //std::cout << "  MOVED BODY pos[" << i << "] to " << pos[i] << std::endl;
    }
  }

  return pos;
}

Vec Body::get_vel() {
  return vel;
}
Vec Body::get_vel(const double _time) {
  // for testing, return a loop
  //vel[0] = cos(2.0*_time);
  //vel[1] = 0.5 * 2.0*sin(2.0*_time);

  // 2-point first derivative estimate
  const double dt = 1.e-5;
  for (size_t i=0; i<Dimensions; ++i) {
    if (pos_func[i]) {
      this_time = _time + dt;
      const double pplus = te_eval(pos_func[i]);
      this_time = _time - dt;
      const double pminus = te_eval(pos_func[i]);
      vel[i] = (pplus - pminus) / (2.0*dt);
      //std::cout << "      VEL " << (0.5 * 2.0*sin(2.0*_time)) << "  AND ACTUAL " << vel[i] << std::endl;
      //std::cout << "        used " << pplus << " and " << pminus << std::endl;
    }
  }

  return vel;
}

double Body::get_orient() {
  return apos;
}
double Body::get_orient(const double _time) {
  // for realsies, get the value or evaluate the expression
  this_time = _time;
  //std::cout << "  ROTATING BODY (" << get_name() << ") at time " << _time << std::endl;
  if (apos_func) {
    apos = te_eval(apos_func);
    apos = std::remainder(apos, 2.0*M_PI);
    //std::cout << "  ROTATED BODY apos to " << apos << std::endl;
  }

  return apos;
}

double Body::get_rotvel() {
  return avel;
}
double Body::get_rotvel(const double _time) {

  // 2-point first derivative estimate
  const double dt = 1.e-5;
  if (apos_func) {
    this_time = _time + dt;
    const double pplus = te_eval(apos_func);
    this_time = _time - dt;
    const double pminus = te_eval(apos_func);
    avel = (pplus - pminus) / (2.0*dt);
  }

  return avel;
}

// compare motion vs another Body
bool Body::relative_motion_vs(std::shared_ptr<Body> _other, const double _last, const double _current) {
  bool motion = false;

  const Vec this_old_pos = get_pos(_last);
  const Vec other_old_pos = _other->get_pos(_last);
  const Vec this_new_pos = get_pos(_current);
  const Vec other_new_pos = _other->get_pos(_current);

  for (size_t i=0; i<Dimensions; ++i) {
    const double relative = this_new_pos[i] - this_old_pos[i] - other_new_pos[i] + other_old_pos[i];
    if (std::abs(relative) > 4.0*std::numeric_limits<double>::epsilon()) motion = true;
  }

  // how do we account for rotation?

  // the correct way - compute the true relative motion between the two objects
  // the hack-y way - just compare positions and orientations separately <- do this one for now

  const double this_old_theta = get_orient(_last);
  const double other_old_theta = _other->get_orient(_last);
  const double this_new_theta = get_orient(_current);
  const double other_new_theta = _other->get_orient(_current);
  const double relative = this_new_theta - this_old_theta - other_new_theta + other_old_theta;
  if (std::abs(relative) > 4.0*std::numeric_limits<double>::epsilon()) motion = true;

  //std::cout << "  relative_motion_vs " << this << " " << _other << " returns " << motion << std::endl;

  return motion;
}

