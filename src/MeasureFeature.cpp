/*
 * MeasureFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "MeasureFeature.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <random>

// write out any object of parent type MeasureFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff) {
  ff.debug(os);
  return os;
}


//
// parse the json and dispatch the constructors
//
void parse_measure_json(std::vector<std::unique_ptr<MeasureFeature>>& _flist,
                        const nlohmann::json _jin) {

  // must have one and only one type
  if (_jin.count("type") != 1) return;

  const std::string ftype = _jin["type"];

  if      (ftype == "tracer") {           _flist.emplace_back(std::make_unique<SinglePoint>()); }
  else if (ftype == "tracer emitter") {   _flist.emplace_back(std::make_unique<TracerEmitter>()); }
  else if (ftype == "tracer blob") {      _flist.emplace_back(std::make_unique<TracerBlob>()); }
  else if (ftype == "tracer line") {      _flist.emplace_back(std::make_unique<TracerLine>()); }
  else if (ftype == "measurement line") { _flist.emplace_back(std::make_unique<MeasurementLine>()); }
  else if (ftype == "measurement grid") { _flist.emplace_back(std::make_unique<GridPoints>()); }
  else {
    std::cout << "  type " << ftype << " does not name an available measurement feature, ignoring" << std::endl;
    return;
  }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  found " << ftype << std::endl;
}


//
// Create a single measurement point
//
std::vector<float>
SinglePoint::init_particles(float _ips) const {
  // created once
  if (this->is_enabled()) return std::vector<float>({m_x, m_y});
  else return std::vector<float>();
}

std::vector<float>
SinglePoint::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
SinglePoint::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SinglePoint::to_string() const {
  std::stringstream ss;
  ss << "single field point at " << m_x << " " << m_y;
  return ss.str();
}

void
SinglePoint::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
SinglePoint::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer";
  j["center"] = {m_x, m_y};
  j["enabled"] = m_enabled;
  return j;
}


//
// Create a single, stable point which emits Lagrangian points
//
std::vector<float>
TracerEmitter::init_particles(float _ips) const {
  // is not a measurement point in itself
  // but if it was, we could use the local velocity to help generate points at any given time
  return std::vector<float>();
}

std::vector<float>
TracerEmitter::step_particles(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> zmean_dist(-0.5, 0.5);

  // emits one per step, jittered slightly
  return std::vector<float>({m_x + _ips*zmean_dist(gen),
                             m_y + _ips*zmean_dist(gen)});
}

void
TracerEmitter::debug(std::ostream& os) const {
  os << to_string();
}

std::string
TracerEmitter::to_string() const {
  std::stringstream ss;
  ss << "tracer emitter at " << m_x << " " << m_y << " spawning tracers every step";
  return ss.str();
}

void
TracerEmitter::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
TracerEmitter::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer emitter";
  j["center"] = {m_x, m_y};
  j["enabled"] = m_enabled;
  return j;
}


//
// Create a circle of tracer points
//
std::vector<float>
TracerBlob::init_particles(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> zmean_dist(-0.5, 0.5);

  // create a new vector to pass on
  std::vector<float> x;

  // what size 2D integer array will we loop over
  int irad = 1 + m_rad / _ips;
  //std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-irad; j<=irad; ++j) {

    // how far from the center are we?
    float dr = sqrt((float)(i*i+j*j)) * _ips;
    if (dr < m_rad) {
      // create a particle here
      x.emplace_back(m_x + _ips*((float)i+zmean_dist(gen)));
      x.emplace_back(m_y + _ips*((float)j+zmean_dist(gen)));
    }
  }
  }

  return x;
}

std::vector<float>
TracerBlob::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
TracerBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
TracerBlob::to_string() const {
  std::stringstream ss;
  ss << "tracer blob at " << m_x << " " << m_y << " with radius " << m_rad;
  return ss.str();
}

void
TracerBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_rad = j["rad"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
TracerBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer blob";
  j["center"] = {m_x, m_y};
  j["rad"] = m_rad;
  j["enabled"] = m_enabled;
  return j;
}


//
// Create a line of tracer points
//
std::vector<float>
TracerLine::init_particles(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

  // how many points do we need?
  float llen = std::sqrt( std::pow(m_xf-m_x, 2) + std::pow(m_yf-m_y, 2) );
  int ilen = 1 + llen / _ips;

  // loop over integer indices
  for (int i=0; i<ilen; ++i) {

    // how far along the line?
    float frac = (float)i / (float)(ilen-1);

    // create a particle here
    x.emplace_back((1.0-frac)*m_x + frac*m_xf);
    x.emplace_back((1.0-frac)*m_y + frac*m_yf);
  }

  return x;
}

std::vector<float>
TracerLine::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
TracerLine::debug(std::ostream& os) const {
  os << to_string();
}

std::string
TracerLine::to_string() const {
  std::stringstream ss;
  ss << "tracer line from " << m_x << " " << m_y << " to " << m_xf << " " << m_yf;
  return ss.str();
}

void
TracerLine::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  const std::vector<float> e = j["end"];
  m_xf = e[0];
  m_yf = e[1];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
TracerLine::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer line";
  j["center"] = {m_x, m_y};
  j["end"] = {m_xf, m_yf};
  j["enabled"] = m_enabled;
  return j;
}


//
// Create a line of static measurement points
//
std::vector<float>
MeasurementLine::init_particles(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

  // how many points do we need?
  float llen = std::sqrt( std::pow(m_xf-m_x, 2) + std::pow(m_yf-m_y, 2) );
  int ilen = 1 + llen / _ips;

  // loop over integer indices
  for (int i=0; i<ilen; ++i) {

    // how far along the line?
    float frac = (float)i / (float)(ilen-1);

    // create a particle here
    x.emplace_back((1.0-frac)*m_x + frac*m_xf);
    x.emplace_back((1.0-frac)*m_y + frac*m_yf);
  }

  return x;
}

std::vector<float>
MeasurementLine::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
MeasurementLine::debug(std::ostream& os) const {
  os << to_string();
}

std::string
MeasurementLine::to_string() const {
  std::stringstream ss;
  ss << "measurement line from " << m_x << " " << m_y << " to " << m_xf << " " << m_yf;
  return ss.str();
}

void
MeasurementLine::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  const std::vector<float> e = j["end"];
  m_xf = e[0];
  m_yf = e[1];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
MeasurementLine::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement line";
  j["center"] = {m_x, m_y};
  j["end"] = {m_xf, m_yf};
  j["enabled"] = j.value("enabled", true);
  return j;
}


//
// Create a grid of static measurement points
//
std::vector<float>
GridPoints::init_particles(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

  // ignore _ips and use m_dx to define grid density

  // loop over integer indices
  for (float xp=m_x+0.5*m_dx; xp<m_xf+0.01*m_dx; xp+=m_dx) {
    for (float yp=m_y+0.5*m_dx; yp<m_yf+0.01*m_dx; yp+=m_dx) {
      // create a field point here
      x.emplace_back(xp);
      x.emplace_back(yp);
    }
  }

  return x;
}

std::vector<float>
GridPoints::step_particles(float _ips) const {
  // does not emit
  return std::vector<float>();
}

void
GridPoints::debug(std::ostream& os) const {
  os << to_string();
}

std::string
GridPoints::to_string() const {
  std::stringstream ss;
  ss << "measurement grid from " << m_x << " " << m_y << " to " << m_xf << " " << m_yf << " widh dx " << m_dx;
  return ss.str();
}

void
GridPoints::from_json(const nlohmann::json j) {
  const std::vector<float> s = j["start"];
  m_x = s[0];
  m_y = s[1];
  const std::vector<float> e = j["end"];
  m_xf = e[0];
  m_yf = e[1];
  m_dx = j["dx"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
GridPoints::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement grid";
  j["start"] = {m_x, m_y};
  j["end"] = {m_xf, m_yf};
  j["dx"] = m_dx;
  j["enabled"] = m_enabled;
  return j;
}

