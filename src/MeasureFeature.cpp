/*
 * MeasureFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "MeasureFeature.h"

#include <cmath>
#include <iostream>
#include <sstream>
//#include <random>

// write out any object of parent type MeasureFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff) {
  ff.debug(os);
  return os;
}

//
// drop a single particle
//
std::vector<float>
SinglePoint::init_particles(float _ips) const {
  return std::vector<float>({m_x, m_y});
}

std::vector<float>
SinglePoint::step_particles(float _ips) const {
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


//
// drop a single particle from the emitter
//
std::vector<float>
TracerEmitter::init_particles(float _ips) const {
  return std::vector<float>();
}

std::vector<float>
TracerEmitter::step_particles(float _ips) const {
  return std::vector<float>({m_x, m_y});
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

