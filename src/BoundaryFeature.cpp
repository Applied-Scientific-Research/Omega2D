/*
 * BoundaryFeature.cpp - GUI-side descriptions of boundary features
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "BoundaryFeature.h"

#include <cmath>
#include <iostream>
#include <sstream>

// write out any object of parent type BoundaryFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, BoundaryFeature const& ff) { 
  ff.debug(os);
  return os;
}

// various debug print methods for the subclasses
void
SolidCircle::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidCircle::to_string() const {
  std::stringstream ss;
  ss << "solid circle at " << m_x << " " << m_y << " with diameter " << m_diam;
  return ss.str();
}


//
// should we handle discretization here?
//

//
// return boundary properties for use in Simulation->Boundaries
//

bdryType
SolidCircle::get_type() const {
  return m_type;
}

std::vector<float>
SolidCircle::get_params() const {
  return std::vector<float>({m_x, m_y, m_diam});
}

