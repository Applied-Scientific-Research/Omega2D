/*
 * BoundaryFeature.cpp - GUI-side descriptions of boundary features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "BoundaryFeature2.h"

#include <cmath>
#include <iostream>
#include <sstream>

// write out any object of parent type BoundaryFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, BoundaryFeature2 const& ff) { 
  ff.debug(os);
  return os;
}

//
// Create a circle (fluid is outside circle)
//
ElementPacket<float>
SolidCircle2::init_elements(const float _ips) const {

  // how many panels?
  const size_t num_panels = std::min(10000, std::max(5, (int)(m_diam * M_PI / _ips)));

  std::cout << "Creating circle with " << num_panels << " panels" << std::endl;

  // created once
  std::vector<float>   x(num_panels*2);
  std::vector<Int>   idx(num_panels*2);
  std::vector<float> val(num_panels);

  // outside is to the left walking from one point to the next
  // so go CW around the circle starting at theta=0 (+x axis)
  for (size_t i=0; i<num_panels; i++) {
    x[2*i]     = m_x + 0.5*m_diam * std::cos(2.0 * M_PI * (float)i / (float)num_panels);
    x[2*i+1]   = m_y - 0.5*m_diam * std::sin(2.0 * M_PI * (float)i / (float)num_panels);
    idx[2*i]   = i;
    idx[2*i+1] = i+1;
    val[i]     = 0.0;
  }

  // correct the final index
  idx[2*num_panels-1] = 0;

  return ElementPacket<float>({x, idx, val});
}

void
SolidCircle2::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidCircle2::to_string() const {
  std::stringstream ss;
  ss << "solid circle at " << m_x << " " << m_y << " with diameter " << m_diam;
  return ss.str();
}

