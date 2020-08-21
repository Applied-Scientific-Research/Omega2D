/*
 * Feature.h - Parent class to all features
 *
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

//#include "ElementPacket.h"

#include <iostream>
#include <vector>

//
// Abstract class for any feature
//
class Feature {
public:
  explicit
  Feature(bool _enabled)
    : m_enabled(_enabled)
    {}
  ~Feature() = default;
  //virtual Feature* copy() const = 0;

  void enable() { m_enabled = true; };
  void disable() { m_enabled = false; };
  bool is_enabled() const { return m_enabled; };
  bool* addr_enabled() { return &m_enabled; };
//  virtual std::string to_string() = 0;
 // virtual void generate_draw_geom() = 0;
#ifdef USE_IMGUI
  //virtual bool draw_info_gui(const std::string, const float);
#endif

protected:
  bool m_enabled;
};

