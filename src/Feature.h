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
    : m_enabled(_enabled)//,
      //m_draw()
    {}

  void enable() { m_enabled = true; };
  void disable() { m_enabled = false; };
  bool is_enabled() const { return m_enabled; };
  bool* addr_enabled() { return &m_enabled; };
  //void generate_draw_geom();

protected:
  bool m_enabled;
  //ElementPacket<float> m_draw;
};

