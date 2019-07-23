/*
 * Feature.h - Parent class to all features
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

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

  void enable() { m_enabled = true; };
  void disable() { m_enabled = false; };
  bool is_enabled() const { return m_enabled; };

protected:
  bool m_enabled;
};

