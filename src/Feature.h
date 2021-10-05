/*
 * Feature.h - Parent class to all features
 *
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
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

#include "ElementPacket.h"

#include <iostream>
#include <vector>

//
// Abstract class for any feature
//
class Feature {
public:
  explicit
  Feature(float _x,
          float _y,
          bool _enabled,
          std::shared_ptr<Body> _bp)
    : m_x(_x),
      m_y(_y),
      m_enabled(_enabled),
      m_bp(_bp)
    {}
  ~Feature() = default;
  virtual Feature* copy() const = 0;

  virtual void debug(std::ostream &) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(float) const = 0;
  virtual void generate_draw_geom() = 0;

  void enable() { m_enabled = true; };
  void disable() { m_enabled = false; };
  bool is_enabled() const { return m_enabled; };
  bool* addr_enabled() { return &m_enabled; };
  ElementPacket<float> get_draw_packet() { return m_draw; }
  std::shared_ptr<Body> get_body() { return m_bp; }
  void set_body(std::shared_ptr<Body> _bp) { m_bp = _bp; }

protected:
  float m_x;
  float m_y;
  bool m_enabled;
  std::shared_ptr<Body> m_bp;
  ElementPacket<float> m_draw;
};

