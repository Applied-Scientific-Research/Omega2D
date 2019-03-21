/*
 * BoundaryFeature.h - GUI-side descriptions of boundary features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Body.h"

#include "json/json.hpp"

#include <iostream>
#include <vector>
#include <memory>

//
// Abstract class for any boundary feature present initially
//
class BoundaryFeature {
public:
  explicit
  BoundaryFeature(std::shared_ptr<Body> _bp,
                  float _x,
                  float _y)
    : m_bp(_bp),
      m_x(_x),
      m_y(_y)
    {}
  virtual ~BoundaryFeature() {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(const float) const = 0;
  //virtual std::vector<float> step_elements(const float) const = 0;
  std::shared_ptr<Body> get_body() { return m_bp; }

protected:
  std::shared_ptr<Body> m_bp;
  float m_x;
  float m_y;
};

std::ostream& operator<<(std::ostream& os, BoundaryFeature const& ff);

/*
//
// Concrete class for a Kutta point (trailing edge point)
//
class KuttaCondition : public BoundaryFeature {
public:
  KuttaCondition(float _x, float _y, float _theta)
    : BoundaryFeature(_x, _y),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_theta;
};
*/

//
// Concrete class for a circle (fluid is outside circle)
//
class SolidCircle : public BoundaryFeature {
public:
  SolidCircle(std::shared_ptr<Body> _bp = nullptr,
              float _x = 0.0,
              float _y = 0.0,
              float _diam = 1.0)
    : BoundaryFeature(_bp, _x, _y),
      m_diam(_diam)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_diam;
};


//
// Concrete class for an oval (fluid is outside)
//
class SolidOval : public SolidCircle {
public:
  SolidOval(std::shared_ptr<Body> _bp = nullptr,
            float _x = 0.0,
            float _y = 0.0,
            float _diam = 1.0,
            float _dmin = 0.5,
            float _theta = 0.0)
    : SolidCircle(_bp, _x, _y, _diam),
      m_dmin(_dmin),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_dmin;
  float m_theta;
};


//
// Concrete class for a square (fluid is outside)
//
class SolidSquare : public BoundaryFeature {
public:
  SolidSquare(std::shared_ptr<Body> _bp = nullptr,
              float _x = 0.0,
              float _y = 0.0,
              float _side = 1.0,
              float _theta = 0.0)
    : BoundaryFeature(_bp, _x, _y),
      m_side(_side),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_side;
  float m_theta;
};


//
// Parser for converting json object to new feature
//
void parse_boundary_json(std::vector<std::unique_ptr<BoundaryFeature>>&,
                         std::shared_ptr<Body>,
                         const nlohmann::json);

