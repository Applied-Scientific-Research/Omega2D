/*
 * BoundaryFeature.h - GUI-side descriptions of boundary features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"

#include <iostream>
#include <vector>


//
// Abstract class for any boundary feature present initially
//
class BoundaryFeature2 {
public:
  explicit
  BoundaryFeature2(float _x, float _y)
    : m_x(_x),
      m_y(_y)
    {}
  virtual ~BoundaryFeature2() {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual ElementPacket<float> init_elements(const float) const = 0;
  //virtual std::vector<float> step_elements(const float) const = 0;

protected:
  float m_x;
  float m_y;
};

std::ostream& operator<<(std::ostream& os, BoundaryFeature2 const& ff);

/*
//
// Concrete class for a Kutta point (trailing edge point)
//
class KuttaCondition : public BoundaryFeature2 {
public:
  KuttaCondition(float _x, float _y, float _theta)
    : BoundaryFeature2(_x, _y),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_theta;
};
*/

//
// Concrete class for a circle (fluid is outside circle)
//
class SolidCircle2 : public BoundaryFeature2 {
public:
  SolidCircle2(float _x, float _y, float _diam)
    : BoundaryFeature2(_x, _y),
      m_diam(_diam)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_diam;
};

/*
//
// Concrete class for an oval (fluid is outside)
//
class SolidOval : public SolidCircle {
public:
  SolidOval(float _x, float _y, float _diam, float _dmin, float _theta)
    : SolidCircle(_x, _y, _diam),
      m_dmin(_dmin),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_dmin;
  float m_theta;
};


//
// Concrete class for a square (fluid is outside)
//
class SolidSquare : public BoundaryFeature2 {
public:
  SolidSquare(float _x, float _y, float _side, float _theta)
    : BoundaryFeature2(_x, _y),
      m_side(_side),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  ElementPacket<float> init_elements(const float) const override;

protected:
  float m_side;
  float m_theta;
};
*/

