/*
 * BoundaryFeature.h - GUI-side descriptions of boundary features
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <iostream>
#include <vector>

enum bdryType {
  circle,
  square
};

//
// Abstract class for any boundary feature present initially
//
class BoundaryFeature {
public:
  explicit
  BoundaryFeature(float _x, float _y, bdryType _t)
    : m_x(_x),
      m_y(_y),
      m_type(_t)
    {}
  virtual ~BoundaryFeature() {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual bdryType get_type() const = 0;
  virtual std::vector<float> get_params() const = 0;

protected:
  float m_x;
  float m_y;
  bdryType m_type;
};

std::ostream& operator<<(std::ostream& os, BoundaryFeature const& ff);


//
// Concrete class for a solid circle
//
class SolidCircle : public BoundaryFeature {
public:
  SolidCircle(float _x, float _y, float _diam)
    : BoundaryFeature(_x, _y, circle),
      m_diam(_diam)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  bdryType get_type() const override;
  std::vector<float> get_params() const override;

protected:
  float m_diam;
};


