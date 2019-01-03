/*
 * MeasureFeature.h - GUI-side descriptions of flow measurement features
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <iostream>
#include <vector>

//
// Abstract class for any measurement feature (streamlines, rakes, tracers, etc.) present initially
//
class MeasureFeature {
public:
  explicit
  MeasureFeature(float _x, float _y, bool _moves)
    : m_x(_x),
      m_y(_y),
      m_is_lagrangian(_moves)
    {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual std::vector<float> init_particles(float) const = 0;
  virtual std::vector<float> step_particles(float) const = 0;

protected:
  float m_x;
  float m_y;
  bool m_is_lagrangian;
};

std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff);


//
// types of measurement features:
//
// single origin point, continuous tracer emitter
// single set of tracer particles
// fixed set of field points
// periodic rake tracer emitter
// grid of fixed field points
// solid block (square, circle) of tracers
//



//
// Concrete class for a single particle
//
class SinglePoint : public MeasureFeature {
public:
  SinglePoint(float _x, float _y, bool _moves)
    : MeasureFeature(_x, _y, _moves)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
  //float m_str;
};


//
// Concrete class for a particle emitter (one per frame)
//
class TracerEmitter : public SinglePoint {
public:
  TracerEmitter(float _x, float _y)
    : SinglePoint(_x, _y, false)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
  // eventually implement frequency, but for now, once per step
  //float m_frequency;
};

