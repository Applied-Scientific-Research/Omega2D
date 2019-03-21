/*
 * FlowFeature.h - GUI-side descriptions of flow features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "json/json.hpp"

#include <iostream>
#include <vector>

//
// Abstract class for any flow feature present initially
//
class FlowFeature {
public:
  explicit
  FlowFeature(float _x,
              float _y)
    : m_x(_x),
      m_y(_y)
    {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual std::vector<float> init_particles(float) const = 0;
  virtual std::vector<float> step_particles(float) const = 0;

  // emit particles as vector of float4

protected:
  float m_x;
  float m_y;
};

std::ostream& operator<<(std::ostream& os, FlowFeature const& ff);


//
// make intermediate abstract classes for flow and solid elements?
//


//
// Concrete class for a single particle
//
class SingleParticle : public FlowFeature {
public:
  SingleParticle(float _x = 0.0,
                 float _y = 0.0,
                 float _str = 1.0)
    : FlowFeature(_x, _y),
      m_str(_str)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
  float m_str;
};


//
// Concrete class for a vortex blob
//
class VortexBlob : public SingleParticle {
public:
  VortexBlob(float _x = 0.0,
             float _y = 0.0,
             float _str = 1.0,
             float _rad = 0.1,
             float _soft = 0.1)
    : SingleParticle(_x, _y, _str),
      m_rad(_rad),
      m_softness(_soft)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
  float m_rad;
  float m_softness;
};


//
// Concrete class for a rectangle of constant-strength particles
//
class UniformBlock : public SingleParticle {
public:
  UniformBlock(float _x = 0.0,
               float _y = 0.0,
               float _xsize = 1.0,
               float _ysize = 1.0,
               float _str = 1.0)
    : SingleParticle(_x, _y, _str),
      m_xsize(_xsize),
      m_ysize(_ysize)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

private:
  float m_xsize;
  float m_ysize;
};


//
// Concrete class for a rectangle of randomly-placed particles
//
class BlockOfRandom : public FlowFeature {
public:
  BlockOfRandom(float _x = 0.0,
                float _y = 0.0,
                float _xsize = 1.0,
                float _ysize = 1.0,
                float _minstr = -0.1,
                float _maxstr = 0.1,
                int   _num = 100)
    : FlowFeature(_x, _y),
      m_xsize(_xsize),
      m_ysize(_ysize),
      m_minstr(_minstr),
      m_maxstr(_maxstr),
      m_num(_num)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
  float m_xsize;
  float m_ysize;
  float m_minstr;
  float m_maxstr;
  int m_num;
};


//
// Concrete class for a particle emitter (one per frame)
//
class ParticleEmitter : public SingleParticle {
public:
  ParticleEmitter(float _x = 0.0,
                  float _y = 0.0,
                  float _str = 0.1)
    : SingleParticle(_x, _y, _str)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

protected:
};


// uniformly-spaced particles

// particles from file


