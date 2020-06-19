/*
 * FlowFeature.h - GUI-side descriptions of flow features
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Feature.h"

#include "json/json.hpp"

#include <iostream>
#include <vector>

//
// Abstract class for any flow feature present initially
//
class FlowFeature : public Feature {
public:
  explicit
  FlowFeature(float _x,
              float _y)
    : Feature(true),
      m_x(_x),
      m_y(_y)
    {}

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
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
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &);
#endif

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
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &, const float);
#endif

protected:
  float m_rad;
  float m_softness;
};


//
// Concrete class for an asymmetric vortex blob
//
class AsymmetricBlob : public VortexBlob {
public:
  AsymmetricBlob(float _x = 0.0,
                 float _y = 0.0,
                 float _str = 1.0,
                 float _majrad = 0.2,
                 float _minrad = 0.1,
                 float _soft = 0.1,
                 float _theta = 0.0)
    : VortexBlob(_x, _y, _str, _majrad, _soft),
      m_minrad(_minrad),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &, const float);
#endif

protected:
  float m_minrad;
  float m_theta;
};


//
// Concrete class for a vortex blob
//
class GaussianBlob : public SingleParticle {
public:
  GaussianBlob(float _x = 0.0,
               float _y = 0.0,
               float _str = 1.0,
               float _stddev = 0.5)
    : SingleParticle(_x, _y, _str),
      m_stddev(_stddev)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &, const float);
#endif

protected:
  float m_stddev;
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
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &, const float);
#endif

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
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &);
#endif

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
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;
#ifdef USE_IMGUI
  static void draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &);
#endif

protected:
};


// uniformly-spaced particles

// particles from file


//
// Parser for converting json object to new feature
//
void parse_flow_json(std::vector<std::unique_ptr<FlowFeature>>&, const nlohmann::json);

