/*
 * FlowFeature.h - GUI-side descriptions of flow features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "Body.h"
#include "Feature.h"
#include "ElementPacket.h"
#include "json/json.hpp"

#include <memory>
#include <iostream>
#include <vector>
#include <string>

//
// Abstract class for any flow feature present initially
//
class FlowFeature : public Feature {
public:
  explicit
  FlowFeature(float _x,
              float _y,
              std::shared_ptr<Body> _bp)
    : Feature(_x, _y, true, _bp)
    {}
  virtual ~FlowFeature() = default; 
  virtual FlowFeature* copy() const = 0;

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(float) const = 0;
  virtual ElementPacket<float> step_elements(float) const = 0;
  virtual void generate_draw_geom() = 0;
#ifdef USE_IMGUI
  virtual bool draw_info_gui(const std::string, const float) = 0;
#endif 
  
#ifdef USE_IMGUI
  static bool draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &, const float);
  static void draw_feature_list(std::vector<std::unique_ptr<FlowFeature>> &, std::unique_ptr<FlowFeature> &, int &,
                                int &, bool &, int &);
#endif

  // emit particles as vector of float4
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
                 float _str = 1.0,
                 std::shared_ptr<Body> _bp = nullptr)
    : FlowFeature(_x, _y, _bp),
      m_str(_str)
    {}
  SingleParticle* copy() const override 
                  { return new SingleParticle(*this); }
  
  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
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
             float _soft = 0.1,
             std::shared_ptr<Body> _bp = nullptr)
    : SingleParticle(_x, _y, _str, _bp),
      m_rad(_rad),
      m_softness(_soft)
    {}
  VortexBlob* copy() const override 
              { return new VortexBlob(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
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
                 float _theta = 0.0,
                 std::shared_ptr<Body> _bp = nullptr)
    : VortexBlob(_x, _y, _str, _majrad, _soft, _bp),
      m_minrad(_minrad),
      m_theta(_theta)
    {}
  AsymmetricBlob* copy() const override 
                  { return new AsymmetricBlob(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
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
               float _stddev = 0.5,
               std::shared_ptr<Body> _bp = nullptr)
    : SingleParticle(_x, _y, _str, _bp),
      m_stddev(_stddev)
    {}
  GaussianBlob* copy() const override 
                { return new GaussianBlob(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
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
               float _str = 1.0,
               std::shared_ptr<Body> _bp = nullptr)
    : SingleParticle(_x, _y, _str, _bp),
      m_xsize(_xsize),
      m_ysize(_ysize)
    {}
  UniformBlock* copy() const override 
                { return new UniformBlock(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
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
                float _minstr = -0.5,
                float _maxstr = 0.5,
                int   _num = 100,
                std::shared_ptr<Body> _bp = nullptr)
    : FlowFeature(_x, _y, _bp),
      m_xsize(_xsize),
      m_ysize(_ysize),
      m_minstr(_minstr),
      m_maxstr(_maxstr),
      m_num(_num)
    {}
  BlockOfRandom* copy() const override 
                 { return new BlockOfRandom(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
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
                  float _str = 0.1,
                  std::shared_ptr<Body> _bp = nullptr)
    : SingleParticle(_x, _y, _str, _bp)
    {}
  ParticleEmitter* copy() const override 
                   { return new ParticleEmitter(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float) override;
#endif

protected:
};

// uniformly-spaced particles

// particles from file

//
// Parser for converting json object to new feature
//
void parse_flow_json(std::vector<std::unique_ptr<FlowFeature>>&, const nlohmann::json);

