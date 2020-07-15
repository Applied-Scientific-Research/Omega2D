/*
 * BoundaryFeature.h - GUI-side descriptions of boundary features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "Body.h"
#include "ElementPacket.h"
#include "Feature.h"
#include "json/json.hpp"
#include "Omega2D.h"
#include "Simulation.h"

#include <iostream>
#include <memory>
#include <vector>

//
// Abstract class for any boundary feature present initially
//
class BoundaryFeature : public Feature {
public:
  explicit
  BoundaryFeature(std::shared_ptr<Body> _bp,
                  bool _ext,
                  float _x,
                  float _y)
    : Feature(true),
      m_bp(_bp),
      m_external(_ext),
      m_x(_x),
      m_y(_y)
    {}
  virtual ~BoundaryFeature() = default;

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(const float) const = 0;
  //virtual std::vector<float> step_elements(const float) const = 0;
  std::shared_ptr<Body> get_body() { return m_bp; }
  virtual void generate_draw_geom() = 0;
  virtual ElementPacket<float> get_draw_packet() { return m_draw; }
#ifdef USE_IMGUI
  static int obj_movement_gui(int &, char* , char* , char* );
  static bool draw_creation_gui(std::vector<std::unique_ptr<BoundaryFeature>> &, Simulation&);
#endif

protected:
  std::shared_ptr<Body> m_bp;
  bool m_external;
  float m_x;
  float m_y;
  ElementPacket<float> m_draw;
};

std::ostream& operator<<(std::ostream& os, BoundaryFeature const& ff);

//
// Parser for converting json object to new feature
//
void parse_boundary_json(std::vector<std::unique_ptr<BoundaryFeature>>&,
                         std::shared_ptr<Body>,
                         const nlohmann::json);

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
// Concrete class for a circle
//
class SolidCircle : public BoundaryFeature {
public:
  SolidCircle(std::shared_ptr<Body> _bp = nullptr,
              bool _ext = true,
              float _x = 0.0,
              float _y = 0.0,
              float _diam = 1.0)
    : BoundaryFeature(_bp, _ext, _x, _y),
      m_diam(_diam)
    {}
  SolidCircle(const SolidCircle&) = default;
  ~SolidCircle() = default;

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui();
#endif
  void generate_draw_geom() override;

protected:
  float m_diam;
};


//
// Concrete class for an oval
//
class SolidOval : public SolidCircle {
public:
  SolidOval(std::shared_ptr<Body> _bp = nullptr,
            bool _ext = true,
            float _x = 0.0,
            float _y = 0.0,
            float _diam = 1.0,
            float _dmin = 0.5,
            float _theta = 0.0)
    : SolidCircle(_bp, _ext, _x, _y, _diam),
      m_dmin(_dmin),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui();
#endif
  void generate_draw_geom() override;

protected:
  float m_dmin;
  float m_theta;
};


//
// Concrete class for a square
//
class SolidSquare : public BoundaryFeature {
public:
  SolidSquare(std::shared_ptr<Body> _bp = nullptr,
              bool _ext = true,
              float _x = 0.0,
              float _y = 0.0,
              float _side = 1.0,
              float _theta = 0.0)
    : BoundaryFeature(_bp, _ext, _x, _y),
      m_side(_side),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui();
#endif
  void generate_draw_geom() override;

protected:
  float m_side;
  float m_theta;
};


//
// Concrete class for a rectangle
//
class SolidRect : public SolidSquare {
public:
  SolidRect(std::shared_ptr<Body> _bp = nullptr,
            bool _ext = true,
            float _x = 0.0,
            float _y = 0.0,
            float _sidex = 1.0,
            float _sidey = 0.5,
            float _theta = 0.0)
    : SolidSquare(_bp, _ext, _x, _y, _sidex, _theta),
      m_sidey(_sidey)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui();
#endif
  void generate_draw_geom() override;

protected:
  float m_sidey;
};


//
// Concrete class for a straight boundary segment
//
class BoundarySegment : public BoundaryFeature {
public:
  BoundarySegment(std::shared_ptr<Body> _bp = nullptr,
            bool _ext = true,
            float _x = 0.0,
            float _y = 0.0,
            float _xend = 1.0,
            float _yend = 0.0,
            float _normflow = 0.0,
            float _tangflow = 0.0)
    : BoundaryFeature(_bp, _ext, _x, _y),
      m_xe(_xend),
      m_ye(_yend),
      m_normflow(_normflow),
      m_tangflow(_tangflow)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui();
#endif
  void generate_draw_geom() override;

protected:
  float m_xe, m_ye;
  float m_normflow, m_tangflow;
};

/*
 Polygon Class
 This is a simple polygon with equa-length sides and angles.
*/
class SolidPolygon : public BoundaryFeature {
public:
  SolidPolygon(std::shared_ptr<Body> _bp = nullptr,
               bool _ext = true,
               float _x = 0.0,
               float _y = 0.0,
               int _numSides = 4,
               float _side = std::sqrt(2),
               float _radius = 1.0,
               float _theta = 0.0)
    : BoundaryFeature(_bp, _ext, _x, _y),
      m_numSides(_numSides),
      m_side(_side),
      m_radius(_radius),
      m_theta(_theta)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui();
#endif
  void generate_draw_geom() override;

protected:
  int m_numSides;
  float m_side;
  float m_radius;
  float m_theta;
};

/*
 NACA Class
 This is a simple polygon with equa-length sides and angles.
*/
class SolidAirfoil: public BoundaryFeature {
public:
  SolidAirfoil(std::shared_ptr<Body> _bp = nullptr,
               bool _ext = true,
               float _x = 0.0,
               float _y = 0.0,
               int _maxCamber = 2,
               int _maxCambLoc = 4,
               int _thickness = 15,
               float _theta = 0.0,
               float _chordLength = 1.0)
    : BoundaryFeature(_bp, _ext, _x, _y),
      m_maxCamber(_maxCamber),
      m_maxCambLoc(_maxCambLoc),
      m_thickness(_thickness),
      m_theta(_theta),
      m_chordLength(_chordLength)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui();
#endif
  void generate_draw_geom() override;
protected:
  int m_maxCamber;
  int m_maxCambLoc;
  int m_thickness;
  float m_theta;
  float m_chordLength;
};

