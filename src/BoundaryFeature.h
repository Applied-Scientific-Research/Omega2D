/*
 * BoundaryFeature.h - GUI-side descriptions of boundary features
 *
 * (c)2017-21 Applied Scientific Research, Inc.
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

#include "Omega2D.h"
#include "Body.h"
#include "ElementPacket.h"
#include "Feature.h"
#include "Simulation.h"

#include "json/json.hpp"

#include <iostream>
#include <memory>
#include <vector>
#include <string>

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
    : Feature(_x, _y, true, _bp),
      m_external(_ext)
    {}

  virtual ~BoundaryFeature() = default;
  virtual BoundaryFeature* copy() const = 0;

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual std::string to_short_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual void create() = 0;
  virtual ElementPacket<float> init_elements(const float) const = 0;
  virtual std::vector<ElementPacket<float>> init_hybrid(const float) const = 0;
  virtual void generate_draw_geom() = 0;

#ifdef USE_IMGUI
  virtual bool draw_info_gui(const std::string) = 0;
  static int obj_movement_gui(int &, char* , char* , char* );
  static int draw_creation_gui(std::vector<std::unique_ptr<BoundaryFeature>> &, Simulation&);
  static void draw_feature_list(std::vector<std::unique_ptr<BoundaryFeature>> &,
                                std::unique_ptr<BoundaryFeature> &,
                                int &, int &, bool &, int &);
#endif

protected:
  bool m_external;
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
            bool _normunif = true,
            float _tangflow = 0.0)
    : BoundaryFeature(_bp, _ext, _x, _y),
      m_xe(_xend),
      m_ye(_yend),
      m_normflow(_normflow),
      m_normisuniform(_normunif),
      m_tangflow(_tangflow)
    {}
  ~BoundarySegment() = default;
  BoundarySegment* copy() const override { return new BoundarySegment(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "segmented boundary"; }
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override { }
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
#endif
  void generate_draw_geom() override;

protected:
  float m_xe, m_ye;
  float m_normflow;
  bool m_normisuniform;
  float m_tangflow;
};

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
  SolidCircle* copy() const override { return new SolidCircle(*this); }
  ~SolidCircle() = default;

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "circular cylinder"; }
  void from_json(nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override { }
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
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
  ~SolidOval() = default;
  SolidOval* copy() const override { return new SolidOval(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "oval cylinder"; }
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override { }
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
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
  ~SolidSquare() = default;
  SolidSquare* copy() const override { return new SolidSquare(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "square cylinder"; }
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override;
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
#endif
  void generate_draw_geom() override;

protected:
  float m_side;
  float m_theta;
  std::list<BoundarySegment> m_bsl;
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
  ~SolidRect() = default;
  SolidRect* copy() const override { return new SolidRect(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "rectangular cylinder"; }
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override;
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
#endif
  void generate_draw_geom() override;

protected:
  float m_sidey;
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
  ~SolidPolygon() = default;
  SolidPolygon* copy() const override { return new SolidPolygon(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "polygon cylinder"; }
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override;
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
#endif
  void generate_draw_geom() override;

protected:
  int m_numSides;
  float m_side;
  float m_radius;
  float m_theta;
  std::list<BoundarySegment> m_bsl;
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
  ~SolidAirfoil() = default;
  SolidAirfoil* copy() const override { return new SolidAirfoil(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "airfoil cylinder"; }
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override { }
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
#endif
  void generate_draw_geom() override;

protected:
  int m_maxCamber;
  int m_maxCambLoc;
  int m_thickness;
  float m_theta;
  float m_chordLength;
};

// This loads points from a .msh file
class FromMsh : public BoundaryFeature {
public:
  FromMsh(std::shared_ptr<Body> _bp = nullptr,
          bool _ext = true,
          float _x = 0.0,
          float _y = 0.0,
          float _inflow = 1.0,
          bool _inunif = true,
          float _tangflow = 0.0,
          std::string _infile = "input.msh")
    : BoundaryFeature(_bp, _ext, _x, _y),
      m_inmeanflow(_inflow),
      m_inisuniform(_inunif),
      m_tangflow(_tangflow),
      m_infile(_infile)
    {}
  ~FromMsh() = default;
  FromMsh* copy() const override { return new FromMsh(*this); }
 
  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::string to_short_string() const override { return "from msh file"; }
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  void create() override { }
  ElementPacket<float> init_elements(const float) const override;
  std::vector<ElementPacket<float>> init_hybrid(const float) const override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string) override;
#endif
  void generate_draw_geom() override;

protected:
  float m_inmeanflow;
  bool m_inisuniform;
  float m_tangflow;
  std::string m_infile;

  std::list<BoundarySegment> m_bsl;
};
