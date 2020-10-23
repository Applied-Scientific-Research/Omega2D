/*
 * MeasureFeature.h - GUI-side descriptions of flow measurement features
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#pragma once

#include "Body.h"
#include "Omega2D.h"
#include "Feature.h"
#include "ElementPacket.h"
#include "json/json.hpp"

#include <iostream>
#include <vector>
#include <string>

//
// Abstract class for any measurement feature (streamlines, rakes, tracers, etc.) present initially
//
class MeasureFeature : public Feature {
public:
  explicit
  MeasureFeature(float _x,
                 float _y,
                 bool _moves,
                 bool _emits,
                 std::shared_ptr<Body> _bp)
    : Feature(_x, _y, true, _bp),
      m_is_lagrangian(_moves),
      m_emits(_emits)
    {}
  virtual ~MeasureFeature() {}
  virtual MeasureFeature* copy() const = 0;

  virtual void debug(std::ostream& os) const = 0;
  virtual std::string to_string() const = 0;
  virtual void from_json(const nlohmann::json) = 0;
  virtual nlohmann::json to_json() const = 0;
  virtual ElementPacket<float> init_elements(float) const = 0;
  virtual ElementPacket<float> step_elements(float) const = 0;
  virtual void generate_draw_geom() = 0;
#ifdef USE_IMGUI
  virtual bool draw_info_gui(const std::string, const float &, const float) = 0;
#endif
 
  bool moves() const { return m_is_lagrangian; }
  bool emits() const { return m_emits; }
  float jitter(const float, const float) const;
  bool get_is_lagrangian() { return m_is_lagrangian; }

#ifdef USE_IMGUI
  static bool draw_creation_gui(std::vector<std::unique_ptr<MeasureFeature>> &, const float, const float &);
  static void draw_feature_list(std::vector<std::unique_ptr<MeasureFeature>> &, std::unique_ptr<MeasureFeature> &,
                                int &, int &, bool &, int &);
#endif

protected:
  bool m_is_lagrangian;
  bool m_emits;
};

std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff);

//
// types of measurement features:
//
// single origin point, continuous tracer emitter
// fixed set of field points
// grid of fixed field points
// solid block (square, circle) of tracers
// single streamline (save all positions of a single point, draw as a line)
//


//
// Concrete class for a single measurement point
//
class SinglePoint : public MeasureFeature {
public:
  SinglePoint(float _x = 0.0,
              float _y = 0.0,
              bool _moves = false,
              bool _emits = false,
              std::shared_ptr<Body> _bp = nullptr)
    : MeasureFeature(_x, _y, _moves, _emits, _bp)
    {}
  SinglePoint* copy() const override { return new SinglePoint(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  //float m_str;
};

//
// Concrete class for a circle of tracer points
//
class MeasurementBlob : public SinglePoint {
public:
  MeasurementBlob(float _x = 0.0,
             float _y = 0.0,
             bool _moves = false,
             bool _emits = false,
             float _rad = 0.1,
             std::shared_ptr<Body> _bp = nullptr)
    : SinglePoint(_x, _y, _moves, _emits, _bp),
      m_rad(_rad)
    {}
  MeasurementBlob* copy() const override { return new MeasurementBlob(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  float m_rad;
};

//
// Concrete class for a line of measurement points
//
class MeasurementLine : public SinglePoint {
public:
  MeasurementLine(float _x = 0.0,
                  float _y = 0.0,
                  bool _moves = false,
                  bool _emits = false,
                  float _xf = 1.0,
                  float _yf = 0.0,
                  float _dx = 0.1,
                  std::shared_ptr<Body> _bp = nullptr)
    : SinglePoint(_x, _y, _moves, _emits, _bp),
      m_xf(_xf),
      m_yf(_yf),
      m_dx(_dx)
    {}
  MeasurementLine* copy() const override { return new MeasurementLine(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  float m_xf, m_yf;
  float m_dx;
};

//
// Concrete class for a grid of measurement points
//
class GridPoints : public MeasureFeature {
public:
  GridPoints(float _xs = -1.0,
             float _ys = -1.0,
             float _xf = 1.0,
             float _yf = 1.0,
             float _dx = 0.1)
    : MeasureFeature(_xs, _ys, false, false, nullptr),
      m_xf(_xf),
      m_yf(_yf),
      m_dx(_dx)
    {}
  GridPoints* copy() const override { return new GridPoints(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  float m_xf, m_yf;
  float m_dx;
};

//
// Concrete class for a grid of measurement points
//
class GridField: public MeasureFeature {
public:
  GridField(float _xs = -1.0,
            float _ys = -1.0,
            float _xf = 1.0,
            float _yf = 1.0,
            float _nx = 10,
            float _ny = 10)
    : MeasureFeature(_xs, _ys, false, false, nullptr),
      m_xf(_xf),
      m_yf(_yf),
      m_nx(_nx),
      m_ny(_ny)
    {}
  GridField* copy() const override { return new GridField(*this); }

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  void from_json(const nlohmann::json) override;
  nlohmann::json to_json() const override;
  ElementPacket<float> init_elements(float) const override;
  ElementPacket<float> step_elements(float) const override;
  void generate_draw_geom() override;
#ifdef USE_IMGUI
  bool draw_info_gui(const std::string, const float&, const float) override;
#endif

protected:
  float m_xf, m_yf;
  Int m_nx, m_ny;
};

//
// Parser for converting json object to new feature
//
void parse_measure_json(std::vector<std::unique_ptr<MeasureFeature>>&, const nlohmann::json);

