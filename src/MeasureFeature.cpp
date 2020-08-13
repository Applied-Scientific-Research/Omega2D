/*
 * MeasureFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2018-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#include "MeasureFeature.h"
#include "BoundaryFeature.h"
#include "imgui/imgui.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <random>
#include <string>

// write out any object of parent type MeasureFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, MeasureFeature const& ff) {
  ff.debug(os);
  return os;
}


//
// parse the json and dispatch the constructors
//
void parse_measure_json(std::vector<std::unique_ptr<MeasureFeature>>& _flist,
                        const nlohmann::json _jin) {

  // must have one and only one type
  if (_jin.count("type") != 1) return;

  const std::string ftype = _jin["type"];

  if      ((ftype == "tracer") || (ftype == "point")) {
    _flist.emplace_back(std::make_unique<SinglePoint>());
  } else if (ftype == "tracer emitter") {   _flist.emplace_back(std::make_unique<SinglePoint>(0.0, 0.0, false, true)); }
  else if (ftype == "tracer blob") {      _flist.emplace_back(std::make_unique<MeasurementBlob>()); }
  else if (ftype == "tracer line") {      _flist.emplace_back(std::make_unique<MeasurementLine>(0.0, 0.0, false, true)); }
  else if (ftype == "measurement line") { _flist.emplace_back(std::make_unique<MeasurementLine>()); }
  else if (ftype == "measurement grid") { _flist.emplace_back(std::make_unique<GridPoints>()); }
  else {
    std::cout << "  type " << ftype << " does not name an available measurement feature, ignoring" << std::endl;
    return;
  }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  found " << ftype << std::endl;
}

#ifdef USE_IMGUI
bool MeasureFeature::draw_creation_gui(std::vector<std::unique_ptr<MeasureFeature>> &mfs, const float ips, const float &tracerScale) {
  static int item = 0;
  static int oldItem = -1;
  const char* items[] = { "single point", "measurement circle", "measurement line", "measurement grid" };
  ImGui::Combo("type", &item, items, 4);

  // show different inputs based on what is selected
  static std::unique_ptr<MeasureFeature> mf = nullptr;
  if (oldItem != item) {
    switch(item) {
      case 0: {
        mf = std::make_unique<SinglePoint>();
      } break;
      case 1: {
        mf = std::make_unique<MeasurementBlob>();
      } break;
      case 2: {
        mf = std::make_unique<MeasurementLine>();
      } break;
      case 3: {
        mf = std::make_unique<GridPoints>();
      } break;
    }
    oldItem = item;
  }

  bool created = false;  
  if (mf->draw_info_gui("Add", tracerScale, ips)) {
    mf->generate_draw_geom();
    mfs.emplace_back(std::move(mf));
    mf = nullptr;
    oldItem = -1;
    created = true;
    ImGui::CloseCurrentPopup();
  }
  
  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    oldItem = -1;
    mf = nullptr;
    ImGui::CloseCurrentPopup();
  }

  ImGui::EndPopup();
  return created;
}
#endif

float MeasureFeature::jitter(const float _z, const float _ips) const {
  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> dist(-0.5, 0.5);
  // emits one per step, jittered slightly
  return _z+_ips*dist(gen);
}

//
// Create a single measurement point
//
ElementPacket<float>
SinglePoint::init_elements(float _ips) const {
  // created once
  std::vector<float> x = {m_x, m_y};
  std::vector<Int> idx;
  std::vector<float> vals;
  ElementPacket<float> packet({x, idx, vals, (size_t)1, (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
SinglePoint::step_elements(float _ips) const {
  if ((m_enabled) && (m_emits)) {
    std::vector<float> x = {jitter(m_x, _ips), jitter(m_y, _ips)};
    std::vector<Int> idx;
    std::vector<float> vals;
    ElementPacket<float> packet({x, idx, vals, (size_t)1, (uint8_t)0});
    if (packet.verify(packet.x.size(), Dimensions)) {
      return packet;
    } else {
      return ElementPacket<float>();
    }
  } else {
    return ElementPacket<float>();
  }
}

void
SinglePoint::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SinglePoint::to_string() const {
  std::stringstream ss;
  if (m_emits) {
    ss << "emiter";
  } else if (m_is_lagrangian) {
    ss << "tracer";
  } else {
    ss << "stationary";
  }
  ss << " point at " << m_x << " " << m_y;
  return ss.str();
}

void
SinglePoint::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits = j.value("emits", m_emits);
}

nlohmann::json
SinglePoint::to_json() const {
  nlohmann::json j;
  j["type"] = "point";
  j["center"] = {m_x, m_y};
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void SinglePoint::generate_draw_geom() {
  const float diam = 0.02;
  std::unique_ptr<SolidCircle> tmp = std::make_unique<SolidCircle>(nullptr, true, m_x, m_y, diam);
  m_draw = tmp->init_elements(diam/25.0);
}

#ifdef USE_IMGUI
bool SinglePoint::draw_info_gui(const std::string action, const float &tracerScale,
                                const float ips) {
  float xc[2] = {m_x, m_y};
  bool add = false;
  const std::string buttonText = action+" single point";

  ImGui::InputFloat2("position", xc);
  if (!m_emits) {
    ImGui::Checkbox("Point follows flow", &m_is_lagrangian);
  }
  if (!m_is_lagrangian) {
    ImGui::Checkbox("Point emits particles", &m_emits);
  }
  ImGui::TextWrapped("\nThis feature will add 1 point");
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  
  return add;
}
#endif

//
// Create a circle of tracer points
//
ElementPacket<float>
MeasurementBlob::init_elements(float _ips) const {

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<float> zmean_dist(-0.5, 0.5);

  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // what size 2D integer array will we loop over
  int irad = 1 + m_rad / _ips;
  //std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-irad; j<=irad; ++j) {

    // how far from the center are we?
    float dr = sqrt((float)(i*i+j*j)) * _ips;
    if (dr < m_rad) {
      // create a particle here
      x.emplace_back(m_x + _ips*((float)i+zmean_dist(gen)));
      x.emplace_back(m_y + _ips*((float)j+zmean_dist(gen)));
    }
  }
  }
  
  ElementPacket<float> packet({x, idx, vals, (size_t)(x.size()/2), (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }

}

ElementPacket<float>
MeasurementBlob::step_elements(float _ips) const {
  if ((m_enabled) && (m_emits)) {
    ElementPacket<float> packet = init_elements(_ips);
    for (size_t i=0; i<packet.x.size(); i++) {
      packet.x[i] = jitter(packet.x[i], _ips);
    }
    if (packet.verify(packet.x.size(), Dimensions)) {
      return packet;
    } else {
      return ElementPacket<float>();
    }
  } else {
    return ElementPacket<float>();
  }
}

void
MeasurementBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
MeasurementBlob::to_string() const {
  std::stringstream ss;
  if (m_emits) {
    ss << "emiter";
  } else if (m_is_lagrangian) {
    ss << "tracer";
  } else {
    ss << "stationary";
  }
  ss << " blob at " << m_x << " " << m_y << " with radius " << m_rad;
  return ss.str();
}

void
MeasurementBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_rad = j["rad"];
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits = j.value("emits", m_emits);
}

nlohmann::json
MeasurementBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "tracer blob";
  j["center"] = {m_x, m_y};
  j["rad"] = m_rad;
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void MeasurementBlob::generate_draw_geom() {
  std::unique_ptr<SolidCircle> tmp = std::make_unique<SolidCircle>(nullptr, true, m_x, m_y, m_rad*2.0);
  m_draw = tmp->init_elements(m_rad/12.5);
}

#ifdef USE_IMGUI
bool MeasurementBlob::draw_info_gui(const std::string action, const float &tracerScale, float ips) {
  float xc[2] = {m_x, m_y};
  bool add = false;
  const std::string buttonText = action+" circle of tracers";
 
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("radius", &m_rad, ips, 1.0f, "%.4f");
  if (!m_emits) {
    ImGui::Checkbox("Point follows flow", &m_is_lagrangian);
  }
  if (!m_is_lagrangian) {
    ImGui::Checkbox("Point emits particles", &m_emits);
  }
  ImGui::TextWrapped("This feature will add about %d field points",
                     (int)(0.785398175*std::pow(2*m_rad/(tracerScale*ips), 2)));
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];

  return add;
}
#endif

//
// Create a line of static measurement points
//
ElementPacket<float>
MeasurementLine::init_elements(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // how many points do we need?
  float llen = std::sqrt( std::pow(m_xf-m_x, 2) + std::pow(m_yf-m_y, 2) );
  int ilen = 1 + llen / m_dx;

  // loop over integer indices
  for (int i=0; i<ilen; ++i) {
    // how far along the line?
    float frac = (float)i / (float)(ilen-1);

    // create a particle here
    x.emplace_back((1.0-frac)*m_x + frac*m_xf);
    x.emplace_back((1.0-frac)*m_y + frac*m_yf);
  }
  
  ElementPacket<float> packet({x, idx, vals, (size_t)(2*ilen), (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }
}

ElementPacket<float>
MeasurementLine::step_elements(float _ips) const {
  if ((m_enabled) && (m_emits)) {
    ElementPacket<float> packet = init_elements(_ips);
    for (size_t i=0; i<packet.x.size(); i++) {
      packet.x[i] = jitter(packet.x[i], _ips);
    }
    if (packet.verify(packet.x.size(), Dimensions)) {
      return packet;
    } else {
      return ElementPacket<float>();
    }
  } else {
    return ElementPacket<float>();
  }
}

void
MeasurementLine::debug(std::ostream& os) const {
  os << to_string();
}

std::string
MeasurementLine::to_string() const {
  std::stringstream ss;
  if (m_emits) {
    ss << "emiter";
  } else if (m_is_lagrangian) {
    ss << "tracer";
  } else {
    ss << "stationary";
  }
  ss << " line from " << m_x << " " << m_y << " to " << m_xf << " " << m_yf << " with dx " << m_dx;
  return ss.str();
}

void
MeasurementLine::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  const std::vector<float> e = j["end"];
  m_xf = e[0];
  m_yf = e[1];
  m_dx = j.value("dx", 0.1);
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits = j.value("emits", m_emits);
}

nlohmann::json
MeasurementLine::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement line";
  j["center"] = {m_x, m_y};
  j["end"] = {m_xf, m_yf};
  j["dx"] = m_dx;
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void MeasurementLine::generate_draw_geom() {
  std::unique_ptr<BoundarySegment> tmp = std::make_unique<BoundarySegment>(nullptr, true, m_x, m_y,
                                                                           m_xf, m_yf, 0.0, 0.0);
  m_draw = tmp->init_elements(1.0);
}

#ifdef USE_IMGUI
bool MeasurementLine::draw_info_gui(const std::string action, const float &tracerScale, float ips) {
  float xc[2] = {m_x, m_y};
  float xf[2] = {m_xf, m_yf};
  bool add = false;
  const std::string buttonText = action+" line of measurement points";
 
  ImGui::InputFloat2("start", xc);
  ImGui::InputFloat2("finish", xf);
  ImGui::SliderFloat("dx", &m_dx, ips, 1.0f, "%.4f");
  if (!m_emits) {
    ImGui::Checkbox("Point follows flow", &m_is_lagrangian);
  }
  if (!m_is_lagrangian) {
    ImGui::Checkbox("Point emits particles", &m_emits);
  }
  ImGui::TextWrapped("This feature will add about %d field points",
		     1+(int)(std::sqrt(std::pow(xf[0]-xc[0],2)+std::pow(xf[1]-xc[1],2))/m_dx));
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_xf = xf[0];
  m_yf = xf[1];

  return add;
}
#endif


//
// Create a grid of static measurement points
//
ElementPacket<float>
GridPoints::init_elements(float _ips) const {

  // create a new vector to pass on
  std::vector<float> x;
  std::vector<Int> idx;
  std::vector<float> vals;

  // ignore _ips and use m_dx to define grid density
  // loop over integer indices
  for (float xp=m_x+0.5*m_dx; xp<m_xf+0.01*m_dx; xp+=m_dx) {
    for (float yp=m_y+0.5*m_dx; yp<m_yf+0.01*m_dx; yp+=m_dx) {
      // create a field point here
      x.emplace_back(xp);
      x.emplace_back(yp);
    }
  }

  ElementPacket<float> packet({x, idx, vals, (size_t)(x.size()/2), (uint8_t)0});
  if (packet.verify(packet.x.size(), Dimensions)) {
    return packet;
  } else {
    return ElementPacket<float>();
  }

}

ElementPacket<float>
GridPoints::step_elements(float _ips) const {
  // does not emit
  return ElementPacket<float>();
}

void
GridPoints::debug(std::ostream& os) const {
  os << to_string();
}

std::string
GridPoints::to_string() const {
  std::stringstream ss;
  ss << "measurement grid from " << m_x << " " << m_y << " to " << m_xf << " " << m_yf << " with dx " << m_dx;
  return ss.str();
}

void
GridPoints::from_json(const nlohmann::json j) {
  const std::vector<float> s = j["start"];
  m_x = s[0];
  m_y = s[1];
  const std::vector<float> e = j["end"];
  m_xf = e[0];
  m_yf = e[1];
  m_dx = j["dx"];
  m_enabled = j.value("enabled", true);
  m_is_lagrangian = j.value("lagrangian", m_is_lagrangian);
  m_emits= j.value("emits", m_emits);
}

nlohmann::json
GridPoints::to_json() const {
  nlohmann::json j;
  j["type"] = "measurement grid";
  j["start"] = {m_x, m_y};
  j["end"] = {m_xf, m_yf};
  j["dx"] = m_dx;
  j["enabled"] = m_enabled;
  j["lagrangian"] = m_is_lagrangian;
  j["emits"] = m_emits;
  return j;
}

void GridPoints::generate_draw_geom() {
  const float xc = (m_x+m_xf)/2;
  const float yc = (m_y+m_yf)/2;
  std::unique_ptr<SolidRect> tmp = std::make_unique<SolidRect>(nullptr, true, xc, yc,
                                                               m_xf-m_x, m_yf-m_y, 0.0);          
  m_draw = tmp->init_elements(m_xf-m_x);
}

#ifdef USE_IMGUI
bool GridPoints::draw_info_gui(const std::string action, const float &tracer_scale, const float ips) {
  float xc[2] = {m_x, m_y};
  float xf[2] = {m_xf, m_yf};
  bool add = false;
  const std::string buttonText = action+" grid of measurement points";
 
  ImGui::InputFloat2("start", xc);
  ImGui::InputFloat2("finish", xf);
  ImGui::SliderFloat("dx", &m_dx, ips, 1.0f, "%.4f");
  ImGui::TextWrapped("This feature will add about %d field points",
		     1+(int)((xf[0]-xc[0])*(xf[1]-xc[1])/(m_dx*m_dx)));
  if (ImGui::Button(buttonText.c_str())) { add = true; }
  m_x = xc[0];
  m_y = xc[1];
  m_xf = xf[0];
  m_yf = xf[1];

  return add;
}
#endif
