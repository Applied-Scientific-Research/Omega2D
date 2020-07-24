/*
 * FlowFeature.cpp - GUI-side descriptions of flow features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake Hillier <blakehillier@mac.com>
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <random>

// write out any object of parent type FlowFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, FlowFeature const& ff) {
  ff.debug(os);
  return os;
}


//
// parse the json and dispatch the constructors
//
void parse_flow_json(std::vector<std::unique_ptr<FlowFeature>>& _flist,
                     const nlohmann::json _jin) {

  // must have one and only one type
  if (_jin.count("type") != 1) return;

  const std::string ftype = _jin["type"];

  if      (ftype == "single particle") {  _flist.emplace_back(std::make_unique<SingleParticle>()); }
  else if (ftype == "vortex blob") {      _flist.emplace_back(std::make_unique<VortexBlob>()); }
  else if (ftype == "gaussian blob") {    _flist.emplace_back(std::make_unique<GaussianBlob>()); }
  else if (ftype == "asymmetric blob") {  _flist.emplace_back(std::make_unique<AsymmetricBlob>()); }
  else if (ftype == "uniform block") {    _flist.emplace_back(std::make_unique<UniformBlock>()); }
  else if (ftype == "block of random") {  _flist.emplace_back(std::make_unique<BlockOfRandom>()); }
  else if (ftype == "particle emitter") { _flist.emplace_back(std::make_unique<ParticleEmitter>()); }
  else {
    std::cout << "  type " << ftype << " does not name an available flow feature, ignoring" << std::endl;
    return;
  }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  std::cout << "  found " << ftype << std::endl;
}

#ifdef USE_IMGUI
bool FlowFeature::draw_creation_gui(std::vector<std::unique_ptr<FlowFeature>> &ffs, const float ips) {
  static int item = 1;
  const char* items[] = { "single particle", "round vortex blob", "Gaussian vortex blob", "asymmetric vortex blob", "block of vorticity", "random particles", "particle emitter" };
  ImGui::Combo("type", &item, items, 7);

  // show different inputs based on what is selected
  static std::unique_ptr<FlowFeature> ff = nullptr;
  switch(item) {
    case 0: {
      // creates a single particle
      ff = std::make_unique<SingleParticle>();
      //SingleParticle::draw_creation_gui(ffs);
    } break;
    case 1: {
        // creates a blob of vorticies
      ff = std::make_unique<VortexBlob>();
        //VortexBlob::draw_creation_gui(ffs, ips);
    } break;
    case 2: {
      // a gaussian blob of multiple vorticies
      ff = std::make_unique<GaussianBlob>();
      //GaussianBlob::draw_creation_gui(ffs, ips);
    } break;
    case 3: {
      // an asymmetric blob of multiple vorticies
      ff = std::make_unique<AsymmetricBlob>();
      //AsymmetricBlob::draw_creation_gui(ffs, ips);
    } break;
    case 4: {
      // particles in a rectangle
      ff = std::make_unique<UniformBlock>();
      //UniformBlock::draw_creation_gui(ffs, ips);
    } break;
    case 5: {
      // random particles in a rectangle
      ff = std::make_unique<BlockOfRandom>();
      //BlockOfRandom::draw_creation_gui(ffs);
    } break;
    case 6: {
      // create a particle emitter
      ff = std::make_unique<ParticleEmitter>();
      //ParticleEmitter::draw_creation_gui(ffs);
    } break;
  }

  bool created = false;
  if (ff->draw_info_gui("Add", ips)) {
    ffs.emplace_back(std::move(ff));
    ff = nullptr;
    created = true;
    ImGui::CloseCurrentPopup();
  }
 
  ImGui::SameLine(); 
  if (ImGui::Button("Cancel", ImVec2(120,0))) { 
    ImGui::CloseCurrentPopup();
  }

  ImGui::EndPopup();
  return created;
}
#endif

//
// important feature: convert flow feature definition into actual float4 particles
//
// each 4 floats is one particle's: x, y, strength, vdelta (radius)
//

//
// drop a single particle
//
std::vector<float>
SingleParticle::init_particles(float _ips) const {
  if (this->is_enabled()) return std::vector<float>({m_x, m_y, m_str, 0.0});
  else return std::vector<float>();
}

std::vector<float>
SingleParticle::step_particles(float _ips) const {
  return std::vector<float>();
}

void
SingleParticle::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SingleParticle::to_string() const {
  std::stringstream ss;
  ss << "single particle at " << m_x << " " << m_y << " with strength " << m_str;
  return ss.str();
}

void
SingleParticle::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_str = j["strength"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
SingleParticle::to_json() const {
  nlohmann::json j;
  j["type"] = "single particle";
  j["center"] = {m_x, m_y};
  j["strength"] = m_str;
  j["enabled"] = m_enabled;
  return j;
}

void SingleParticle::generate_draw_geom() {
  const float diam = 0.01;
  std::unique_ptr<SolidCircle> tmp = std::make_unique<SolidCircle>(nullptr, true, m_x, m_y, diam);
  m_draw = tmp->init_elements(diam/25.0);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
bool SingleParticle::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  // a single vortex particle
  static float xc[2] = {m_x, m_y};
  static float stren = m_str;
  const std::string buttonText = action+" single particle";
 
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("strength", &stren, -1.0f, 1.0f, "%.4f");
  ImGui::TextWrapped("This feature will add 1 particle");
  if (ImGui::Button(buttonText.c_str())) {
    m_x = xc[0];
    m_y = xc[1];
    m_str = stren;
    generate_draw_geom();
    add = true;
  }
  return add;
}
#endif

//
// make a circular vortex blob with soft transition
//
std::vector<float>
VortexBlob::init_particles(float _ips) const {
  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

  // what size 2D integer array will we loop over
  int irad = 1 + (m_rad + 0.5*m_softness) / _ips;
  std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // and a counter for the total circulation
  double tot_circ = 0.0;

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-irad; j<=irad; ++j) {

    // how far from the center are we?
    float dr = std::sqrt((float)(i*i+j*j)) * _ips;
    if (dr < m_rad + 0.5*m_softness) {

      // create a particle here
      x.emplace_back(m_x + _ips*(float)i);
      x.emplace_back(m_y + _ips*(float)j);

      // figure out the strength from another check
      double this_str = 1.0;
      if (dr > m_rad - 0.5*m_softness) {
        // create a weaker particle
        this_str = 0.5 - 0.5*std::sin(M_PI * (dr - m_rad) / m_softness);
      }
      x.emplace_back((float)this_str);
      tot_circ += this_str;

      // this is the radius - still zero for now
      x.emplace_back(0.0f);
    }
  }
  }

  // finally, normalize all particle strengths so that the whole blob
  //   has exactly the right strength
  std::cout << "  blob had " << tot_circ << " initial circulation" << std::endl;
  double str_scale = (double)m_str / tot_circ;
  for (size_t i=2; i<x.size(); i+=4) {
    x[i] = (float)((double)x[i] * str_scale);
  }

  return x;
}

std::vector<float>
VortexBlob::step_particles(float _ips) const {
  return std::vector<float>();
}

void
VortexBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
VortexBlob::to_string() const {
  std::stringstream ss;
  ss << "vortex blob at " << m_x << " " << m_y << ", radius " << m_rad << ", softness " << m_softness << ", and strength " << m_str;
  return ss.str();
}

void
VortexBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_rad = j["radius"];
  m_softness = j["softness"];
  m_str = j["strength"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
VortexBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "vortex blob";
  j["center"] = {m_x, m_y};
  j["radius"] = m_rad;
  j["softness"] = m_softness;
  j["strength"] = m_str;
  j["enabled"] = m_enabled;
  return j;
}

void VortexBlob::generate_draw_geom() {
  std::unique_ptr<SolidCircle> tmp = std::make_unique<SolidCircle>(nullptr, true, m_x, m_y, m_rad*2);
  m_draw = tmp->init_elements(m_rad/12.5);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
bool VortexBlob::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  static float xc[2] = {m_x, m_y};
  static float stren = m_str;
  static float rad = m_rad;
  static float soft = m_softness;
  const std::string buttonText = action+" vortex blob";

  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("strength", &stren, -5.0f, 5.0f, "%.4f");
  ImGui::SliderFloat("radius", &rad, ips, 1.0f, "%.4f");
  ImGui::SliderFloat("softness", &soft, ips, 1.0f, "%.4f");
  ImGui::TextWrapped("This feature will add about %d particles", (int)(0.785398175*std::pow((2 * rad + soft) / ips, 2)));
  if (ImGui::Button(buttonText.c_str())) {
    m_x = xc[0];
    m_y = xc[1];
    m_str = stren;
    m_rad = rad;
    m_softness = soft;
    generate_draw_geom();
    add = true;
  }
  return add;
}
#endif

//
// make an anymmetric vortex blob with soft transition
//
std::vector<float>
AsymmetricBlob::init_particles(float _ips) const {
  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

  // what size 2D integer array will we loop over
  int irad = 1 + (m_rad    + 0.5*m_softness) / _ips;
  int jrad = 1 + (m_minrad + 0.5*m_softness) / _ips;
  std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // and a counter for the total circulation
  double tot_circ = 0.0;

  const float st = std::sin(M_PI * m_theta / 180.0);
  const float ct = std::cos(M_PI * m_theta / 180.0);

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-jrad; j<=jrad; ++j) {

    const float dx = (float)i * _ips;
    const float dy = (float)j * _ips;

    // reproject to m_rad circular before calculating the distance to center
    float dr = std::sqrt(dx*dx + std::pow(dy*m_rad/m_minrad,2));
    if (dr < m_rad + 0.5*m_softness) {

      // create a particle here
      x.emplace_back(m_x + dx*ct - dy*st);
      x.emplace_back(m_y + dx*st + dy*ct);

      // figure out the strength from another check
      double this_str = 1.0;
      if (dr > m_rad - 0.5*m_softness) {
        // create a weaker particle
        this_str = 0.5 - 0.5*std::sin(M_PI * (dr - m_rad) / m_softness);
      }
      x.emplace_back((float)this_str);
      tot_circ += this_str;

      // this is the radius - still zero for now
      x.emplace_back(0.0f);
    }
  }
  }

  // finally, normalize all particle strengths so that the whole blob
  //   has exactly the right strength
  std::cout << "  asym blob had " << tot_circ << " initial circulation" << std::endl;
  double str_scale = (double)m_str / tot_circ;
  for (size_t i=2; i<x.size(); i+=4) {
    x[i] = (float)((double)x[i] * str_scale);
  }

  return x;
}

std::vector<float>
AsymmetricBlob::step_particles(float _ips) const {
  return std::vector<float>();
}

void
AsymmetricBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
AsymmetricBlob::to_string() const {
  std::stringstream ss;
  ss << "asymmetric blob at " << m_x << " " << m_y << ", radii " << m_rad << " " << m_minrad << ", softness " << m_softness << ", and strength " << m_str;
  return ss.str();
}

void
AsymmetricBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_softness = j["softness"];
  m_str = j["strength"];
  const std::vector<float> sc = j["scale"];
  m_rad = sc[0];
  m_minrad = sc[1];
  m_theta = j.value("rotation", 0.0);
  m_enabled = j.value("enabled", true);
}

nlohmann::json
AsymmetricBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "asymmetric blob";
  j["center"] = {m_x, m_y};
  j["softness"] = m_softness;
  j["strength"] = m_str;
  j["scale"] = {m_rad, m_minrad};;
  j["rotation"] = m_theta;
  j["enabled"] = m_enabled;
  return j;
}

void AsymmetricBlob::generate_draw_geom() {
  std::unique_ptr<SolidOval> tmp = std::make_unique<SolidOval>(nullptr, true, m_x, m_y, m_rad*2, m_minrad*2);
  m_draw = tmp->init_elements(m_rad/12.5);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
bool AsymmetricBlob::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  static float xc[2] = {m_x, m_y};
  static float stren = m_str;
  static float minrad = m_minrad;
  static float rotdeg = m_theta;
  static float soft = m_softness;
  static float rad = m_rad;
  const std::string buttonText = action+" asymmetric vortex blob";

  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("strength", &stren, -5.0f, 5.0f, "%.4f");
  ImGui::SliderFloat("major radius", &rad, ips, 1.0f, "%.4f");
  ImGui::SliderFloat("minor radius", &minrad, ips, 1.0f, "%.4f");
  ImGui::SliderFloat("softness", &soft, ips, 1.0f, "%.4f");
  ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 179.0f, "%.0f");
  ImGui::TextWrapped("This feature will add about %d particles", (int)(0.785398175*std::pow((2*rad+soft)/ips, 2)));
  if (ImGui::Button(buttonText.c_str())) {
    m_x = xc[0];
    m_y = xc[1];
    m_str = stren;
    m_minrad = minrad;
    m_theta = rotdeg;
    m_softness = soft;
    m_rad = rad;
    generate_draw_geom();
    add = true;
  }
  return add;
}
#endif

//
// make a Gaussian vortex blob
//
std::vector<float>
GaussianBlob::init_particles(float _ips) const {
  // create a new vector to pass on
  std::vector<float> x;

  if (not this->is_enabled()) return x;

  // what size 2D integer array will we loop over
  int irad = 1 + (3.0*m_stddev) / _ips;
  std::cout << "blob needs " << (-irad) << " to " << irad << " spaces" << std::endl;

  // and a counter for the total circulation
  double tot_circ = 0.0;

  // loop over integer indices
  for (int i=-irad; i<=irad; ++i) {
  for (int j=-irad; j<=irad; ++j) {

    // how far from the center are we?
    float dr = std::sqrt((float)(i*i+j*j)) * _ips;
    if (dr < 3.0*m_stddev) {

      // create a particle here
      x.emplace_back(m_x + _ips*(float)i);
      x.emplace_back(m_y + _ips*(float)j);

      // figure out the strength from another check
      double this_str = std::exp(-std::pow(dr/m_stddev, 2.0));
      x.emplace_back((float)this_str);
      tot_circ += this_str;

      // this is the radius - still zero for now
      x.emplace_back(0.0f);
    }
  }
  }

  // finally, normalize all particle strengths so that the whole blob
  //   has exactly the right strength
  std::cout << "  blob had " << tot_circ << " initial circulation" << std::endl;
  double str_scale = (double)m_str / tot_circ;
  for (size_t i=2; i<x.size(); i+=4) {
    x[i] = (float)((double)x[i] * str_scale);
  }

  return x;
}

std::vector<float>
GaussianBlob::step_particles(float _ips) const {
  return std::vector<float>();
}

void
GaussianBlob::debug(std::ostream& os) const {
  os << to_string();
}

std::string
GaussianBlob::to_string() const {
  std::stringstream ss;
  ss << "gaussian blob at " << m_x << " " << m_y << ", stddev " << m_stddev << ", and strength " << m_str;
  return ss.str();
}

void
GaussianBlob::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_stddev = j["stddev"];
  m_str = j["strength"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
GaussianBlob::to_json() const {
  nlohmann::json j;
  j["type"] = "gaussian blob";
  j["center"] = {m_x, m_y};
  j["stddev"] = m_stddev;
  j["strength"] = m_str;
  j["enabled"] = m_enabled;
  return j;
}

void GaussianBlob::generate_draw_geom() {
  std::unique_ptr<SolidCircle> tmp = std::make_unique<SolidCircle>(nullptr, true, m_x, m_y, m_stddev*6);
  m_draw = tmp->init_elements(m_stddev*6/25);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
bool GaussianBlob::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  static float stren = m_str;
  static float stddev = m_stddev;
  static float xc[2] = {m_x, m_y};
  const std::string buttonText = action+" gaussian blob";
  
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("strength", &stren, -5.0f, 5.0f, "%.4f");
  ImGui::SliderFloat("std dev", &stddev, ips, 1.0f, "%.4f");
  ImGui::TextWrapped("This feature will add about %d particles", (int)(0.785398175*std::pow((6*stddev) / ips, 2)));
  if (ImGui::Button(buttonText.c_str())) {
    m_str = stren;
    m_stddev = stddev;
    m_x = xc[0];
    m_y = xc[1];
    generate_draw_geom();
    add = true;
  }
  return add;
}
#endif

//
// make the block of regular, and uniform-strength particles
//
std::vector<float>
UniformBlock::init_particles(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // what size 2D integer array will we loop over
  int isize = 1 + m_xsize / _ips;
  int jsize = 1 + m_ysize / _ips;
  std::cout << "block needs " << isize << " by " << jsize << " particles" << std::endl;

  // create a new vector to pass on
  std::vector<float> x(4*isize*jsize);

  const float each_str = m_str / (float)(isize*jsize);

  // initialize the particles' locations and strengths, leave radius zero for now
  size_t iptr = 0;
  for (int i=0; i<isize; ++i) {
  for (int j=0; j<jsize; ++j) {
    x[iptr++] = m_x + m_xsize * (((float)i + 0.5)/(float)isize - 0.5);
    x[iptr++] = m_y + m_ysize * (((float)j + 0.5)/(float)jsize - 0.5);
    x[iptr++] = each_str;
    x[iptr++] = 0.0f;
  }
  }
  return x;
}

std::vector<float>
UniformBlock::step_particles(float _ips) const {
  return std::vector<float>();
}

void
UniformBlock::debug(std::ostream& os) const {
  os << to_string();
}

std::string
UniformBlock::to_string() const {
  std::stringstream ss;
  ss << "block of particles in [" << (m_x-0.5*m_xsize) << " " << (m_x+0.5*m_xsize) << "] ["
                                  << (m_y-0.5*m_ysize) << " " << (m_y+0.5*m_ysize) <<
                             "] with strength " << m_str;
  return ss.str();
}

void
UniformBlock::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  const std::vector<float> s = j["size"];
  m_xsize = s[0];
  m_ysize = s[1];
  m_str = j["strength"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
UniformBlock::to_json() const {
  nlohmann::json j;
  j["type"] = "uniform block";
  j["center"] = {m_x, m_y};
  j["size"] = {m_xsize, m_ysize};
  j["strength"] = m_str;
  j["enabled"] = m_enabled;
  return j;
}

void UniformBlock::generate_draw_geom() {
  std::unique_ptr<SolidRect> tmp = std::make_unique<SolidRect>(nullptr, true, m_x, m_y, m_xsize, m_ysize);
  m_draw = tmp->init_elements(m_xsize);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
bool UniformBlock::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  static float stren = m_str;
  static float xs[2] = {m_xsize, m_ysize};
  static float xc[2] = {m_x, m_y};
  const std::string buttonText = action+" block of vorticies";

  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("strength", &stren, -5.0f, 5.0f, "%.4f");
  ImGui::SliderFloat2("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
  ImGui::TextWrapped("This feature will add %d particles", (int)(xs[0]*xs[1]/std::pow(ips,2)));
  if (ImGui::Button(buttonText.c_str())) {
    m_str = stren;
    m_xsize = xs[0];
    m_ysize = xs[1];
    m_x = xc[0];
    m_y = xc[1];
    generate_draw_geom();
    add = true;
  }
  return add;
}
#endif

//
// make the block of randomly-placed and random-strength particles
//
std::vector<float>
BlockOfRandom::init_particles(float _ips) const {

  if (not this->is_enabled()) return std::vector<float>();

  // set up the random number generator
  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> zmean_dist(-1.0, 1.0);
  static std::uniform_real_distribution<> zo_dist(0.0, 1.0);

  std::vector<float> x(4*m_num);
  // initialize the particles' locations and strengths, leave radius zero for now
  for (size_t i=0; i<(size_t)m_num; ++i) {
    size_t idx = 4*i;
    x[idx+0] = m_x + m_xsize*zmean_dist(gen);
    x[idx+1] = m_y + m_ysize*zmean_dist(gen);
    x[idx+2] = m_minstr + (m_maxstr-m_minstr)*zo_dist(gen);
    x[idx+3] = 0.0f;
  }
  return x;
}

std::vector<float>
BlockOfRandom::step_particles(float _ips) const {
  return std::vector<float>();
}

void
BlockOfRandom::debug(std::ostream& os) const {
  os << to_string();
}

std::string
BlockOfRandom::to_string() const {
  std::stringstream ss;
  ss << "block of " << m_num << " particles in [" << (m_x-0.5*m_xsize) << " " << (m_x+0.5*m_xsize) << "] ["
                                                  << (m_y-0.5*m_ysize) << " " << (m_y+0.5*m_ysize) <<
                             "] with strengths [" << m_minstr << " " << m_maxstr << "]";
  return ss.str();
}

void
BlockOfRandom::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  const std::vector<float> s = j["size"];
  m_xsize = s[0];
  m_ysize = s[1];
  const std::vector<float> sr = j["strength range"];
  m_minstr = sr[0];
  m_maxstr = sr[1];
  m_num = j["num"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
BlockOfRandom::to_json() const {
  nlohmann::json j;
  j["type"] = "block of random";
  j["center"] = {m_x, m_y};
  j["size"] = {m_xsize, m_ysize};
  j["strength range"] = {m_minstr, m_maxstr};
  j["num"] = m_num;
  j["enabled"] = m_enabled;
  return j;
}

void BlockOfRandom::generate_draw_geom() {
  std::unique_ptr<SolidRect> tmp = std::make_unique<SolidRect>(nullptr, true, m_x, m_y, m_xsize*2, m_ysize*2);
  m_draw = tmp->init_elements(m_xsize*2);
  std::fill(m_draw.val.begin(), m_draw.val.end(), (m_minstr+m_maxstr)/2);
}

#ifdef USE_IMGUI
bool BlockOfRandom::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  static int npart = m_num;
  static float xs[2] = {m_xsize, m_ysize};
  static float strenLo = m_minstr;
  static float strenHi = m_maxstr;
  static float xc[2] = {m_x, m_y};
  const std::string buttonText = action+" random vorticies";

  ImGui::InputFloat2("center", xc);
  ImGui::SliderInt("number", &npart, 1, 10000);
  ImGui::SliderFloat2("box size", xs, 0.01f, 10.0f, "%.4f", 2.0f);
  ImGui::DragFloatRange2("strength range", &strenLo, &strenHi, 0.001f, -0.1f, 0.1f);
  ImGui::TextWrapped("This feature will add %d particles", npart);
  if (ImGui::Button(buttonText.c_str())) {
    m_num = npart;
    m_xsize = xs[0];
    m_ysize = xs[1];
    m_minstr = strenLo;
    m_maxstr = strenHi;
    m_x = xc[0];
    m_y = xc[1];
    generate_draw_geom();
    add = true;
  }
  return add;
}
#endif

//
// drop a single particle from the emitter
//
std::vector<float>
ParticleEmitter::init_particles(float _ips) const {
  return std::vector<float>();
}

std::vector<float>
ParticleEmitter::step_particles(float _ips) const {
  if (this->is_enabled()) return std::vector<float>({m_x, m_y, m_str, 0.0});
  else return std::vector<float>();
}

void
ParticleEmitter::debug(std::ostream& os) const {
  os << to_string();
}

std::string
ParticleEmitter::to_string() const {
  std::stringstream ss;
  ss << "particle emitter at " << m_x << " " << m_y << " spawning particles with strength " << m_str;
  return ss.str();
}

void
ParticleEmitter::from_json(const nlohmann::json j) {
  const std::vector<float> c = j["center"];
  m_x = c[0];
  m_y = c[1];
  m_str = j["strength"];
  m_enabled = j.value("enabled", true);
}

nlohmann::json
ParticleEmitter::to_json() const {
  nlohmann::json j;
  j["type"] = "particle emitter";
  j["center"] = {m_x, m_y};
  j["strength"] = m_str;
  j["enabled"] = m_enabled;
  return j;
}

void ParticleEmitter::generate_draw_geom() {
  const float diam = 0.01;
  std::unique_ptr<SolidCircle> tmp = std::make_unique<SolidCircle>(nullptr, true, m_x, m_y, diam);
  m_draw = tmp->init_elements(diam/25.0);
  std::fill(m_draw.val.begin(), m_draw.val.end(), m_str);
}

#ifdef USE_IMGUI
bool ParticleEmitter::draw_info_gui(const std::string action, const float ips) {
  bool add = false;
  static float eStren = m_str;
  static float xc[2] = {m_x, m_y};
  const std::string buttonText = action+" particle emitter";

  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("strength", &eStren, -0.1f, 0.1f, "%.4f");
  ImGui::TextWrapped("This feature will add 1 particle per time step");
  if (ImGui::Button(buttonText.c_str())) {
    m_str = eStren;
    m_x = xc[0];
    m_y = xc[1];
    generate_draw_geom();
    add = true;
  }
  return add;
}
#endif
