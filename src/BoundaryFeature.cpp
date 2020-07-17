/*
 * BoundaryFeature.cpp - GUI-side descriptions of boundary features
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
 */

#include "BoundaryFeature.h"
#include "GuiHelper.h"
#include "imgui/imgui.h"
#include "imgui/imgui_stdlib.h"
#include "Simulation.h"

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>

// write out any object of parent type BoundaryFeature by dispatching to appropriate "debug" method
std::ostream& operator<<(std::ostream& os, BoundaryFeature const& ff) { 
  ff.debug(os);
  return os;
}


//
// parse the json and dispatch the constructors
//
void parse_boundary_json(std::vector<std::unique_ptr<BoundaryFeature>>& _flist,
                         std::shared_ptr<Body> _bp,
                         const nlohmann::json _jin) {

  // must have one and only one type
  if (_jin.count("geometry") != 1) return;

  const std::string ftype = _jin["geometry"];

  if      (ftype == "circle") { _flist.emplace_back(std::make_unique<SolidCircle>(_bp)); }
  else if (ftype == "oval") {   _flist.emplace_back(std::make_unique<SolidOval>(_bp)); }
  else if (ftype == "square") { _flist.emplace_back(std::make_unique<SolidSquare>(_bp)); }
  else if (ftype == "rectangle") { _flist.emplace_back(std::make_unique<SolidRect>(_bp)); }
  else if (ftype == "segment") { _flist.emplace_back(std::make_unique<BoundarySegment>(_bp)); }
  else if (ftype == "polygon") { _flist.emplace_back(std::make_unique<SolidPolygon>(_bp)); }
  else if (ftype == "airfoil") { _flist.emplace_back(std::make_unique<SolidAirfoil>(_bp)); }
  else {
    std::cout << "  type " << ftype << " does not name an available boundary feature, ignoring" << std::endl;
    return;
  }

  // and pass the json object to the specific parser
  _flist.back()->from_json(_jin);

  // finally, generate the draw information
  _flist.back()->generate_draw_geom();

  std::cout << "  found " << _flist.back()->to_string() << std::endl;
}

#ifdef USE_IMGUI
int BoundaryFeature::obj_movement_gui(int &mitem, char* strx, char* stry, char* strrad) {
  // fixed to ground      - this geometry is fixed (attached to inertial)
  // attached to previous - this geometry is attached to the previous geometry
  // according to formula - this geometry is attached to a new moving body
  const char* mitems[] = { "fixed to ground", "attached to previous", "according to formula" };
  int changed = 0;
  static int tmp = -1;
  //const char* mitems[] = { "fixed", "attached to previous", "according to formula", "dynamic" };
  ImGui::Combo("movement", &mitem, mitems, 3);
  if (tmp != mitem) { 
    tmp = mitem;
    changed += 1;
  }
  // show different inputs based on what is selected
  if (mitem == 2) {
    changed += ImGui::InputText("x position", strx, 512);
    ImGui::SameLine();
    ShowHelpMarker("Use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
    changed += ImGui::InputText("y position", stry, 512);
    ImGui::SameLine();
    ShowHelpMarker("Use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
    changed += ImGui::InputText("angular position", strrad, 512);
    ImGui::SameLine();
    ShowHelpMarker("In radians, use C-style expressions, t is time\n+ - / * % ^ ( ) pi e\nabs, sin, cos, tan, exp, log, log10, sqrt, floor, pow");
  }
  
  return changed;
}

bool BoundaryFeature::draw_creation_gui(std::vector<std::unique_ptr<BoundaryFeature>>& bfs, Simulation& sim) {
  // define movement first
  static int mitem = 0;
  static char strx[512] = "0.0*t";
  static char stry[512] = "0.0*t";
  static char strrad[512] = "0.0*t";
  int changed = BoundaryFeature::obj_movement_gui(mitem, strx, stry, strrad);

  // define geometry second
  static int item = 0;
  static int numItems = 7;
  const char* items[] = { "circle", "square", "oval", "rectangle", "segment", "polygon", "NACA 4-digit" };

  ImGui::Combo("geometry type", &item, items, numItems);

  // static bp prevents a bunch of pointers from being created during the same boundary creation
  // The switch prevents constant assignment (mainly to prevent the terminal from being flooded from messages)
  static std::shared_ptr<Body> bp = nullptr;
  if (changed) {
    switch(mitem) {
      case 0:
         // this geometry is fixed (attached to inertial)
         bp = sim.get_pointer_to_body("ground");
         break;
      case 1:
         // this geometry is attached to the previous geometry (or ground)
         bp = sim.get_last_body();
         break;
      case 2:
         // this geometry is attached to a new moving body
         bp = std::make_shared<Body>();
         bp->set_pos(0, std::string(strx));
         bp->set_pos(1, std::string(stry));
         bp->set_rot(std::string(strrad));
         break;
    }
  }

  // show different inputs based on what is selected
  bool created = false;
  static std::unique_ptr<BoundaryFeature> bf = nullptr; 
  switch(item) {
    case 0: {
      bf = std::make_unique<SolidCircle>(bp);
    } break;
    case 1: {
      bf = std::make_unique<SolidSquare>(bp);
    } break;
    case 2: {
      bf = std::make_unique<SolidOval>(bp);
    } break;
    case 3: {
      bf = std::make_unique<SolidRect>(bp);
    } break;
    case 4: {
      bf = std::make_unique<BoundarySegment>(bp);
    } break;
    case 5: {
      bf = std::make_unique<SolidPolygon>(bp);
    } break;
    case 6: {
      bf = std::make_unique<SolidAirfoil>(bp);
    }
  } // end switch for geometry

  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    ImGui::CloseCurrentPopup();
    created = false;
  }
  
  if (bf->draw_info_gui()) {
    if (mitem == 2) {
      bp->set_name(bf->to_short_string());
      sim.add_body(bp);
    }
    bfs.emplace_back(std::move(bf));
    bf = nullptr;
    created = true;
  }
  
  ImGui::SameLine();
  if (ImGui::Button("Cancel", ImVec2(120,0))) {
    ImGui::CloseCurrentPopup();
    created = false;
  }
 
  return created;
}
#endif

//
// Create a circle
//
ElementPacket<float>
SolidCircle::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // how many panels?
  const size_t num_panels = std::min(40000, std::max(5, (int)(m_diam * M_PI / _ips)));

  std::cout << "Creating circle with " << num_panels << " panels" << std::endl;

  // created once
  std::vector<float>   x(num_panels*2);
  std::vector<Int>   idx(num_panels*2);
  std::vector<float> val(num_panels);

  // outside is to the left walking from one point to the next
  // so go CW around the circle starting at theta=0 (+x axis)
  // internal flow needs the opposite
  for (size_t i=0; i<num_panels; i++) {
    x[2*i]     = m_x + 0.5*m_diam * std::cos(2.0 * M_PI * (float)i / (float)num_panels);
    x[2*i+1]   = m_y - 0.5*m_diam * std::sin(2.0 * M_PI * (float)i / (float)num_panels);
    idx[2*i]   = i;
    idx[2*i+1] = i+1;
    val[i]     = 0.0;
  }

  // correct the final index
  idx[2*num_panels-1] = 0;

  // flip the orientation of the panels
  if (not m_external) {
    for (size_t i=0; i<num_panels; i++) {
      std::swap(idx[2*i], idx[2*i+1]);
    }
  }

  return ElementPacket<float>({x, idx, val});
}

void
SolidCircle::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidCircle::to_string() const {
  std::stringstream ss;
  if (m_external) {
    ss << "solid circle";
  } else {
    ss << "circular hole";
  }
  ss << " at " << m_x << " " << m_y << " with diameter " << m_diam;
  return ss.str();
}

void
SolidCircle::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  m_diam = j["scale"];
  m_external = j.value("external", true);
  m_enabled = j.value("enabled", true);
}

nlohmann::json
SolidCircle::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "circle";
  mesh["translation"] = {m_x, m_y};
  mesh["scale"] = m_diam;
  mesh["external"] = m_external;
  mesh["enabled"] = m_enabled;
  return mesh;
}

#ifdef USE_IMGUI
bool SolidCircle::draw_info_gui() {
  static bool external = m_external;
  static float xc[2] = {m_x, m_y};
  static float diam = m_diam;
  bool add = false;

  ImGui::Checkbox("Object is in flow", &external);
  ImGui::SameLine();
  ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("diameter", &diam, 0.01f, 10.0f, "%.4f", 2.0);
  ImGui::TextWrapped("This feature will add a solid circular boundary centered at the given coordinates");
  if (ImGui::Button("Add circular boundary")) {
    m_external = external;
    m_x = xc[0];
    m_y = xc[1];
    m_diam = diam;
    add = true;
    ImGui::CloseCurrentPopup();
  }

  return add;
}
#endif

void SolidCircle::generate_draw_geom() {
  m_draw = init_elements(m_diam/25.0);
}


//
// Create an oval
//
ElementPacket<float>
SolidOval::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // use the formula for equidistant nodes?
  // from https://math.stackexchange.com/questions/172766/calculating-equidistant-points-around-an-ellipse-arc
  const bool equidistant = true;

  float circum = 0.0;
  if (equidistant) {
#if defined(__APPLE__) && defined(__clang__)
    // do it the dumb way, because XCode doesn't support comp_ellint_2
    circum = m_diam * M_PI;
#else
    circum = 4.0*0.5*m_diam*std::comp_ellint_2(1.0-std::pow(m_dmin/m_diam,2));
    std::cout << "analytic circumference is " << circum << std::endl;
#endif
  } else {
    circum = m_diam * M_PI;
  }

  // how many panels?
  size_t num_panels = std::min(40000, std::max(5, (int)(circum / _ips)));

  std::cout << "Creating oval with " << num_panels << " panels" << std::endl;

  // these are for rotating the oval
  const float st = std::sin(M_PI * m_theta / 180.0);
  const float ct = std::cos(M_PI * m_theta / 180.0);

  // arrays created once
  std::vector<float>   x(num_panels*2);
  std::vector<Int>   idx(num_panels*2);
  std::vector<float> val(num_panels);

  static float phi = 0.0;

  // outside is to the left walking from one point to the next
  // so go CW around the circle starting at theta=0 (+x axis)
  for (size_t i=0; i<num_panels; i++) {
    float theta = 2.0 * M_PI * (float)i / (float)num_panels;

    // attempt to make uniform-sized panels by adjusting theta
    if (equidistant and i>0) {
      // one method: numerically solve for the correct new angle phi
      //const float m = 1.0-std::pow(m_dmin/m_diam,2);
      // theta = 0.5*m_diam*std::ellint_2(m, phi);
      //std::cout << "  i is " << i << " phi is " << phi << " and ips is " << _ips << std::endl;

      // previous point
      const float lastx = x[2*i-2];
      const float lasty = x[2*i-1];
      // first test point
      float tdx = 0.5*m_diam * std::cos(phi);
      float tdy = -0.5*m_dmin * std::sin(phi);
      float testx = m_x + tdx*ct - tdy*st;
      float testy = m_y + tdx*st + tdy*ct;
      // distance
      float dist = std::sqrt(std::pow(testx-lastx,2)+std::pow(testy-lasty,2));
      float lastdist = dist;
      float lastphi = phi;

      // easier: march forward with small steps until the panel is the correct length
      while (dist < _ips) {
        lastdist = dist;
        lastphi = phi;
        // increment theta and test again
        phi += 1.e-3;
        tdx = 0.5*m_diam * std::cos(phi);
        tdy = -0.5*m_dmin * std::sin(phi);
        testx = m_x + tdx*ct - tdy*st;
        testy = m_y + tdx*st + tdy*ct;
        dist = std::sqrt(std::pow(testx-lastx,2)+std::pow(testy-lasty,2));
        //std::cout << "    phi " << phi << " gives dist " << dist << std::endl;
      }

      // linear interpolate to find best phi
      phi = lastphi + (phi-lastphi) * (_ips-lastdist) / (dist-lastdist);
      //std::cout << "    solution is " << phi << std::endl;
      theta = phi;

      if (theta > 2.0 * M_PI) {
        // we've gone too far! break out of the for loop!
        num_panels = i-1;
        break;
      }
    }

    // now we have a usable theta, make the node
    const float dx =  0.5*m_diam * std::cos(theta);
    const float dy = -0.5*m_dmin * std::sin(theta);
    x[2*i]     = m_x + dx*ct - dy*st;
    x[2*i+1]   = m_y + dx*st + dy*ct;
    idx[2*i]   = i;
    idx[2*i+1] = i+1;
    val[i]     = 0.0;
    //std::cout << "  node " << i << " at " << x[2*i] << " " << x[2*i+1] << std::endl;
  }

  // reset phi
  phi = 0.0;

  // resize the arrays (num_panels may have changed)
  x.resize(2*num_panels);
  idx.resize(2*num_panels);
  val.resize(num_panels);

  // correct the final index
  idx[2*num_panels-1] = 0;

  // flip the orientation of the panels
  if (not m_external) {
    for (size_t i=0; i<num_panels; i++) {
      std::swap(idx[2*i], idx[2*i+1]);
    }
  }

  return ElementPacket<float>({x, idx, val});
}

void
SolidOval::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidOval::to_string() const {
  std::stringstream ss;
  if (m_external) {
    ss << "solid oval";
  } else {
    ss << "oval hole";
  }
  ss << " at " << m_x << " " << m_y << " with diameters " << m_diam << " " << m_dmin << " rotated " << m_theta << " deg";
  return ss.str();
}

void
SolidOval::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  const std::vector<float> sc = j["scale"];
  m_diam = sc[0];
  m_dmin = sc[1];
  m_theta = j.value("rotation", 0.0);
  m_external = j.value("external", true);
  m_enabled = j.value("enabled", true);
}

nlohmann::json
SolidOval::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "oval";
  mesh["translation"] = {m_x, m_y};
  mesh["scale"] = {m_diam, m_dmin};;
  mesh["rotation"] = m_theta;
  mesh["external"] = m_external;
  mesh["enabled"] = m_enabled;
  return mesh;
}

#ifdef USE_IMGUI
bool SolidOval::draw_info_gui() {
  static bool external = m_external;
  static float xc[2] = {m_x, m_y};
  static float diam = m_diam;
  static float minordiam = m_dmin;
  static float rotdeg = m_theta;
  bool add = false;

  ImGui::Checkbox("Object is in flow", &external);
  ImGui::SameLine();
  ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("major diameter", &diam, 0.01f, 10.0f, "%.4f", 2.0);
  ImGui::SliderFloat("minor diameter", &minordiam, 0.01f, 10.0f, "%.4f", 2.0);
  ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 179.0f, "%.0f");
  ImGui::TextWrapped("This feature will add a solid oval boundary centered at the given coordinates");
  if (ImGui::Button("Add oval boundary")) {
    m_external = external;
    m_x = xc[0];
    m_y = xc[1];
    m_diam = diam;
    m_dmin = minordiam;
    m_theta = rotdeg;
    add = true;
    ImGui::CloseCurrentPopup();
  }

  return add;
}
#endif

void SolidOval::generate_draw_geom() {
  m_draw = init_elements(m_diam/60.0);
}


//
// Create a square
//
ElementPacket<float>
SolidSquare::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // how many panels?
  const size_t num_panels = 4 * std::min(10000, std::max(1, (int)(m_side / _ips)));

  std::cout << "Creating square with " << num_panels << " panels" << std::endl;

  // created once
  std::vector<float>   x(num_panels*2);
  std::vector<Int>   idx(num_panels*2);
  std::vector<float> val(num_panels);

  const float st = std::sin(M_PI * m_theta / 180.0);
  const float ct = std::cos(M_PI * m_theta / 180.0);

  // outside is to the left walking from one point to the next
  // so go CW around the body
  size_t icnt = 0;
  for (size_t i=0; i<num_panels/4; i++) {
    const float px = m_side * -0.5;
    const float py = m_side * (-0.5 + (float)i / (float)(num_panels/4));
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }
  for (size_t i=0; i<num_panels/4; i++) {
    const float px = m_side * (-0.5 + (float)i / (float)(num_panels/4));
    const float py = m_side * 0.5;
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }
  for (size_t i=0; i<num_panels/4; i++) {
    const float px = m_side * 0.5;
    const float py = m_side * (0.5 - (float)i / (float)(num_panels/4));
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }
  for (size_t i=0; i<num_panels/4; i++) {
    const float px = m_side * (0.5 - (float)i / (float)(num_panels/4));
    const float py = m_side * -0.5;
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }

  // outside is to the left walking from one point to the next
  for (size_t i=0; i<num_panels; i++) {
    idx[2*i]   = i;
    idx[2*i+1] = i+1;
    val[i]     = 0.0;
  }

  // correct the final index
  idx[2*num_panels-1] = 0;

  // flip the orientation of the panels
  if (not m_external) {
    for (size_t i=0; i<num_panels; i++) {
      std::swap(idx[2*i], idx[2*i+1]);
    }
  }

  return ElementPacket<float>({x, idx, val});
}

void
SolidSquare::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidSquare::to_string() const {
  std::stringstream ss;
  if (m_external) {
    ss << "solid square";
  } else {
    ss << "square hole";
  }
  ss << " at " << m_x << " " << m_y << " with side " << m_side << " rotated " << m_theta << " deg";
  return ss.str();
}

void
SolidSquare::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  m_side = j["scale"];
  m_theta = j.value("rotation", 0.0);
  m_external = j.value("external", true);
  m_enabled = j.value("enabled",true);
}

nlohmann::json
SolidSquare::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "square";
  mesh["translation"] = {m_x, m_y};
  mesh["scale"] = m_side;
  mesh["rotation"] = m_theta;
  mesh["external"] = m_external;
  mesh["enabled"] = m_enabled;
  return mesh;
}

#ifdef USE_IMGUI
bool SolidSquare::draw_info_gui() {
  static bool external = m_external;
  static float xc[2] = {m_x, m_y};
  static float side = m_side;
  static float rotdeg = m_theta;
  bool add = false;
  
  ImGui::Checkbox("Object is in flow", &external);
  ImGui::SameLine();
  ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("side length", &side, 0.1f, 10.0f, "%.4f");
  ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 89.0f, "%.0f");
  //ImGui::SliderAngle("orientation", &rotdeg);
  ImGui::TextWrapped("This feature will add a solid square boundary centered at the given coordinates");
  if (ImGui::Button("Add square boundary")) {
    m_external = external;
    m_x = xc[0];
    m_y = xc[1];
    m_side = side;
    m_theta = rotdeg;
    add = true;
    ImGui::CloseCurrentPopup();
  }

  return add;
}
#endif

void SolidSquare::generate_draw_geom() {
  m_draw = init_elements(m_side);
  // transform according to body position at t=0?
}

//
// Create a rectangle
//
ElementPacket<float>
SolidRect::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // how many panels?
  const size_t np_x = std::min(10000, std::max(1, (int)(m_side / _ips)));
  const size_t np_y = std::min(10000, std::max(1, (int)(m_sidey / _ips)));
  const size_t num_panels = 2 * (np_x + np_y);

  std::cout << "Creating rectangle with " << num_panels << " panels" << std::endl;

  // created once
  std::vector<float>   x(num_panels*2);
  std::vector<Int>   idx(num_panels*2);
  std::vector<float> val(num_panels);

  const float st = std::sin(M_PI * m_theta / 180.0);
  const float ct = std::cos(M_PI * m_theta / 180.0);

  // outside is to the left walking from one point to the next
  // so go CW around the body
  size_t icnt = 0;
  for (size_t i=0; i<np_y; i++) {
    const float px = m_side * -0.5;
    const float py = m_sidey * (-0.5 + (float)i / (float)(np_y));
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }
  for (size_t i=0; i<np_x; i++) {
    const float px = m_side * (-0.5 + (float)i / (float)(np_x));
    const float py = m_sidey * 0.5;
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }
  for (size_t i=0; i<np_y; i++) {
    const float px = m_side * 0.5;
    const float py = m_sidey * (0.5 - (float)i / (float)(np_y));
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }
  for (size_t i=0; i<np_x; i++) {
    const float px = m_side * (0.5 - (float)i / (float)(np_x));
    const float py = m_sidey * -0.5;
    x[icnt++] = m_x + px*ct - py*st;
    x[icnt++] = m_y + px*st + py*ct;
  }

  // outside is to the left walking from one point to the next
  for (size_t i=0; i<num_panels; i++) {
    idx[2*i]   = i;
    idx[2*i+1] = i+1;
    val[i]     = 0.0;
  }

  // correct the final index
  idx[2*num_panels-1] = 0;

  // flip the orientation of the panels
  if (not m_external) {
    for (size_t i=0; i<num_panels; i++) {
      std::swap(idx[2*i], idx[2*i+1]);
    }
  }

  return ElementPacket<float>({x, idx, val});
}

void
SolidRect::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidRect::to_string() const {
  std::stringstream ss;
  if (m_external) {
    ss << "solid rect";
  } else {
    ss << "rect hole";
  }
  ss << " at " << m_x << " " << m_y << " with sides " << m_side << " " << m_sidey << " rotated " << m_theta << " deg";
  return ss.str();
}

void
SolidRect::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  const std::vector<float> sc = j["scale"];
  m_side = sc[0];
  m_sidey = sc[1];
  m_theta = j.value("rotation", 0.0);
  m_external = j.value("external", true);
  m_enabled = j.value("enabled",true);
}

nlohmann::json
SolidRect::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "rectangle";
  mesh["translation"] = {m_x, m_y};
  mesh["scale"] = {m_side, m_sidey};
  mesh["rotation"] = m_theta;
  mesh["external"] = m_external;
  mesh["enabled"] = m_enabled;
  return mesh;
}

#ifdef USE_IMGUI
bool SolidRect::draw_info_gui() {
  static bool external = m_external;
  static float xc[2] = {m_x, m_y};
  static float side = m_side;
  static float rectside = m_sidey;
  static float rotdeg = m_theta;
  bool add = false;

  ImGui::Checkbox("Object is in flow", &external);
  ImGui::SameLine();
  ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
  ImGui::InputFloat2("center", xc);
  ImGui::SliderFloat("horizontal size", &side, 0.1f, 10.0f, "%.4f");
  ImGui::SliderFloat("vertical size", &rectside, 0.1f, 10.0f, "%.4f");
  ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 89.0f, "%.0f");
  //ImGui::SliderAngle("orientation", &rotdeg);
  ImGui::TextWrapped("This feature will add a solid rectangular boundary centered at the given coordinates");
  if (ImGui::Button("Add rectangular boundary")) {
    m_external = external;
    m_x = xc[0];
    m_y = xc[1];
    m_side = side;
    m_sidey = rectside;
    m_theta = rotdeg;
    add = true;
    ImGui::CloseCurrentPopup();
  }

  return add;
}
#endif

void SolidRect::generate_draw_geom() {
  m_draw = init_elements(m_side);
}


//
// Create a segment of a solid boundary
//
ElementPacket<float>
BoundarySegment::init_elements(const float _ips) const {

  if (not this->is_enabled()) return ElementPacket<float>();

  // how many panels?
  const float seg_length = std::sqrt(std::pow(m_xe-m_x, 2) + std::pow(m_ye-m_y, 2));
  const size_t num_panels = std::min(10000, std::max(1, (int)(seg_length / _ips)));

  std::cout << "Creating segment with " << num_panels << " panels" << std::endl;
  std::cout << "  " << to_string() << std::endl;

  // created once
  std::vector<float>   x((num_panels+1)*2);
  std::vector<Int>   idx(num_panels*2);
  std::vector<float> val(num_panels);

  // outside is to the left walking from one point to the next
  // so go CW around the body
  size_t icnt = 0;
  for (size_t i=0; i<num_panels+1; i++) {
    const float s = (float)i / (float)num_panels;
    x[icnt++] = (1.0-s)*m_x + s*m_xe;
    x[icnt++] = (1.0-s)*m_y + s*m_ye;
  }

  // outside is to the left walking from one point to the next
  for (size_t i=0; i<num_panels; i++) {
    idx[2*i]   = i;
    idx[2*i+1] = i+1;
    val[i]     = m_tangflow;
  }

  // flip the orientation of the panels
  if (not m_external) {
    for (size_t i=0; i<num_panels; i++) {
      std::swap(idx[2*i], idx[2*i+1]);
    }
  }

  return ElementPacket<float>({x, idx, val});
}

void
BoundarySegment::debug(std::ostream& os) const {
  os << to_string();
}

std::string
BoundarySegment::to_string() const {
  std::stringstream ss;
  ss << "segment from " << m_x << " " << m_y << " to " << m_xe << " " << m_ye;
  if (std::abs(m_tangflow) > std::numeric_limits<float>::epsilon()) {
    ss << " with boundary vel " << m_tangflow;
  }
  return ss.str();
}

void
BoundarySegment::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["startpt"];
  m_x = tr[0];
  m_y = tr[1];
  const std::vector<float> ep = j["endpt"];
  m_xe = ep[0];
  m_ye = ep[1];
  m_normflow = j.value("normalVel", 0.0);
  m_tangflow = j.value("tangentialVel", 0.0);
  m_external = true;//j.value("external", true);
  m_enabled = j.value("enabled",true);
}

nlohmann::json
BoundarySegment::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "segment";
  mesh["startpt"] = {m_x, m_y};
  mesh["endpt"] = {m_xe, m_ye};
  mesh["normalVel"] = m_normflow;
  mesh["tangentialVel"] = m_tangflow;
  //mesh["external"] = m_external;
  mesh["enabled"] = m_enabled;
  return mesh;
}

#ifdef USE_IMGUI
bool BoundarySegment::draw_info_gui() {
  static float xc[2] = {m_x, m_y};
  static float xe[2] = {m_xe, m_ye};
  static float tangbc = m_tangflow;
  bool add = false;

  ImGui::InputFloat2("start", xc);
  ImGui::InputFloat2("end", xe);
  ImGui::SliderFloat("force tangential flow", &tangbc, -2.0f, 2.0f, "%.1f");
  ImGui::TextWrapped("This feature will add a solid boundary segment from start to end, where fluid is on the left when marching from start to end, and positive tangential flow is as if segment is moving along vector from start to end. Make sure enough segments are created to fully enclose a volume.");
  if (ImGui::Button("Add boundary segment")) {
    m_x = xc[0];
    m_y = xc[1];
    m_xe = xe[0];
    m_ye = xe[1];
    m_tangflow = tangbc;
    add = true;
    ImGui::CloseCurrentPopup();
  }

  return add;
}
#endif

void BoundarySegment::generate_draw_geom() {
  m_draw = init_elements(1.0);
}


// Create a Polygon 
ElementPacket<float>
SolidPolygon::init_elements(const float _ips) const {
  // If object has been removed, return no elements?
  if (not this->is_enabled()) return ElementPacket<float>();

  // how many panels
  const size_t panlsPerSide = std::min(10000, std::max(1, (int)(m_side / _ips)));
  const size_t num_panels = panlsPerSide * m_numSides;
 
  std::cout << "Creating " << m_numSides << "-sided polygon with " << num_panels << " panels" << std::endl;

  // created once
  std::vector<float>   x(num_panels*2);
  std::vector<Int>   idx(num_panels*2);
  std::vector<float> val(num_panels);

  // Turning degrees into radians: deg*pi/180
  float phi = 0;
  if (m_numSides % 2 == 0) { phi = 360/(m_numSides*2.0); }
  const float st = std::sin(M_PI * (m_theta - phi) / 180.0);
  const float ct = std::cos(M_PI * (m_theta - phi) / 180.0);

  // outside is to the left walking from one point to the next
  // so go CW around the body
  // m_side * i / panlsPerSide reflects distance between two adjacent panels
  // If m_numSides is even, it seems to rotate CCW by 360/m_numSides/2 degrees in place
  float vx = 0.0;
  float vy = m_radius;
  size_t icnt = 0;
  for (int j=0; j<m_numSides; j++) {
    // Find next vertex
    const float nxtVx = m_radius * std::sin(2*M_PI*(float)(j+1)/(float)m_numSides);
    const float nxtVy = m_radius * std::cos(2*M_PI*(float)(j+1)/(float)m_numSides);
    // std::cout << '(' << vx << ',' << vy << ") -> (" << nxtVx << ',' << nxtVy << ')' << std::endl;
    for (size_t i=0; i<panlsPerSide; i++) {
      const float px = vx+(nxtVx-vx)*(float)i/(float)panlsPerSide;
      const float py = vy+(nxtVy-vy)*(float)i/(float)panlsPerSide;
      x[icnt++] = m_x + px*ct - py*st;
      x[icnt++] = m_y + px*st + py*ct;
    }
    vx = nxtVx;
    vy = nxtVy;
  }

  // outside is to the left walking from one point to the next
  for (size_t i=0; i<num_panels; i++) {
    idx[2*i]   = i;
    idx[2*i+1] = i+1;
    val[i]     = 0.0;
  }

  // correct the final index
  idx[2*num_panels-1] = 0;

  // flip the orientation of the panels
  if (not m_external) {
    for (size_t i=0; i<num_panels; i++) {
      std::swap(idx[2*i], idx[2*i+1]);
    }
  }

  return ElementPacket<float>({x, idx, val});
}

void
SolidPolygon::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidPolygon::to_string() const {
  std::stringstream ss;
  if (m_external) {
    ss << "solid polygon";
  } else {
    ss << "polygon hole";
  }
  ss << " at " << m_x << " " << m_y << " with " << m_numSides << " sides length " << m_side << " rotated " << m_theta << " deg";
  return ss.str();
}

void
SolidPolygon::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  m_numSides = j["numberSides"];
  const std::vector<float>sc = j["scale"];
  m_side = sc[0];
  m_radius = sc[1];
  m_theta = j.value("rotation", 0.0);
  m_external = j.value("external", true);
  m_enabled = j.value("enabled",true);
}

nlohmann::json
SolidPolygon::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "polygon";
  mesh["translation"] = {m_x, m_y};
  mesh["numberSides"] = m_numSides;
  mesh["scale"] = {m_side, m_radius};
  mesh["rotation"] = m_theta;
  mesh["external"] = m_external;
  mesh["enabled"] = m_enabled;
  return mesh;
}

#ifdef USE_IMGUI
bool SolidPolygon::draw_info_gui() {
  static bool external = m_external;
  static float xc[2] = {m_x, m_y};
  static float rotdeg = m_theta;
  static int numSides = m_numSides;
  static float side = m_side;
  static float rad = m_radius;
  bool add = false;
  
  ImGui::Checkbox("Object is in flow", &external);
  ImGui::SameLine();
  ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
  if (ImGui::InputInt("number of sides", &numSides)) {
    // Must have at least 3 sides
    if (numSides < 3) {
      numSides = 3;
      // Currently crashes if there are more than 17 sides
    } else if (numSides > 17) {
      numSides = 17;
    }
    // Set initial radius to 1 for number of sides
    rad = 1.0;
    // Set side length st radius is 1
    side = std::sqrt(2*(1-std::cos(M_PI*2/numSides)));
  }
  ImGui::InputFloat2("center", xc);
  if (ImGui::SliderFloat("side length", &side, 0.1f, 10.0f, "%.4f")) {
    rad = side/std::sqrt(2*(1-std::cos(M_PI*2/numSides)));
  }
  if (ImGui::SliderFloat("radius", &rad, 0.1f, 10.0f, "%.4f")) { side = rad*std::sqrt(2*(1-std::cos(M_PI*2/numSides))); }
  ImGui::SameLine();
  ShowHelpMarker("Polygons with different numbers of sides will appear similar in size with the same radius");
  ImGui::SliderFloat("orientation", &rotdeg, 0.0f, 359.0f, "%.0f");
  //ImGui::SliderAngle("orientation", &rotdeg);
  ImGui::TextWrapped("This feature will add a solid polygon boundary with n sides centered at the given coordinates");
  if (ImGui::Button("Add polygon boundary")) {
    m_external = external;
    m_x = xc[0];
    m_y = xc[1];
    m_theta = rotdeg;
    m_numSides = numSides;
    m_side = side;
    m_radius = rad;
    add = true;
    ImGui::CloseCurrentPopup();
  }

  return add;
}
#endif

void SolidPolygon::generate_draw_geom() {
  m_draw = init_elements(m_radius);
}

// node distribution with dense points at each end
double chebeshev_node(double a, double b, double k, double n) {
  return (a+b)*0.5+(b-a)*0.5*std::cos((2*(n-k)-1)*M_PI*0.5/n);
}

// node distribution with tunable density at each end
// first, a measure of how many panels are needed given densities at each end and relative panel size
size_t chebeshev_node_2_count(const double dens_left, const double dens_right, const double rel_pan_size) {
  // first, compute theta bounds given node densities
  assert(dens_left > 0.0 && dens_left < 1.0 && dens_right > 0.0 && dens_right < 1.0 && "Panel densities out of range");
  const double theta_left = M_PI - std::asin(dens_left);
  const double theta_right = std::asin(dens_right);
  // theta range
  const double theta_range = theta_left - theta_right;
  // x range from these
  const double x_range = std::cos(theta_right) - std::cos(theta_left);
  // now, size of the middle panel
  const double x_panel = x_range * rel_pan_size;
  // convert that to an angle
  const double theta_panel = 2.0 * std::asin(0.5*x_panel);
  // finally, find panel count
  return std::max((int)3, (int)(0.5 + theta_range/theta_panel));
}
// then the mathematics itself to generate the node positions
double chebeshev_node_2(const double dens_left, const double dens_right, const size_t k, const size_t n) {
  // first, compute theta bounds given node densities
  assert(dens_left > 0.0 && dens_left < 1.0 && dens_right > 0.0 && dens_right < 1.0 && "Panel densities out of range");
  const double theta_left = M_PI - std::asin(dens_left);
  const double theta_right = std::asin(dens_right);
  // theta range
  const double theta_range = theta_left - theta_right;
  // x range from these
  const double x_range = std::cos(theta_right) - std::cos(theta_left);

  // now, find test theta from input indices
  const double theta = theta_right + theta_range*k/(double)n;

  // finally, compute and scale the return value (from 0..1)
  return (std::cos(theta_right) - std::cos(theta)) / x_range;
}

// Create a NACA 4-digit airfoil
ElementPacket<float>
SolidAirfoil::init_elements(const float _ips) const {
  // If object has been removed, return no elements?
  if (not this->is_enabled()) return ElementPacket<float>();

  // created once
  // This should be equivalent to the num_panels param other boundaries use 
  const float m = m_maxCamber/100.0f;
  const float p = m_maxCambLoc/10.0f;
  const float t = m_thickness/100.0f;
  
  // number of panels on top surface
  //const size_t numX = std::ceil(m_chordLength*M_PI/_ips);
  const size_t numX = chebeshev_node_2_count(0.75, 0.1, _ips/m_chordLength);
  std::cout << "Creating NACA airfoil " << m_maxCamber << m_maxCambLoc << m_thickness << " with " << 2*numX << " panels" << std::endl;
  std::vector<float> x(4*numX);
  std::vector<Int> idx(4*numX);
  
  // first node (leading edge)
  x[0] = 0.0;
  x[1] = 0.0;
  // middle node (trailing edge)
  x[2*numX] = 1.0;
  x[2*numX+1] = 0.0;
  // formula thickness at TE (xol=1.0)
  const float yt_trail = (t/0.2)*(0.2969-0.1260-0.3516+0.2843-0.1015);

  // march along the chord, generating nodes and panels
  for (size_t i=1; i<=numX; i++) {
    //const float xol = chebeshev_node(0.0, 1.0, i, numX);
    const float xol = chebeshev_node_2(0.75, 0.1, i, numX);
    float yc;
    float dyc;
    if (xol < p) {
      yc = (m/std::pow(p,2))*(2*p*xol-std::pow(xol,2));
      dyc = (m/std::pow(p,2))*2*(p-xol);
    } else {
      yc = (m/std::pow(1-p,2))*(1-2*p+2*p*xol-std::pow(xol,2));
      dyc = (m/std::pow(1-p,2))*2*(p-xol);
    }
    float yt = (t/0.2)*(0.2969*std::sqrt(xol)-0.1260*xol-0.3516*std::pow(xol,2)+0.2843*std::pow(xol,3)-0.1015*std::pow(xol,4));
    // correct yt so that top and bottom meet at a single point
    yt -= xol*yt_trail;
    const float theta = std::tan(dyc);
    // add nodes
    if (i < numX) {
      // The top half
      x[2*i]   = xol-yt*std::sin(theta);
      x[2*i+1] = yc +yt*std::cos(theta);
      // The bottom half
      x[2*(2*numX-i)]   = xol+yt*std::sin(theta);
      x[2*(2*numX-i)+1] = yc -yt*std::cos(theta);
    }
    // Indices (arranged CW around the section because the it is solid)
    idx[2*(i-1)]   = i-1;
    idx[2*(i-1)+1] = i;
    idx[2*(2*numX-i)]   = 2*numX-i;
    idx[2*(2*numX-i)+1] = 2*numX-i+1;
  }

  // the last panel needs to point to the first node
  idx[4*numX-1] = 0;

  // scale, translate, and rotate into place
  const float st = std::sin(-M_PI * m_theta / 180.0);
  const float ct = std::cos(-M_PI * m_theta / 180.0);
  for (size_t i=0; i<2*numX; i++) {
    const float px = x[2*i];
    const float py = x[2*i+1];
    x[2*i]   = m_x + m_chordLength * (px*ct - py*st);
    x[2*i+1] = m_y + m_chordLength * (px*st + py*ct);
  }

  // val is bc, which is 0.0
  std::vector<float> val(idx.size()/2, 0.0);

  return ElementPacket<float>({x, idx, val});
}

void
SolidAirfoil::debug(std::ostream& os) const {
  os << to_string();
}

std::string
SolidAirfoil::to_string() const {
  std::stringstream ss;
  ss << " NACA " << m_maxCamber << m_maxCambLoc << m_thickness << " at (" << m_x << ", " << m_y << ") with chord length " << m_chordLength << " and " << m_theta << " aoa";
  return ss.str();
}

void
SolidAirfoil::from_json(const nlohmann::json j) {
  const std::vector<float> tr = j["translation"];
  m_x = tr[0];
  m_y = tr[1];
  m_maxCamber = j["maxCamber"];
  m_maxCambLoc = j["maxCambLoc"];
  m_thickness = j["thickness"];
  m_theta = j.value("rotation", 0.0);
  m_external = j.value("external", true);
  m_enabled = j.value("enabled",true);
  m_chordLength = j.value("chordLength", 1.0f);
}

nlohmann::json
SolidAirfoil::to_json() const {
  // make an object for the mesh
  nlohmann::json mesh = nlohmann::json::object();
  mesh["geometry"] = "airfoil";
  mesh["translation"] = {m_x, m_y};
  mesh["maxCamber"] = m_maxCamber;
  mesh["maxCambLoc"] = m_maxCambLoc;
  mesh["thickness"] = m_thickness;
  mesh["rotation"] = m_theta;
  mesh["chordLength"] = m_chordLength;
  mesh["external"] = m_external;
  mesh["enabled"] = m_enabled;
  return mesh;
}

#ifdef USE_IMGUI
bool SolidAirfoil::draw_info_gui() {
  static bool external = m_external;
  static float xc[2] = {m_x, m_y};
  static float rotdeg = m_theta;
  static std::string naca = std::to_string(m_maxCamber)+std::to_string(m_maxCambLoc)+std::to_string(m_thickness);
  static float chordLength = m_chordLength;
  bool add = false;
 
  ImGui::Checkbox("Object is in flow", &external);
  ImGui::SameLine();
  ShowHelpMarker("Keep checked if object is immersed in flow,\nuncheck if flow is inside of object");
  if (ImGui::InputText("NACA 4-digit Number", &naca)) {
    if (naca.size() > 4) {
        naca = naca.substr(0, 4);
    }
  }
  ImGui::InputFloat2("Leading Edge Position", xc);
  ImGui::SliderFloat("Angle of Attack", &rotdeg, -180.0f, 180.0f, "%.0f");
  ImGui::SliderFloat("Chord Length", &chordLength, 0.1f, 5.0f, "%.1f");
  ImGui::TextWrapped("This feature will add a NACA airfoil with LE at the given coordinates");
  if (ImGui::Button("Add NACA airfoil")) {
    m_external = external;
    m_x = xc[0];
    m_y = xc[1];
    m_theta = rotdeg;
    m_maxCamber = std::stoi(naca.substr(0, 1));
    m_maxCambLoc = std::stoi(naca.substr(1, 1));
    m_thickness = std::stoi(naca.substr(2, 2));
    m_chordLength = chordLength;
    add = true;
    ImGui::CloseCurrentPopup();
  }

  return add;
}
#endif

void SolidAirfoil::generate_draw_geom() {
  m_draw = init_elements(m_chordLength/20.0);
}
