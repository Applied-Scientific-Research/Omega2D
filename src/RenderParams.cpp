/*
 * RenderParams.cpp - Methods for rendering/drawing parameters
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#include "RenderParams.h"

#include <vector>
//#include <iostream>


// set parameters from a json object
void
RenderParams::from_json(const nlohmann::json j) {

  if (j.find("backgroundColor") != j.end()) {
    std::vector<float> new_color = j["backgroundColor"];
    for (size_t i=0; i<4; ++i) clear_color[i] = new_color[i];
  }
  if (j.find("featureColor") != j.end()) {
    std::vector<float> new_color = j["featureColor"];
    for (size_t i=0; i<4; ++i) default_color[i] = new_color[i];
  }
  if (j.find("positiveColor") != j.end()) {
    std::vector<float> new_color = j["positiveColor"];
    for (size_t i=0; i<4; ++i) pos_circ_color[i] = new_color[i];
  }
  if (j.find("negativeColor") != j.end()) {
    std::vector<float> new_color = j["negativeColor"];
    for (size_t i=0; i<4; ++i) neg_circ_color[i] = new_color[i];
  }

  if (j.find("density") != j.end()) {
    circ_density = j["density"];
  }
  if (j.find("tracerScale") != j.end()) {
    tracer_scale = j["tracerScale"];
  }
  if (j.find("viewPoint") != j.end()) {
    std::vector<float> new_vec = j["viewPoint"];
    vcx = new_vec[0];
    vcy = new_vec[1];
  }
  if (j.find("viewScale") != j.end()) {
    vsize = j["viewScale"];
  }
  if (j.find("windowSize") != j.end()) {
    std::vector<int> new_vec = j["windowSize"];
    width = new_vec[0];
    height = new_vec[1];
  }

}


// create and write a json object for this object's members
nlohmann::json
RenderParams::to_json() const {
  nlohmann::json j;

  j["density"] = circ_density;
  j["tracerScale"] = tracer_scale;
  j["viewPoint"] = {vcx, vcy};
  j["viewScale"] = vsize;
  j["windowSize"] = {width, height};

  // colors have to be arrays - probably an easier way to do this
  nlohmann::json jcol = nlohmann::json::array();
  for (size_t i=0; i<4; ++i) { jcol.push_back(clear_color[i]); }
  j["backgroundColor"] = jcol;

  for (size_t i=0; i<4; ++i) { jcol[i] = default_color[i]; }
  j["featureColor"] = jcol;

  for (size_t i=0; i<4; ++i) { jcol[i] = pos_circ_color[i]; }
  j["positiveColor"] = jcol;

  for (size_t i=0; i<4; ++i) { jcol[i] = neg_circ_color[i]; }
  j["negativeColor"] = jcol;

  for (size_t i=0; i<4; ++i) { jcol[i] = default_color[i]; }
  j["featureColor"] = jcol;

  return j;
}

