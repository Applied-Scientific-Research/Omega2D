/*
 * JsonHelper.cpp - Methods to help with reading and writing JSON files
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#include "JsonHelper.h"

#include "json/json.hpp"

//#include <string>
#include <vector>
#include <iostream>	// for cout,endl
#include <iomanip>	// for setw
#include <fstream>	// for ofstream

using json = nlohmann::json;
 

//
// create a json object from the file and apply it to the simulation
//
void read_json (Simulation& sim,
                std::vector<std::unique_ptr<FlowFeature>>& ffeatures,
                std::vector<std::unique_ptr<BoundaryFeature>>& bfeatures,
                std::vector<std::unique_ptr<MeasureFeature>>& mfeatures,
                const std::string filename) {

  // read a JSON file
  std::ifstream json_in(filename);
  json j;
  json_in >> j;

  // march through the parameters and apply them

  // must do this first, as we need to set viscous before reading Re
  if (j.count("simparams") == 1) {
    json params = j["simparams"];
    if (params.find("nominalDt") != params.end()) {
      float* dt = sim.addr_dt();
      *dt = params["nominalDt"];
    }
    if (params.find("viscous") != params.end()) {
      std::string viscous = params["viscous"];
      if (viscous == "vrm") { sim.set_diffuse(true); }
      else { sim.set_diffuse(false); }
    }
  }

  // now we can read Re
  if (j.count("flowparams") == 1) {
    json params = j["flowparams"];
    if (params.find("Re") != params.end()) {
      float* re = sim.addr_re();
      *re = params["Re"];
    }
    if (params.find("Uinf") != params.end()) {
      float* fs = sim.addr_fs();
      std::vector<float> new_fs = params["Uinf"];
      fs[0] = new_fs[0];
      fs[1] = new_fs[1];
    }
  }

  // read the flow features, if any
  if (j.count("flowstructures") == 1) {
    json params = j["flowstructures"];
    //std::vector<json> = 
  }
}


//
// create a json object representing the simulation and write it to the given file
//
// the inputs should be const
//
void write_json(Simulation& sim,
                std::vector<std::unique_ptr<FlowFeature>>& ffeatures,
                std::vector<std::unique_ptr<BoundaryFeature>>& bfeatures,
                std::vector<std::unique_ptr<MeasureFeature>>& mfeatures,
                const std::string filename) {

  json j;

  j["description"] = "Simulation created by Omega2D";

  j["version"] = { {"Omega2D", 1}, {"jsonInput", 1} };

  float re = *(sim.addr_re());
  float* fs = sim.addr_fs();
  j["flowparams"] = { {"Re", re},
                      {"Uinf", {fs[0], fs[1]} } };

  float dt = *(sim.addr_dt());
  std::string viscous = sim.get_diffuse() ? "vrm" : "none";
  j["simparams"] = { {"nominalDt", dt},
                     {"viscous", viscous } };

  // assemble a vector of flow features
  std::vector<json> jflows;
  for (auto const& ff: ffeatures) {
    jflows.push_back(ff->to_json());
  }
  j["flowstructures"] = jflows;

  // assemble a vector of boundary features
  std::vector<json> jbounds;
  for (auto const& bf: bfeatures) {
    jbounds.push_back(bf->to_json());
  }
  j["bodies"] = jbounds;

  // assemble a vector of measurement features
  std::vector<json> jmeas;
  for (auto const& mf: mfeatures) {
    jmeas.push_back(mf->to_json());
  }
  j["measurements"] = jmeas;

  // write prettified JSON to the given file
  std::ofstream json_out(filename);
  json_out << std::setw(2) << j << std::endl;
}

