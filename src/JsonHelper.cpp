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

  std::cout << std::endl << "Loading simulation from " << filename << std::endl;

  // march through the parameters and apply them

  // must do this first, as we need to set viscous before reading Re
  if (j.count("simparams") == 1) {
    json params = j["simparams"];
    if (params.find("nominalDt") != params.end()) {
      float* dt = sim.addr_dt();
      *dt = params["nominalDt"];
      std::cout << "  setting dt= " << *dt << std::endl;
    }
    if (params.find("viscous") != params.end()) {
      std::string viscous = params["viscous"];
      if (viscous == "vrm") { sim.set_diffuse(true); }
      else { sim.set_diffuse(false); }
      std::cout << "  setting is_viscous= " << sim.get_diffuse() << std::endl;
    }
  }

  // now we can read Re
  if (j.count("flowparams") == 1) {
    json params = j["flowparams"];
    if (params.find("Re") != params.end()) {
      float* re = sim.addr_re();
      *re = params["Re"];
      std::cout << "  setting re= " << *re << std::endl;
    }
    if (params.find("Uinf") != params.end()) {
      float* fs = sim.addr_fs();
      std::vector<float> new_fs = params["Uinf"];
      fs[0] = new_fs[0];
      fs[1] = new_fs[1];
      std::cout << "  setting freestream to " << fs[0] << " " << fs[1] << std::endl;
    }
  }

  // read the flow features, if any
  if (j.count("flowstructures") == 1) {
    std::vector<json> ff_json = j["flowstructures"];
    // iterate through vector of flow features
    for (auto const& ff: ff_json) {
      //if (ff.count("type") == 1) {
      std::string ftype = ff["type"];
      if (ftype == "single particle") {
        std::cout << "  found single particle" << std::endl;
        const float str = ff["strength"];
        const std::vector<float> c = ff["center"];
        ffeatures.emplace_back(std::make_unique<SingleParticle>(c[0], c[1], str));
      } else if (ftype == "vortex blob") {
        std::cout << "  found vortex blob" << std::endl;
        const float rad = ff["radius"];
        const float soft = ff["softness"];
        const float str = ff["strength"];
        const std::vector<float> c = ff["center"];
        ffeatures.emplace_back(std::make_unique<VortexBlob>(c[0], c[1], str, rad, soft));
      } else if (ftype == "block of random") {
        std::cout << "  found block of random" << std::endl;
        const std::vector<float> c = ff["center"];
        const std::vector<float> sz = ff["size"];
        const std::vector<float> sr = ff["strength range"];
        const int num = ff["num"];
        ffeatures.emplace_back(std::make_unique<BlockOfRandom>(c[0], c[1], sz[0], sz[1], sr[0], sr[1], num));
      } else if (ftype == "particle emitter") {
        std::cout << "  found particle emitter" << std::endl;
        const float str = ff["strength"];
        const std::vector<float> c = ff["center"];
        ffeatures.emplace_back(std::make_unique<ParticleEmitter>(c[0], c[1], str));
      }
    }
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

