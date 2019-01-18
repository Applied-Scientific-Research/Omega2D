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
// see if the library works
//
void read_json_test() {

  // read a JSON file
  std::ifstream json_in("file.json");
  json j_object;
  json_in >> j_object;

  // create JSON values
  //json j_object = {{"one", 1}, {"two", 2}};
  //json j_array = {1, 2, 4, 8, 16};
 
  // example for an object
  for (auto& x : j_object.items()) {
    std::cout << "key: " << x.key() << ", value: " << x.value() << '\n';
  }
 
  // example for an array
  //for (auto& x : j_array.items()) {
  //  std::cout << "key: " << x.key() << ", value: " << x.value() << '\n';
  //}
}

//
// create a json object and write it to a test file
//
void write_json_test() {

  // create an empty structure (null)
  json j;

  // add a number that is stored as double (note the implicit conversion of j to an object)
  j["pi"] = 3.141;

  // add a Boolean that is stored as bool
  j["happy"] = true;

  // add a string that is stored as std::string
  j["name"] = "Niels";

  // add another null object by passing nullptr
  j["nothing"] = nullptr;

  // add an object inside the object
  j["answer"]["everything"] = 42;

  // add an array that is stored as std::vector (using an initializer list)
  j["list"] = { 1, 0, 2 };

  // add another object (using an initializer list of pairs)
  j["object"] = { {"currency", "USD"}, {"value", 42.99} };

  // instead, you could also write (which looks very similar to the JSON above)
  json j2 = {
    {"pi", 3.141},
    {"happy", true},
    {"name", "Niels"},
    {"nothing", nullptr},
    {"answer", {
      {"everything", 42}
    }},
    {"list", {1, 0, 2}},
    {"object", {
      {"currency", "USD"},
      {"value", 42.99}
    }}
  };

  // write prettified JSON to another file
  std::ofstream json_out("pretty.json");
  json_out << std::setw(2) << j << std::endl;
}

//
// create a json object representing the simulation and write it to the given file
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

  //jflows.push_back({{"name", "myname"}, {"value", 3.0}});
  //jflows.push_back({{"name", "othername"}, {"value", 13.0}});
  //j["testing"] = jflows;

  // assemble a vector of boundary features
  std::vector<json> jbounds;
  for (auto const& bf: bfeatures) {
    jbounds.push_back(bf->to_json());
  }
  j["bodies"] = jbounds;

  // assemble a vector of measurement features
  //std::vector<json> jmeas;
  //j["measurements"] = jmeas;

  // write prettified JSON to the given file
  std::ofstream json_out(filename);
  json_out << std::setw(2) << j << std::endl;
}

