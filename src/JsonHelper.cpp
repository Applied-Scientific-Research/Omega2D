/*
 * JsonHelper.cpp - Methods to help with reading and writing JSON files
 *
 * (c)2019-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
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

#include "JsonHelper.h"
#include "Body.h"

#include "json/json.hpp"

#include <string>
#include <vector>
#include <iostream>	// for cout,endl
#include <iomanip>	// for setw
#include <fstream>	// for ofstream
#include <exception>

using json = nlohmann::json;
 

//
// create a json object from the file and apply it to the simulation
//
// We should load the sim in a separate function. We pass values by reference and maybe have it
// return a bool to indicate success (that might be overkill)
json read_json (const std::string filename) {
  // read a JSON file
  std::cout << "  loading simulation " << filename << std::endl;
  std::ifstream json_in(filename);
  json j;
  // TODO: this next line is where we crash if the file is not present!
  try {
    json_in >> j;
  } catch (...) {
    std::cout << "The json file " << filename << " does not exist." << std::endl;
  }
  return j;
}

void parse_json(Simulation& sim,
                std::vector<std::unique_ptr<FlowFeature>>& ffeatures,
                std::vector<std::unique_ptr<BoundaryFeature>>& bfeatures,
                std::vector<std::unique_ptr<MeasureFeature>>& mfeatures,
                RenderParams& rp,
                nlohmann::json j) {

  // march through the parameters and apply them
  if (j.count("description") == 1) {
    std::string this_desc = j["description"];
    sim.set_description(this_desc);
  }

  if (j.count("version") == 1) {
    json params = j["version"];
    if (params.find("Omega2D") != params.end()) {
      int o2dv = params["Omega2D"];
      std::cout << "  is an Omega2D file, version " << o2dv << std::endl;
    }
    if (params.find("Omega3D") != params.end()) {
      std::cout << "  is an Omega3D file" << std::endl;
      throw "Input file is not compatible.";
    }
    if (params.find("OmegaFlow") != params.end()) {
      std::cout << "  is an OmegaFlow file" << std::endl;
      throw "Input file is not compatible.";
    }
    if (params.find("jsonInput") != params.end()) {
      int jiv = params["jsonInput"];
      std::cout << "  json input version " << jiv << std::endl;
    }
  }

  // must do this first, as we need to set viscous before reading Re
  if (j.count("simparams") == 1) {
    sim.from_json(j["simparams"]);
  }

  // now we can read Re, freestream
  if (j.count("flowparams") == 1) {
    sim.flow_from_json(j["flowparams"]);
  }

  // finish up with runtime stuff
  if (j.count("runtime") == 1) {
    sim.runtime_from_json(j["runtime"]);
  }

  // ask RenderParams to read itself
  if (j.count("drawparams") == 1) {
    rp.from_json(j["drawparams"]);
  }

  // read the flow features, if any
  if (j.count("flowstructures") == 1) {
    //std::cout << "  found flowstructures" << std::endl;
    const std::vector<json> ff_json = j["flowstructures"];

    // iterate through vector of flow features
    for (auto const& ff : ff_json) {
      // pass ff into a function in FlowFeature to generate the object
      parse_flow_json(ffeatures, ff);
    }
  }

  // read the boundary features, if any
  if (j.count("bodies") == 1) {
    //std::cout << "  found bodies" << std::endl;
    std::vector<json> bdy_json = j["bodies"];

    // iterate through vector of bodies and generate boundary features
    for (auto const& bdy : bdy_json) {
      // make a new Body (later) and attach these geometries to it
      //std::cout << "  for this body" << std::endl;
      auto bp = std::make_shared<Body>();

      // get the body name ("ground" is default)
      if (bdy.find("name") != bdy.end()) {
        // by making a variable out of it, we implicitly cast whatever is there to a string
        const std::string bname = bdy["name"];
        // so we don't screw this part up
        bp->set_name(bname);
      }

      // get the parent name, if any ("none" is default)
      if (bdy.find("parent") != bdy.end()) {
        const std::string pname = bdy["parent"];
        bp->set_parent_name(pname);
      }

      // get the motions of this body w.r.t. parent
      if (bdy.find("translation") != bdy.end()) {
        // look for an array of possibly dissimilar terms
        const std::vector<json> trans_json = bdy["translation"];
        size_t idx = 0;
        for (auto const& term : trans_json) {
          // each entry can be a float or a string
          if (term.is_number()) {
            const double val = term.get<double>();
            //std::cout << "    component " << idx << " of position is constant (" << val << ")" << std::endl;
            bp->set_pos(idx, val);
          } else if (term.is_string()) {
            const std::string expr = term.get<std::string>();
            //std::cout << "    component " << idx << " of position is expression (" << expr << ")" << std::endl;
            bp->set_pos(idx, expr);
          }
          ++idx;
        }
      }
      if (bdy.find("rotation") != bdy.end()) {
        // look for a float or a string
        if (bdy["rotation"].is_number()) {
          const double val = bdy["rotation"].get<double>();
          bp->set_rot(val);
        } else if (bdy["rotation"].is_string()) {
          const std::string expr = bdy["rotation"].get<std::string>();
          bp->set_rot(expr);
        }
      }

      // see if there are meshes (there don't have to be - a Body can just act as a virtual joint
      if (bdy.count("meshes") == 1) {
        //std::cout << "  found the meshes" << std::endl;
        std::vector<json> bf_json = bdy["meshes"];

        // iterate through vector of meshes, all on this body
        for (auto const& bf : bf_json) {
          // pass bf into a function in BoundaryFeature to generate the object
          parse_boundary_json(bfeatures, bp, bf);

          // this turns off augmentation!
          //parse_boundary_json(bfeatures, nullptr, bf);
        }
      }

      // and add the Body pointer bp to the master list
      sim.add_body(bp);
    }
  }

  // read the measurement features, if any
  if (j.count("measurements") == 1) {
    //std::cout << "  found measurements" << std::endl;
    const std::vector<json> mf_json = j["measurements"];

    // iterate through vector of measurement features
    for (auto const& mf : mf_json) {
      // pass mf into a function in MeasureFeature to generate the object
      parse_measure_json(mfeatures, mf);
    }
  }

  for (auto const& ff : ffeatures) { ff->generate_draw_geom(); }
  for (auto const& bf : bfeatures) { bf->generate_draw_geom(); }
  for (auto const& mf : mfeatures) { mf->generate_draw_geom(); }
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
                const RenderParams& rp,
                const std::string filename) {

  json j;

  const std::string desc = sim.get_description();
  if (desc.empty()) {
    j["description"] = "Simulation created by Omega2D";
  } else {
    j["description"] = desc;
  }

  j["version"] = { {"Omega2D", 1}, {"jsonInput", 1} };

  j["runtime"] = sim.runtime_to_json();

  j["flowparams"] = sim.flow_to_json();

  j["simparams"] = sim.to_json();

  j["drawparams"] = rp.to_json();

  if (ffeatures.size() > 0) {
    // assemble a vector of flow features
    std::vector<json> jflows;
    for (auto const& ff : ffeatures) {
      jflows.push_back(ff->to_json());
    }
    j["flowstructures"] = jflows;
  }

  // assemble a vector of bodies, each with 0 or more boundary features
  std::vector<json> jbods;
  auto b_iter = sim.bodies_begin();
  for (; b_iter != sim.bodies_end(); ++b_iter) {

    // generate the body information
    auto bp = *b_iter;
    nlohmann::json jb = bp->to_json();
    //std::cout << "FOR BODY (" << (*b_iter)->get_name() << ")" << std::endl;

    // meshes has to be an array
    nlohmann::json meshes = nlohmann::json::array();
    // now add all accompanying geometries (if any)
    for (auto const& bf : bfeatures) {
      //std::cout << "  BDRY FEATURE (" << (&bf) << ")";
      if (bf->get_body() == bp) {
        //std::cout << " MATCHES" << std::endl;
        meshes.push_back(bf->to_json());
      } else {
        //std::cout << " DOES NOT MATCH" << std::endl;
      }
    }
    jb["meshes"] = meshes;
    jbods.push_back(jb);
  }
  j["bodies"] = jbods;

  if (mfeatures.size() > 0) {
    // assemble a vector of measurement features
    std::vector<json> jmeas;
    for (auto const& mf : mfeatures) {
      jmeas.push_back(mf->to_json());
    }
    j["measurements"] = jmeas;
  }

  // write prettified JSON to the given file
  std::ofstream json_out(filename);
  json_out << std::setw(2) << j << std::endl;
}

