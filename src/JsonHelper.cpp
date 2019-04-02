/*
 * JsonHelper.cpp - Methods to help with reading and writing JSON files
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
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
void read_json (Simulation& sim,
                std::vector<std::unique_ptr<FlowFeature>>& ffeatures,
                std::vector<std::unique_ptr<BoundaryFeature>>& bfeatures,
                std::vector<std::unique_ptr<MeasureFeature>>& mfeatures,
                RenderParams& rp,
                const std::string filename) {

  // read a JSON file
  std::ifstream json_in(filename);
  json j;
  json_in >> j;

  std::cout << std::endl << "Loading simulation from " << filename << std::endl;

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

  if (j.count("runtime") == 1) {
    json params = j["runtime"];
    // no runtime parameters used yet
    if (params.find("statusFile") != params.end()) {
      std::string sfile = params["statusFile"];
      sim.set_status_file_name(sfile);
      std::cout << "  status file name= " << sfile << std::endl;
    }
  }

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
    if (params.find("outputDt") != params.end()) {
      float odt = params["outputDt"];
      sim.set_output_dt(odt);
      std::cout << "  setting output_dt= " << odt << std::endl;
    }
    if (params.find("endTime") != params.end()) {
      float et = params["endTime"];
      sim.set_end_time(et);
      std::cout << "  setting end_time= " << et << std::endl;
    }
    if (params.find("maxSteps") != params.end()) {
      size_t ms = params["maxSteps"];
      sim.set_max_steps(ms);
      std::cout << "  setting max_steps= " << ms << std::endl;
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
      // eventually support an expression for Uinf instead of just a single float
      std::vector<float> new_fs = {0.0, 0.0, 0.0};
      new_fs.resize(Dimensions);
      if (params["Uinf"].is_array()) {
        new_fs = params["Uinf"].get<std::vector<float>>();
      } else if (params["Uinf"].is_number()) {
        new_fs[0] = params["Uinf"].get<float>();
      }
      float* fs = sim.addr_fs();
      for (size_t i=0; i<Dimensions; ++i) fs[i] = new_fs[i];
      std::cout << "  setting freestream to " << fs[0] << " " << fs[1] << std::endl;
    }
  }

  // must do this first, as we need to set viscous before reading Re
  if (j.count("drawparams") == 1) {
    json params = j["drawparams"];

    if (params.find("backgroundColor") != params.end()) {
      std::vector<float> new_color = params["backgroundColor"];
      //std::cout << "  setting bg color to " << new_color[0] << " " << new_color[1] << " " << new_color[2] << " " << new_color[3] << std::endl;
      for (size_t i=0; i<4; ++i) rp.clear_color[i] = new_color[i];
    }
    if (params.find("featureColor") != params.end()) {
      std::vector<float> new_color = params["featureColor"];
      //std::cout << "  setting def color to " << new_color[0] << " " << new_color[1] << " " << new_color[2] << " " << new_color[3] << std::endl;
      for (size_t i=0; i<4; ++i) rp.default_color[i] = new_color[i];
    }
    if (params.find("positiveColor") != params.end()) {
      std::vector<float> new_color = params["positiveColor"];
      //std::cout << "  setting pos color to " << new_color[0] << " " << new_color[1] << " " << new_color[2] << " " << new_color[3] << std::endl;
      for (size_t i=0; i<4; ++i) rp.pos_circ_color[i] = new_color[i];
    }
    if (params.find("negativeColor") != params.end()) {
      std::vector<float> new_color = params["negativeColor"];
      //std::cout << "  setting neg color to " << new_color[0] << " " << new_color[1] << " " << new_color[2] << " " << new_color[3] << std::endl;
      for (size_t i=0; i<4; ++i) rp.neg_circ_color[i] = new_color[i];
    }
    if (params.find("density") != params.end()) {
      rp.circ_density = params["density"];
    }
    if (params.find("tracerScale") != params.end()) {
      rp.tracer_scale = params["tracerScale"];
    }
    if (params.find("viewPoint") != params.end()) {
      std::vector<float> new_vec = params["viewPoint"];
      rp.vcx = new_vec[0];
      rp.vcy = new_vec[1];
    }
    if (params.find("viewScale") != params.end()) {
      rp.vsize = params["viewScale"];
    }
    if (params.find("windowSize") != params.end()) {
      std::vector<int> new_vec = params["windowSize"];
      rp.width = new_vec[0];
      rp.height = new_vec[1];
    }
  }

  // read the flow features, if any
  if (j.count("flowstructures") == 1) {
    //std::cout << "  found flowstructures" << std::endl;
    const std::vector<json> ff_json = j["flowstructures"];

    // iterate through vector of flow features
    for (auto const& ff: ff_json) {
      // pass ff into a function in FlowFeature to generate the object
      parse_flow_json(ffeatures, ff);
    }
  }

  // read the boundary features, if any
  if (j.count("bodies") == 1) {
    //std::cout << "  found bodies" << std::endl;
    std::vector<json> bdy_json = j["bodies"];

    // iterate through vector of bodies and generate boundary features
    for (auto const& bdy: bdy_json) {
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
        for (auto const& term: trans_json) {
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
        for (auto const& bf: bf_json) {
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
    for (auto const& mf: mf_json) {
      // pass mf into a function in MeasureFeature to generate the object
      parse_measure_json(mfeatures, mf);
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

  const std::string sfile = sim.get_status_file_name();
  if (not sfile.empty()) {
    j["runtime"] = { {"statusFile", sfile} };
  }

  float re = *(sim.addr_re());
  float* fs = sim.addr_fs();
  j["flowparams"] = { {"Re", re},
                      {"Uinf", {fs[0], fs[1]} } };

  float dt = *(sim.addr_dt());
  std::string viscous = sim.get_diffuse() ? "vrm" : "none";
  j["simparams"] = { {"nominalDt", dt},
                     {"viscous", viscous} };
  if (sim.using_max_steps()) j["simparams"].push_back( {"maxSteps", sim.get_max_steps()} );
  if (sim.using_end_time()) j["simparams"].push_back( {"endTime", sim.get_end_time()} );

  j["drawparams"] = { {"backgroundColor", {rp.clear_color[0], rp.clear_color[1], rp.clear_color[2], rp.clear_color[3]} },
                      {"featureColor", {rp.default_color[0], rp.default_color[1], rp.default_color[2], rp.default_color[3]} },
                      {"positiveColor", {rp.pos_circ_color[0], rp.pos_circ_color[1], rp.pos_circ_color[2], rp.pos_circ_color[3]} },
                      {"negativeColor", {rp.neg_circ_color[0], rp.neg_circ_color[1], rp.neg_circ_color[2], rp.neg_circ_color[3]} },
                      {"density", rp.circ_density},
                      {"tracerScale", rp.tracer_scale},
                      {"viewPoint", {rp.vcx, rp.vcy} },
                      {"viewScale", rp.vsize},
                      {"windowSize", {rp.width, rp.height} } };

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
    for (auto const& bf: bfeatures) {
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

