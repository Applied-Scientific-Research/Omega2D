/*
 * main_batch.cpp - Driver code for Omega2D + Vc vortex particle method
 *                  and boundary element method solver, batch version
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "Simulation.h"
#include "JsonHelper.h"
#include "RenderParams.h"

#ifdef _WIN32
  // for glad
  #define APIENTRY __stdcall
  // for C++11 stuff
  #include <ciso646>
#endif

#include <iostream>
#include <vector>


// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << std::endl << "Omega2D Batch" << std::endl;

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;
  std::vector< std::unique_ptr<MeasureFeature> > mfeatures;
  RenderParams rparams;

  // a string to hold any error messages
  std::string sim_err_msg;

  // load a simulation from a JSON file - check command line for file name
  if (argc == 2) {
    std::string infile = argv[1];
    nlohmann::json j = read_json(infile);
    parse_json(sim, ffeatures, bfeatures, mfeatures, rparams, j);
  } else {
    std::cout << std::endl << "Usage:" << std::endl;
    std::cout << "  " << argv[0] << " filename.json" << std::endl << std::endl;
    return -1;
  }


  std::cout << std::endl << "Initializing simulation" << std::endl;

  // initialize particle distributions
  for (auto const& ff: ffeatures) {
    //if (ff->is_enabled()) sim.add_particles( ff->init_particles(sim.get_ips()) );
    if (ff->is_enabled()) {
      ElementPacket<float> newpacket = ff->init_elements(sim.get_ips());
      // echo any errors
      /*if (good)*/ sim.add_elements( newpacket, active, lagrangian, ff->get_body() );
    }
    //if (ff->is_enabled()) sim.add_elements( ff->init_elements(sim.get_ips()), active, lagrangian, ff->get_body() );
  }

  // initialize solid objects
  for (auto const& bf : bfeatures) {
    //if (bf->is_enabled()) sim.add_boundary( bf->get_body(), bf->init_elements(sim.get_ips()) );
    if (bf->is_enabled()) {
      ElementPacket<float> newpacket = bf->init_elements(sim.get_ips());
      const move_t newmovetype = (bf->get_body() ? bodybound : fixed);
      sim.add_elements(newpacket, reactive, newmovetype, bf->get_body() );
    }
  }

  // initialize measurement features
  for (auto const& mf: mfeatures) {
    //if (mf->is_enabled()) sim.add_fldpts( mf->init_particles(rparams.tracer_scale*sim.get_ips()), mf->moves() );
    if (mf->is_enabled()) {
      const move_t newMoveType = (mf->get_is_lagrangian() ? lagrangian : fixed);
      sim.add_elements( mf->init_elements(rparams.tracer_scale*sim.get_ips()), inert, newMoveType, mf->get_body() );
    }
  }

  sim.set_initialized();

  // check init for blow-up or errors
  sim_err_msg = sim.check_initialization();

  if (not sim_err_msg.empty()) {
    // the initialization had some difficulty
    std::cout << std::endl << "ERROR: " << sim_err_msg;
    // stop the run
    return 1;
  }


  //
  // Main loop
  //

  while (true) {

    // check flow for blow-up or errors
    sim_err_msg = sim.check_simulation();

    if (sim_err_msg.empty()) {
      // the last simulation step was fine, OK to continue
      // generate new particles from emitters
      for (auto const& ff: ffeatures) {
        if (ff->is_enabled()) {
          ElementPacket<float> newpacket = ff->step_elements(sim.get_ips());
          // echo any errors
           sim.add_elements( newpacket, active, lagrangian, ff->get_body() );
        }
      }

      for (auto const& mf: mfeatures) {
        //if (mf->is_enabled()) sim.add_fldpts( mf->init_particles(rparams.tracer_scale*sim.get_ips()), mf->moves() );
        if (mf->is_enabled()) {
          const move_t newMoveType = (mf->get_is_lagrangian() ? lagrangian : fixed);
          sim.add_elements( mf->step_elements(rparams.tracer_scale*sim.get_ips()), inert, newMoveType, mf->get_body() );
      }

      // begin a new dynamic step: convection and diffusion
      sim.step();
    } else {
      // the last step had some difficulty
      std::cout << std::endl << "ERROR: " << sim_err_msg;

      // stop the run
      break;
    }

    // export data files at this step?

    // check vs. stopping conditions
    if (sim.test_vs_stop()) break;

  } // end step


  // Save final step if desired
  if (false) {
    // just make up a file name and write it
    std::string outfile = "output.json";
    // write and echo
    write_json(sim, ffeatures, bfeatures, mfeatures, rparams, outfile);
    std::cout << std::endl << "Wrote simulation to " << outfile << std::endl;
  }

  sim.reset();
  std::cout << "Quitting" << std::endl;

  return 0;
}

