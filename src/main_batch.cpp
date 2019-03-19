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
  size_t nsteps = 0;
  RenderParams rparams;

  // a string to hold any error messages
  std::string sim_err_msg;

  // load a simulation from a JSON file - check command line for file name
  if (argc == 2) {
    std::string infile = argv[1];
    read_json(sim, ffeatures, bfeatures, mfeatures, rparams, infile);
    std::cout << std::endl << "Loaded simulation from " << infile << std::endl;
  } else {
    std::cout << std::endl << "Usage:" << std::endl;
    std::cout << "  " << argv[0] << " filename.json" << std::endl << std::endl;
    return -1;
  }


  std::cout << std::endl << "Initializing simulation" << std::endl;

  // initialize particle distributions
  for (auto const& ff: ffeatures) {
    sim.add_particles( ff->init_particles(sim.get_ips()) );
  }

  // initialize solid objects
  for (auto const& bf : bfeatures) {
    sim.add_boundary( bf->get_body(), bf->init_elements(sim.get_ips()) );
  }

  // initialize measurement features
  for (auto const& mf: mfeatures) {
    sim.add_fldpts( mf->init_particles(0.1*sim.get_ips()), mf->moves() );
  }

  sim.set_initialized();

  //
  // Main loop
  //

  while (true) {

    // check flow for blow-up or errors
    sim_err_msg = sim.check_simulation(ffeatures.size(), bfeatures.size());

    if (sim_err_msg.empty()) {
      // the last simulation step was fine, OK to continue

      // generate new particles from emitters
      for (auto const& ff : ffeatures) {
        sim.add_particles( ff->step_particles(sim.get_ips()) );
      }
      for (auto const& mf: mfeatures) {
        sim.add_fldpts( mf->step_particles(0.1*sim.get_ips()), mf->moves() );
      }

      // begin a new dynamic step: convection and diffusion
      sim.step();

    } else {
      // the last step had some difficulty
      std::cout << std::endl << "ERROR: " << sim_err_msg;

      // stop the run
      break;
    }

    nsteps++;

    // for testing: always break after a few steps
    //if (nsteps == 2) break;

    // check vs. stopping conditions
    if (sim.using_max_steps() and sim.get_max_steps() == nsteps) break;
    if (sim.using_end_time() and sim.get_end_time() <= sim.get_time()) break;

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

