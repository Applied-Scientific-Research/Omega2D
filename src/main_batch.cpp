/*
 * main_batch.cpp - Driver code for Omega2D + Vc vortex particle method
 *                  and boundary element method solver, batch version
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
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

#include "FlowFeature.h"
#include "BoundaryFeature.h"
#include "MeasureFeature.h"
#include "Simulation.h"
#include "JsonHelper.h"
#include "RenderParams.h"

#ifdef _WIN32
  // for glad
  #ifndef APIENTRY
    #define APIENTRY __stdcall
  #endif
  // for C++11 stuff that Windows can't get right
  #include <ciso646>
#endif

#include <iostream>
#include <vector>


// execution starts here

int main(int argc, char const *argv[]) {

  std::cout << std::endl << "Omega2D Batch" << std::endl;
  if (VERBOSE) { std::cout << "  VERBOSE is on" << std::endl; }

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
    if (ff->is_enabled()) {
      ElementPacket<float> newpacket = ff->init_elements(sim.get_ips());
      sim.add_elements( newpacket, active, lagrangian, ff->get_body() );
    }
  }

  // initialize solid objects
  for (auto const& bf : bfeatures) {
    if (bf->is_enabled()) {
      ElementPacket<float> newpacket = bf->init_elements(sim.get_ips());
      const move_t newmovetype = (bf->get_body() ? bodybound : fixed);
      sim.add_elements(newpacket, reactive, newmovetype, bf->get_body() );
    }
  }

  // initialize measurement features
  for (auto const& mf: mfeatures) {
    if (mf->is_enabled()) {
      const move_t newMoveType = (mf->get_is_lagrangian() ? lagrangian : fixed);
      sim.add_elements( mf->init_elements(rparams.tracer_scale*sim.get_ips()), inert, newMoveType, mf->get_body() );
    }
  }

  // initialize hybrid features
  for (auto const& bf : bfeatures) {
    if (bf->is_enabled()) {
      // get interior elems and then boundaries
      // this is a noop if hybrid is not enabled
      sim.add_hybrid(bf->init_hybrid(1.0), bf->get_body() );
    }
  }

  sim.set_initialized();

  // check init for blow-up or errors
  sim_err_msg = sim.check_initialization();

  if (sim_err_msg.empty()) {

    // take the 0 step (find state at t=0)
    sim.first_step();

  } else {

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
        if (mf->is_enabled()) {
          move_t newMoveType = fixed;
          if (mf->moves() or mf->emits()) {
            newMoveType = lagrangian;
          }
          sim.add_elements( mf->step_elements(rparams.tracer_scale*sim.get_ips()), inert, newMoveType, mf->get_body() );
        }
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

  // save data at final step
  if (true) {
    sim.write_vtk();
  }

  // Cleanup
  std::cout << "Starting shutdown procedure" << std::endl;
  sim.reset();
  std::cout << "Quitting" << std::endl;

  return 0;
}

