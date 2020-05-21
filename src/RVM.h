/*
 * RVM.h - the Random Vortex Method for diffusion in 2D
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Core.h"
#include "VectorHelper.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <array>
#include <vector>
#include <random>

//
// Class to hold RVM parameters and temporaries
//
// templatized on storage type ST, compute type CT, and max number of moments to solve for
//
template <class ST>
class RVM {
public:
  RVM();

  // all-to-all diffuse; can change array sizes
  void diffuse_all(std::array<Vector<ST>,2>&,
                   const Vector<ST>&,
                   const Vector<ST>&,
                   const ST);

  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

private:

};

// primary constructor
template <class ST>
RVM<ST>::RVM() {}


//
// Apply the random vortex method to the particles
//
template <class ST>
void RVM<ST>::diffuse_all(std::array<Vector<ST>,2>& pos,
                          const Vector<ST>& str,
                          const Vector<ST>& rad,
                          const ST h_nu) {

  // make sure all vector sizes are identical
  assert(pos[0].size()==pos[1].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==str.size() && "Input arrays are not uniform size");
  assert(str.size()==rad.size() && "Input arrays are not uniform size");
  const size_t n = rad.size();

  std::cout << "  Running RVM with n " << n << std::endl;

  // start timer
  auto start = std::chrono::system_clock::now();

  // init the random number generator
  std::random_device rd{};
  std::mt19937 gen{rd()};

  // create a normal distribution rng with mean 0 and std deviation h_nu
  std::normal_distribution<ST> diffuse{0.0, h_nu};

  // for each particle (can parallelize this part)
  // note that an OpenMP loop here will need to use int32_t as the counter variable type
  for (size_t i=0; i<n; ++i) {
    // apply the random walk
    pos[0][i] += diffuse(gen);
    pos[1][i] += diffuse(gen);
  }

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    rvm.diffuse_all:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}

//
// read/write parameters to json
//

// create and write a json object for all diffusion parameters
template <class ST>
void RVM<ST>::from_json(const nlohmann::json simj) {

  /*
  if (simj.find("RVM") != simj.end()) {
    nlohmann::json j = simj["RVM"];

    if (j.find("ignoreBelow") != j.end()) {
      ignore_thresh = j["ignoreBelow"];
      std::cout << "  setting ignore_thresh= " << ignore_thresh << std::endl;
    }

    if (j.find("relativeThresholds") != j.end()) {
      thresholds_are_relative = j["relativeThresholds"];
      std::cout << "  setting thresholds_are_relative= " << thresholds_are_relative << std::endl;
    }
  }
  */
}

// create and write a json object for all diffusion parameters
template <class ST>
void RVM<ST>::add_to_json(nlohmann::json& simj) const {

  /*
  // set rvm-specific parameters
  nlohmann::json j;
  j["ignoreBelow"] = ignore_thresh;
  j["relativeThresholds"] = thresholds_are_relative;
  simj["RVM"] = j;
  */
}

