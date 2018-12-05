/*
 * Merge.h - library code for a two-dimensional particle merging scheme
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "nanoflann.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>


//
// Find close particles and merge them, maintaining 0,1 moments and approximating
// new second moment
//
// templated on storage class ST (typically float or double)
//
template <class ST>
size_t merge_close_particles(std::vector<ST>& xin,
                             std::vector<ST>& uin,
                             const ST         particle_overlap,
                             const ST         threshold,
                             const bool       adapt_radii) {

  // start timer
  auto start = std::chrono::system_clock::now();

  // make sure all vector sizes are identical
  assert(xin.size()==uin.size());
  const size_t n = xin.size() / 4;

  std::cout << "  Merging close particles with n " << n << std::endl;

  // generate the local set of vectors
  std::vector<ST> x, y, r, s;
  x.resize(n);
  y.resize(n);
  r.resize(n);
  s.resize(n);

  // decompose the input vectors into compatible vectors
  for (size_t i=0; i<n; ++i) {
    size_t idx = 4*i;
    x[i] = xin[idx];
    y[i] = xin[idx+1];
    s[i] = xin[idx+2];
    r[i] = xin[idx+3];
  }

  // convert particle positions into something nanoflann can understand
  Eigen::Matrix<ST, Eigen::Dynamic, 2> xp;
  xp.resize(n,2);
  xp.col(0) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(x.data(), n);
  xp.col(1) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(y.data(), n);
  
  // generate the searchable data structure
  typedef nanoflann::KDTreeEigenMatrixAdaptor< Eigen::Matrix<ST, Eigen::Dynamic, 2> >  my_kd_tree_t;
  my_kd_tree_t mat_index(xp, 20 /* max leaf */ );
  mat_index.index->buildIndex();
  std::vector<std::pair<long int,ST> > ret_matches;
  ret_matches.reserve(16);
  nanoflann::SearchParams params;
  params.sorted = true;

  //
  // merge co-located particles with identical radii
  //
  
  // prepare a vector of particles to remove
  std::vector<bool> erase_me;
  erase_me.resize(n);
  std::fill(erase_me.begin(), erase_me.end(), false);

  // now, for every particle, search for a co-located and identical-radius particle!
  for (size_t i=0; i<n; ++i) {

    // nominal separation for this particle
    const ST nom_sep = r[i] / particle_overlap;
    const ST search_rad = nom_sep * threshold;
    const ST distsq_thresh = std::pow(search_rad, 2);

    // tree-based search with nanoflann
    const ST query_pt[2] = { x[i], y[i] };
    const size_t nMatches = mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);

    // one match should be self, but we don't know which one
    // if there are more than one, check the radii
    if (nMatches > 1) {
      for (size_t j=0; j<ret_matches.size(); ++j) {
        const size_t iother = ret_matches[j].first;
        if (i != iother) {
          // make sure distance is also less than target particle's threshold
          // note that distance returned from radiusSearch is already squared
          if (std::sqrt(ret_matches[j].second) < threshold*r[iother]/particle_overlap) {
            if (erase_me[i] or erase_me[iother]) {
              // we've already account for this one
            } else {
              //std::cout << "  particles " << i << " and " << iother << " will merge" << std::endl;
              //std::cout << "    first at " << x[i] << " " << y[i] << " with str " << s[i] << " and rad " << r[i] << std::endl;
              //std::cout << "    other at " << x[iother] << " " << y[iother] << " with str " << s[iother] << " and rad " << r[iother] << std::endl;
              const ST si = std::abs(s[i]) + std::numeric_limits<ST>::epsilon();
              const ST so = std::abs(s[iother]) + std::numeric_limits<ST>::epsilon();
              const ST strength_mag = si + so;
              // find center of strength
              const ST newx = (x[i]*si + x[iother]*so) / strength_mag;
              const ST newy = (y[i]*si + y[iother]*so) / strength_mag;
              // move strengths to particle i
              x[i] = newx;
              y[i] = newy;
              if (adapt_radii) {
                r[i] = std::sqrt((si*r[i]*r[i] + so*r[iother]*r[iother])/strength_mag);
              }
              s[i] = s[i] + s[iother];
              //std::cout << "    result   " << x[i] << " " << y[i] << " with str " << s[i] << " and rad " << r[i] << std::endl;
              // flag other particle for deletion
              erase_me[iother] = true;
            }
          }
        }
      }
    }
  }

  // how many to be erased?
  const size_t num_removed = std::count (erase_me.begin(), erase_me.end(), true);

  if (num_removed > 0) {

    // now march through the arrays and compress them
    size_t copyto = 0;
    for (size_t i=0; i<n; ++i) {
      if (not erase_me[i]) {
        if (i != copyto) {
          //x[copyto] = x[i];
          //y[copyto] = y[i];
          //r[copyto] = r[i];
          //s[copyto] = s[i];
          size_t ito = 4*copyto;
          xin[ito+0] = x[i];
          xin[ito+1] = y[i];
          xin[ito+2] = s[i];
          xin[ito+3] = r[i];
          size_t ifr = 4*i;
          uin[ito+0] = uin[ifr+0];
          uin[ito+1] = uin[ifr+1];
          uin[ito+2] = uin[ifr+2];
          uin[ito+3] = uin[ifr+3];
        }
        copyto++;
      }
    }
    // reset n and resize the arrays
    const size_t new_n = copyto;
    xin.resize(4*new_n);
    uin.resize(4*new_n);
    x.clear();
    y.clear();
    r.clear();
    s.clear();

    std::cout << "    merge removed " << num_removed << " particles" << std::endl;
  }

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    merging time:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());

  return num_removed;
}

