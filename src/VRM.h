/*
 * VRM.h - the Vorticity Redistribution Method for 2D vortex particles
 *                 with particle size adaptivity
 *
 * (c)2017-23 Applied Scientific Research, Inc.
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

#pragma once

#include "Core.h"
#include "VectorHelper.h"
#include "GuiHelper.h"
#include "nanoflann/nanoflann.hpp"
#ifdef PLUGIN_SIMPLEX
#include "simplex.h"
#endif
#include "eigen-nnls/nnls.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

enum SolverType { nnls, simplex, lbfgsb };

//
// Class to hold VRM parameters and temporaries
//
// templatized on storage type ST, compute type CT, and max number of moments to solve for
//
template <class ST, class CT, uint8_t MAXMOM>
class VRM {
public:
  VRM();

  void set_adaptive_radii(const bool);
  void set_relative(const bool _in) { thresholds_are_relative = _in; }
  void set_ignore(const float _in) { ignore_thresh = _in; }
  void set_simplex(const bool _in) { use_solver = (_in ? simplex : nnls); }
  void set_radgrad(const float _in) { radius_lapse = _in; }
  void set_adapt(const float _in) { adapt_thresh = _in; }

  const bool get_adaptive_radii() const { return adapt_radii; }
  const bool get_relative() const { return thresholds_are_relative; }
  const float get_ignore() const { return ignore_thresh; }
  const bool get_simplex() const { return (use_solver==simplex); }
  const float get_radgrad() const { return radius_lapse; }
  const float get_adapt() const { return adapt_thresh; }

  // all-to-all diffuse; can change array sizes
  void diffuse_all(std::array<Vector<ST>,2>&,
                   Vector<ST>&,
                   Vector<ST>&,
                   const ST,
                   const CoreType,
                   const ST);

  // other functions to eventually support:
  // two-to-one merge (when particles are close to each other)
  // one-to-many elongate (re-sphericalize a stretched particle)

  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

#ifdef USE_IMGUI
  void draw_advanced();
#endif

protected:

  // class-local point list
  enum PtOrigin { globalorig, globalnew, threadlocal };
  struct TestPt {
    ST x;
    ST y;
    ST r;
    ST newr;
    ST ds;
    int32_t idx;
    PtOrigin from;
    // constructor needs everything
    TestPt(const ST _x, const ST _y, const ST _r, const ST _newr, const ST _ds, const int32_t _idx, const PtOrigin _from)
      : x(_x), y(_y), r(_r), newr(_newr), ds(_ds), idx(_idx), from(_from)
      {}
    ~TestPt() = default;
  };

private:
  // solve VRM to how many moments?
  static const int32_t num_moments = MAXMOM;
  static constexpr int32_t num_rows = (num_moments+1) * (num_moments+2) / 2;
  // we needed 16 here for static solutions, 32 for dynamic, and 64 for dynamic with adaptivity
  static constexpr int32_t max_near = 64 * num_moments;

  // new point insertion sites (normalized to h_nu and centered around origin)
  static constexpr size_t num_sites = 30 * ((MAXMOM>2) ? 2 : 1);
  std::array<ST,num_sites> xsite,ysite;
  void initialize_sites();

  // search for a new test point location
  std::pair<ST,ST> fill_neighborhood_search(const ST, const ST,
                                            const std::vector<TestPt>&,
                                            const ST);

  // set up and call the solver
  bool attempt_solution(const ST, const ST, const ST,
                        const std::vector<TestPt>&,
                        const ST,
                        const CoreType,
                        Eigen::Matrix<CT, num_rows, Eigen::Dynamic, 0, num_rows, max_near>&,
                        Eigen::Matrix<CT, num_rows, 1>&,
                        Eigen::Matrix<CT, Eigen::Dynamic, 1, 0, max_near, 1>&);

  // for standard VRM

  // do not perform VRM if source particle strength is less than
  //   this fraction of max particle strength
  ST ignore_thresh = 1.e-5;
  // are thresholds absolute or relative to strongest particle?
  bool thresholds_are_relative = true;

  // for any adaptive particle size VRM
  bool adapt_radii = false;
  ST radius_lapse = 0.15;

  // adapt size by local strength
  bool adapt_by_vort = true;
  // only adapt particles if their strength is less than this
  //   fraction of max particle strength
  ST adapt_thresh = 1.e-3;

  // adapt size by distance from the origin
  bool adapt_by_distance = false;
  ST adapt_dist_ratio_max = 0.02;
  ST adapt_dist_ratio_min = 0.01;
  //std::array<ST,2> adapt_dist_center(0.0,0.0);

  // use nanoflann for nearest-neighbor searching? false uses direct search
  const bool use_tree = true;

  SolverType use_solver = nnls;
  //SolverType use_solver = simplex;
};

// primary constructor
template <class ST, class CT, uint8_t MAXMOM>
VRM<ST,CT,MAXMOM>::VRM() {
  initialize_sites();
}

//
// initialize the particle placement sites array
//
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::initialize_sites() {

  if (MAXMOM <= 2) {
    //std::cout << "Creating one ring of " << num_sites << " insertion sites for VRM" << std::endl;

    // for first and second-moment VRM, one ring only
    if (false) {
      // original - CCW from +x
      for (size_t i=0; i<num_sites; ++i) {
        xsite[i] = std::cos(2.0*M_PI*(ST)(i)/(ST)num_sites);
        ysite[i] = std::sin(2.0*M_PI*(ST)(i)/(ST)num_sites);
      }
    } else if (false) {
      // next try - random-ish
      for (size_t i=0; i<num_sites; ++i) {
        xsite[i] = std::cos(2.0*M_PI*(ST)((17*i)%num_sites)/(ST)num_sites);
        ysite[i] = std::sin(2.0*M_PI*(ST)((17*i)%num_sites)/(ST)num_sites);
      }
    } else {
      // latest try - perfect hexagons, then CCW
      for (size_t i=0; i<num_sites/6; ++i) {
        ST theta = (ST)(i) / (ST)num_sites;
        for (size_t j=0; j<6; ++j) {
          theta += 1.0 / 6.0;
          xsite[6*i+j] = std::cos(2.0*M_PI*theta);
          ysite[6*i+j] = std::sin(2.0*M_PI*theta);
        }
      }
    }

  } else {
    std::cout << "Creating two rings of " << num_sites << " insertion sites" << std::endl;
    const size_t nring = num_sites/2;

    // to solve for 3rd-4th moments stably, we need more points
    for (size_t i=0; i<num_sites/2; ++i) {
      xsite[i] = std::cos(2.0*M_PI*(ST)((17*i)%nring)/(ST)nring);
      ysite[i] = std::sin(2.0*M_PI*(ST)((17*i)%nring)/(ST)nring);
    }

    // generate a second layer farther away
    for (size_t i=num_sites/2; i<num_sites; ++i) {
      xsite[i] = 1.85 * std::cos(2.0*M_PI*(ST)((17*i)%nring)/(ST)nring);
      ysite[i] = 1.85 * std::sin(2.0*M_PI*(ST)((17*i)%nring)/(ST)nring);
    }
  }
}

template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::set_adaptive_radii(const bool _doamr) {
  if (!adapt_radii and _doamr) std::cout << "Particle radii will adapt to solution" << std::endl;
  if (adapt_radii and !_doamr) std::cout << "Particle radii will not adapt to solution" << std::endl;
  adapt_radii = _doamr;
}

//
// use a ring of sites to determine the location of a new particle
//   searching against a unified list of test points
//
template <class ST, class CT, uint8_t MAXMOM>
std::pair<ST,ST> VRM<ST,CT,MAXMOM>::fill_neighborhood_search(const ST xc,
                                                             const ST yc,
                                                             const std::vector<TestPt>& pts,
                                                             const ST nom_sep) {

  // create array of potential sites
  std::array<ST,num_sites> tx,ty,nearest;
  for (size_t i=0; i<num_sites; ++i) {
    tx[i] = xc + nom_sep * xsite[i];
    ty[i] = yc + nom_sep * ysite[i];
  }

  // test all points vs. all sites
  size_t iopen = 0;
  ST maxmindist = 0.0;
  for (size_t i=0; i<num_sites; ++i) {
    // find the nearest particle to this site
    ST mindistsq = nom_sep * nom_sep;
    for (size_t j=0; j<pts.size(); ++j) {
      ST distsq = std::pow(pts[j].x-tx[i], 2) + std::pow(pts[j].y-ty[i], 2);
      if (distsq < mindistsq) mindistsq = distsq;
    }
    nearest[i] = mindistsq;
    if (mindistsq > maxmindist) {
      maxmindist = mindistsq;
      iopen = i;
    }
  }

  //std::cout << "    creating particle at " << tx[iopen] << " " << ty[iopen]
  //          << ", now there are " << n << " particles" << std::endl;

  return std::pair<ST,ST>(tx[iopen], ty[iopen]);
}


//
// Find the change in strength and radius that would occur over one dt and apply it
//
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::diffuse_all(std::array<Vector<ST>,2>& pos,
                                    Vector<ST>& str,
                                    Vector<ST>& rad,
                                    const ST h_nu,
                                    const CoreType core_func,
                                    const ST particle_overlap) {

  // make sure all vector sizes are identical
  assert(pos[0].size()==pos[1].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==str.size() && "Input arrays are not uniform size");
  assert(str.size()==rad.size() && "Input arrays are not uniform size");
  size_t n = rad.size();

  std::cout << "  Running VRM with n " << n << std::endl;

  // start timer
  auto start = std::chrono::steady_clock::now();

  // reference or generate the local set of vectors
  Vector<ST>& x = pos[0];
  Vector<ST>& y = pos[1];
  Vector<ST>& r = rad;
  Vector<ST>& s = str;

  Vector<ST> newr, ds;
  newr.resize(n);
  ds.resize(n);

  // the adaptivity criterion
  Vector<ST> crit;
  crit.resize(n);

  // save whether we will VRM any given particle
  std::vector<bool> tovrm;
  tovrm.resize(n);
  std::fill(tovrm.begin(), tovrm.end(), true);

  // zero out delta vector
  std::fill(ds.begin(), ds.end(), 0.0);

  // and copy the new radius
  std::copy(r.begin(), r.end(), newr.begin());

  // fill the adaptivity criterion
  ST maxCrit = 0.0;
  for (size_t i=0; i<n; ++i) {
    // just circulation magnitude
    crit[i] = std::abs(s[i]);
    // circulation times radius
    //crit[i] = r[i] * std::abs(s[i]);

    if (crit[i] > maxCrit) maxCrit = crit[i];
  }

  const int32_t minNearby = 7;
  const int32_t maxNewParts = num_moments*8 - 4;

  // the Ixx and Iyy moments of these core functions is half of the 2nd radial moment
  const CT core_second_mom = get_core_second_mom<CT>(core_func);

  // what is maximum strength of all particles?
  const ST maxStr = s[std::max_element(s.begin(), s.end()) - s.begin()];
  const ST minStr = s[std::min_element(s.begin(), s.end()) - s.begin()];
  const ST maxAbsStr = std::max(maxStr, -1.f*minStr);
  //std::cout << "maxAbsStr " << maxAbsStr << std::endl;

  // convert particle positions into something nanoflann can understand
  Eigen::Matrix<ST, Eigen::Dynamic, Dimensions> xp;
  xp.resize(n,Dimensions);
  xp.col(0) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(x.data(), n);
  xp.col(1) = Eigen::Map<Eigen::Matrix<ST, Eigen::Dynamic, 1> >(y.data(), n);
  
  // generate the searchable data structure
  typedef typename Eigen::Matrix<ST, Eigen::Dynamic, 2> EigenMatType;
  typedef typename EigenMatType::Index EigenIndexType;
  typedef nanoflann::KDTreeEigenMatrixAdaptor< EigenMatType >  my_kd_tree_t;
  my_kd_tree_t mat_index(Dimensions, std::cref(xp));
  if (use_tree) mat_index.index->buildIndex();

  nanoflann::SearchParams params;
  params.sorted = true;

  // Adaptive radius pre-search
  if (adapt_radii) {
    int32_t nshrink = 0;
    int32_t nvrm = 0;
    int32_t nspread = 0;
    int32_t nignore = 0;
    auto astart = std::chrono::steady_clock::now();

    #pragma omp parallel
    {

    std::vector<std::pair<EigenIndexType,ST> > ret_matches;
    ret_matches.reserve(max_near);

    // note that OpenMP loops need int32_t as the counter variable type
    #pragma omp for reduction(+:nshrink,nvrm,nspread,nignore)
    for (int32_t i=0; i<(int32_t)n; ++i) {
      // compute all new radii before performing any VRM

      // nominal separation for this particle
      const ST nom_sep = r[i] / particle_overlap;

      // what are h_nu (global), search radius, and new particle distance?
      const ST search_rad = nom_sep * ((num_moments > 2) ? 2.5 : 1.6);
      const ST distsq_thresh = std::pow(search_rad, 2);

      // initialize vector of indexes of nearest particles
      std::vector<int32_t> inear;

      // switch on search method
      if (use_tree) {
        // tree-based search with nanoflann
        const ST query_pt[2] = { x[i], y[i] };
        (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);

        // copy the indexes into my vector
        for (size_t j=0; j<ret_matches.size(); ++j) inear.push_back(ret_matches[j].first);

      } else {
        // direct search: look for all neighboring particles
        for (size_t j=0; j<n; ++j) {
          ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2);
          if (distsq < distsq_thresh) inear.push_back(j);
        }
      }

      // lets figure out how much larger we can make this particle

      // given our limit to the spatial gradient of particle size
      ST lapserad = 2.f * r[i];  // start with a large radius
      for (size_t j=0; j<inear.size(); ++j) {
        const int32_t jdx = inear[j];
        if (jdx != i) {
          ST dist = std::sqrt( std::pow(x[jdx]-x[i], 2) + std::pow(y[jdx]-y[i], 2) );
          // find the the largest that the current particle could become given a smaller nearby particle
          ST thisrad = r[jdx] + radius_lapse * dist;
          // and update if this is smaller than our current guess
          if (thisrad < lapserad) lapserad = thisrad;
        }
      }
      // lapserad is now an upper bound on particle size
      //std::cout << " can lapse to " << lapserad;

      // and how much should core-spreading grow this particle given no other constraints?
      //const ST growrad = std::sqrt( std::pow(r[i], 2) + 4.0*std::pow(h_nu, 2));
      const ST growrad = std::sqrt( std::pow(r[i], 2) + (2.0/core_second_mom)*std::pow(h_nu, 2));
      // growrad is another upper bound on particle size
      //std::cout << " and spread to " << growrad;

      // grow, but never grow more than the above limits dictate
      ST newrad = std::min(lapserad, growrad);

      // only go half as far as this, to allow for subtle changes in the distribution
      //   without having to shrink the particle later
      if (newrad > r[i]) {
        //newrad = 0.5*newrad + 0.5*r[i];
        newrad = 0.25*newrad + 0.75*r[i];
      } else {
        // if there are smaller nearby particles, we will need to shrink - just not too fast
        // or if they are equal, we continue with the same radius
        newrad = std::max(newrad, (ST)(0.9*r[i]));
        nshrink++;
        //std::cout << "shrink " << x[i] << " " << y[i] << std::endl;
      }

      // now use this to set the new particle size!
      bool allow_dist_to_adapt = true;
      bool allow_vort_to_adapt = true;

      // top priority: if the particle is too weak to VRM at all, grow it
      if (std::abs(s[i]) < ignore_thresh * (thresholds_are_relative ? maxAbsStr : 1.0)) {
      //if (crit[i] < ignore_thresh * (thresholds_are_relative ? maxCrit : 1.0)) {

        // some of this particle's diffusion goes into core-spreading
        newr[i] = newrad;
        // or ALL of it goes into core-spreading
        //newr[i] = growrad;

        tovrm[i] = false;
        nignore++;
        // don't allow any other adaptivity criteria
        allow_dist_to_adapt = false;
        allow_vort_to_adapt = false;
      }

      // second priority: adapt by distance from origin
      if (adapt_by_distance and allow_dist_to_adapt) {
        const ST dist = std::sqrt(std::pow(x[i], 2) + std::pow(y[i], 2));

        // particle is too small for its distance - core-spread it
        if (r[i]/dist < adapt_dist_ratio_min) {
          // some of this particle's diffusion goes into core-spreading
          newr[i] = newrad;
          // or ALL of it goes into core-spreading
          //newr[i] = growrad;
          tovrm[i] = false;
          nignore++;
          // this overrules adaptivity by strength
          allow_vort_to_adapt = false;

        // particle is in the sweet zone, allow it to grow and VRM
        } else if (r[i]/dist < adapt_dist_ratio_max) {
          // grow, but never grow more than the above limits dictate
          newr[i] = newrad;
          tovrm[i] = true;
          if (not adapt_by_vort) nspread++;
          // but if also adapting by vorticity, let it dictate the behavior
          allow_vort_to_adapt = true;

        // particle is too large for its distance - shrink it via VRM
        } else {
          // may change core size, but only to shrink a little
          newr[i] = std::min(r[i], newrad);
          tovrm[i] = true;
          if (not adapt_by_vort) nvrm++;
          // but if also adapting by vorticity, let it dictate the behavior
          allow_vort_to_adapt = true;
        }
      }

      // adapt by strength
      if (adapt_by_vort and allow_vort_to_adapt) {

        // particle is so weak that it shouldn't VRM
        //if (std::abs(s[i]) < ignore_thresh * (thresholds_are_relative ? maxAbsStr : 1.0)) {
        if (crit[i] < ignore_thresh * (thresholds_are_relative ? maxCrit : 1.0)) {

          // some of this particle's diffusion goes into core-spreading
          newr[i] = newrad;
          // or ALL of it goes into core-spreading
          //newr[i] = growrad;

          tovrm[i] = false;
          nignore++;
          //std::cout << "ignore " << x[i] << " " << y[i] << std::endl;
          //std::cout << "ignore r=" << r[i] << " dist=" << dist << std::endl;

        // particle is weak enough to begin considering
        //} else if (std::abs(s[i]) < adapt_thresh * (thresholds_are_relative ? maxAbsStr : 1.0)) {
        } else if (crit[i] < adapt_thresh * (thresholds_are_relative ? maxCrit : 1.0)) {

          //std::cout << "  radius will change to " << newrad << std::endl;
          // grow, but never grow more than the above limits dictate
          newr[i] = newrad;
          tovrm[i] = true;
          nspread++;
          //std::cout << "spread " << x[i] << " " << y[i] << std::endl;
          //std::cout << "spread r=" << r[i] << " dist=" << dist << std::endl;

        // particle is stronger than adapt_thresh, perform standard VRM
        } else {

          // may change core size, but only to shrink a little
          newr[i] = std::min(r[i], newrad);
          tovrm[i] = true;
          nvrm++;
          //std::cout << "vrm    " << x[i] << " " << y[i] << std::endl;
          //std::cout << "vrm r=" << r[i] << " dist=" << dist << std::endl;
        }
      }

    } // end loop i over particles

    } // end omp parallel

    std::cout << "    adapt: shrink/vrm/spread/ignore " << nshrink << "/" << nvrm << "/" << nspread << "/" << nignore << std::endl;

    // finish timer and report
    auto aend = std::chrono::steady_clock::now();
    std::chrono::duration<double> a_seconds = aend-astart;
    printf("    vrm.adapt_radii:\t[%.4f] seconds\n", (float)a_seconds.count());

  } else {
    // do not adapt radii -- copy current to new
    newr = r;
    for (size_t i=0; i<n; ++i) {
      // flag very weak particles to be ignored - always use strength here
      if (std::abs(s[i]) < ignore_thresh * (thresholds_are_relative ? maxAbsStr : 1.0)) {
        tovrm[i] = false;
      }
    }
  }

  // for each particle (can parallelize this part)
  const size_t initial_n = n;
  size_t nsolved = 0;
  size_t nneibs = 0;
  //size_t ntooclose = 0;
  size_t minneibs = 999999;
  size_t maxneibs = 0;

  #pragma omp parallel
  {

  // create the matrix elements, reusable, and dynamically allocated
  //Eigen::Matrix<CT, Eigen::Dynamic, 1> fractions;

  // data for the search results
  std::vector<std::pair<EigenIndexType,ST> > ret_matches;
  ret_matches.reserve(max_near);

  // the matricies that we will repeatedly work on
  Eigen::Matrix<CT, num_rows, Eigen::Dynamic, 0, num_rows, max_near> A;
  Eigen::Matrix<CT, num_rows, 1> b;
  Eigen::Matrix<CT, Eigen::Dynamic, 1, 0, max_near, 1> fractions;
  // can we instead pass these in, AND upon re-calling this function, only add the needed row?!?

  // note that OpenMP loops need int32_t as the counter variable type
  #pragma omp for reduction(+:nsolved,nneibs) reduction(min:minneibs) reduction(max:maxneibs)
  for (int32_t i=0; i<(int32_t)initial_n; ++i) {

    // find the nearest neighbor particles
    //std::cout << "\nDiffusing particle " << i << " with strength " << s[i] << std::endl;
    //std::cout << "  part " << i << " with strength " << s[i] << std::endl;

    // if current particle will not vrm, skip out (particle can still core-spread)
    if (!tovrm[i]) continue;

    nsolved++;

    // nominal separation for this particle (insertion distance)
    const ST nom_sep = r[i] / particle_overlap;

    // what is search radius?
    const ST search_rad = nom_sep * ((num_moments > 2) ? 2.5 : 1.6);

    // initialize vector of test points
    std::vector<TestPt> pts;
    pts.reserve(20);

    // switch on search method
    if (use_tree) {
      // tree-based search with nanoflann
      const ST distsq_thresh = std::pow(search_rad, 2);
      const ST query_pt[2] = { x[i], y[i] };
      (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);
      //if (ret_matches.size() > 20) std::cout << "part " << i << " at " << x[i] << " " << y[i] << " " << z[i] << " has " << ret_matches.size() << " matches" << std::endl;

      // copy the indexes into my vector
      for (size_t j=0; j<ret_matches.size(); ++j) {
        const int32_t jj = ret_matches[j].first;
        pts.push_back(TestPt(x[jj],y[jj],r[jj],newr[jj],0.0,jj,globalorig));
      }

      // now direct search over all newer particles
      #pragma omp critical (UpdateGlobalArrays)
      {
        // cannot do this while someone else is modifying the arrays
        for (size_t j=initial_n; j<n; ++j) {
          ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2);
          if (distsq < distsq_thresh) pts.push_back(TestPt(x[j],y[j],r[j],newr[j],0.0,(int32_t)j,globalnew));
        }
      }
      //std::cout << " and " << (inear.size()-nMatches) << " matches";

      // is the closest of these "too close"? Start at 0.5 if you want to catch the most egregious.
      //if (std::sqrt(ret_matches[1].second) < 0.5*nom_sep) ntooclose++;

      //std::cout << "radiusSearch(): radius " << search_rad << " found " << inear.size();
      //std::cout << std::endl;
      //for (size_t j=0; j<ret_matches.size(); ++j) std::cout << "   " << ret_matches[j].first << "\t" << ret_matches[j].second << std::endl;
    } else {
      // direct search: look for all neighboring particles, include newly-created ones
      //   ideally this would be a tree search
      const ST distsq_thresh = std::pow(search_rad, 2);
      for (size_t j=0; j<initial_n; ++j) {
        ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2);
        if (distsq < distsq_thresh) pts.push_back(TestPt(x[j],y[j],r[j],newr[j],0.0,(int32_t)j,globalorig));
      }
      #pragma omp critical (UpdateGlobalArrays)
      {
        for (size_t j=initial_n; j<n; ++j) {
          ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2);
          if (distsq < distsq_thresh) pts.push_back(TestPt(x[j],y[j],r[j],newr[j],0.0,(int32_t)j,globalnew));
        }
      }
      //std::cout << "radiusSearch(): radius " << search_rad << " found " << inear.size() << " matches";

    }
    //std::cout << "  found " << inear.size() << " particles close to particle " << i << std::endl;
    //std::cout << " :";
    //for (size_t j=0; j<inear.size(); ++j) std::cout << " " << inear[j];
    //std::cout << std::endl;

    size_t numNewParts = 0;

    // if there are less than, say, 6 neighbors, we should just add some to the local list now
    while (pts.size() < minNearby) {
      auto newpt = fill_neighborhood_search(x[i], y[i], pts, nom_sep);
      // add to the participating list as a thread-local new particle
      pts.push_back(TestPt(newpt.first,newpt.second,r[i],newr[i],0.0,-1,threadlocal));
      ++numNewParts;
      //std::cout << "  inear is";
      //for (int32_t j=0; j<inear.size(); ++j) std::cout << " " << inear[j];
      //std::cout << std::endl;
      //std::cout << "    fill neib with new part at " << newpt[0] << " " << newpt[1] << std::endl;
    }

    // now remove close parts if we have more than max_near
    while (pts.size() > max_near) {
      // look for the part closest to the diffusing particle
      int32_t jclose = 0;
      ST distnearest = std::numeric_limits<ST>::max();
      for (size_t j=0; j<pts.size(); ++j) {
        if (pts[j].idx != i) {
          const ST distsq = std::pow(x[i]-pts[j].x, 2) + std::pow(y[i]-pts[j].y, 2);
          if (distsq < distnearest) {
            distnearest = distsq;
            jclose = j;
          }
        }
      }
      // and remove it
      //std::cout << "removing pt near " << x[i] << " " << y[i] << std::endl;
      pts.erase(pts.begin()+jclose);
    }

    bool haveSolution = false;

    // assemble the underdetermined system
    while (not haveSolution and numNewParts < maxNewParts) {
      //std::cout << "  attempt solution with " << inear.size() << " close particles" << std::endl;

      // this does the heavy lifting - assemble and solve the VRM equations for the 
      //   diffusion from particle i to test points in pts
      haveSolution = attempt_solution(x[i], y[i], r[i], pts, h_nu, core_func, A, b, fractions);

      // if that didn't work, add a particle and try again
      if (not haveSolution) {
        auto newpt = fill_neighborhood_search(x[i], y[i], pts, nom_sep);

        // either make a new local test point or replace an existing entry in pts
        if (pts.size() == max_near) {
          // replace an existing point with this newly created one
          size_t ireplace = 1;	// default is 1 because diffusing particle is probably position 0
          for (size_t j=0; j<pts.size(); ++j) {
            if (pts[j].idx != i and pts[j].from != threadlocal) {
              ireplace = j;
              break;
            }
          }
          // we are removing an original particle from the test list, but not the global list
          pts[ireplace].x = newpt.first;
          pts[ireplace].y = newpt.second;
          pts[ireplace].r = r[i];
          pts[ireplace].newr = newr[i];
          pts[ireplace].ds = 0.0;
          pts[ireplace].idx = -1;
          pts[ireplace].from = threadlocal;
          //std::cout << "  replacing pt " << ireplace << " in the list of " << inear.size() << std::endl;
          // we still iterate this because we want to quit out if we just can't find a solution
          ++numNewParts;
        } else {
          // add a new one to the local list of test points
          // add it to the local list
          pts.push_back(TestPt(newpt.first,newpt.second,r[i],newr[i],0.0,-1,threadlocal));
          ++numNewParts;
        }

        //std::cout << "    no solution, added part at " << newpt[0] << " " << newpt[1] << " " << newpt[2] << std::endl;
      }
    }

    // did we eventually reach a solution?
    if (numNewParts >= maxNewParts) {
      std::cout << "Something went wrong" << std::endl;
      std::cout << "  at " << x[i] << " " << y[i] << std::endl;
      std::cout << "  with " << pts.size() << " near neibs" << std::endl;
      std::cout << "  needed numNewParts= " << numNewParts << std::endl;
      // ideally, in this situation, we would create 6 new particles around the original particle with optimal fractions,
      //   ignoring every other nearby particle - let merge take care of the higher density later
      exit(0);
    }

    nneibs += pts.size();
    if (pts.size() < minneibs) minneibs = pts.size();
    if (pts.size() > maxneibs) maxneibs = pts.size();

    // apply those fractions to the delta vector
    //std::cout << "Added strengths" << std::endl;
    // MAKE THIS THE CRITICAL SECTION, where local parts get added to global lists
    //std::cout << "\nDiffused part " << i << " over " << pts.size() << " neibs\n";
    for (size_t j=0; j<pts.size(); ++j) {
      if (pts[j].idx == i) {
        // self-influence
        #pragma omp atomic update
        ds[i] += s[i] * (fractions(j) - 1.0);
        //std::cout << "  self added " << (fractions(j) - 1.0) << "\n";
      } else if (pts[j].from == threadlocal) {
        // test point gets elevated to full particle
        #pragma omp critical (UpdateGlobalArrays)
        {
          x.push_back(pts[j].x);
          y.push_back(pts[j].y);
          r.push_back(pts[j].r);
          newr.push_back(pts[j].newr);
          s.push_back(0.0);
          ds.push_back(s[i]*fractions(j));
          ++n;
        }
        //std::cout << "  new added " << fractions(j) << "\n";
      } else {
        // test point was an existing particle
        const size_t idx = pts[j].idx;
        #pragma omp atomic update
        ds[idx] += s[i] * fractions(j);
        //std::cout << "  existing added " << fractions(j) << "\n";
      }
      //std::cout << "  " << (s[i]*fractions(j)) << " to particle " << idx << std::endl;
      //std::cout << "  now have " << n << " particles" << std::endl;
    }

  } // end loop over all current particles

  }
  // end omp parallel

  std::cout << "    neighbors: min/avg/max " << minneibs << "/" << ((ST)nneibs / (ST)nsolved) << "/" << maxneibs << std::endl;
  //std::cout << "  number of close pairs " << (ntooclose/2) << std::endl;
  std::cout << "    after VRM, n is " << n << std::endl;

  // apply the changes to the master vectors
  assert(n==s.size() and ds.size()==s.size() && "Array size mismatch in VRM");
  assert(n==r.size() and newr.size()==r.size() && "Array size mismatch in VRM");
  for (size_t i=0; i<n; ++i) {
    s[i] += ds[i];
    r[i] = newr[i];
  }

  // finish timer and report
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    vrm.diffuse_all:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}

//
// Set up and solve the VRM equations
//
template <class ST, class CT, uint8_t MAXMOM>
bool VRM<ST,CT,MAXMOM>::attempt_solution(const ST xc, const ST yc, const ST rc,
                                         const std::vector<TestPt>& pts,
                                         const ST h_nu,
                                         const CoreType core_func,
                                         Eigen::Matrix<CT, num_rows, Eigen::Dynamic, 0, num_rows, max_near>& A,
                                         Eigen::Matrix<CT, num_rows, 1>& b,
                                         Eigen::Matrix<CT, Eigen::Dynamic, 1, 0, max_near, 1>& fractions) {

  bool haveSolution = false;

  // second moment in each direction
  // one dt should generate 4 hnu^2 of second moment, or when distances
  //   are normalized by hnu, 2.0 in each direction
  // core moments are the coefficients on the moments based on core type
  static const CT second_moment = 2.0;
  static const CT fourth_moment = 12.0;

  const CT oohnu = 1.0 / h_nu;

  // the Ixx and Iyy moments of these core functions is half of the 2nd radial moment
  static const CT core_second_mom = get_core_second_mom<CT>(core_func);
  // the Ixxxx and Iyyyy moments of these core functions is 3/8th of the 4th radial moment
  static const CT core_fourth_mom = get_core_fourth_mom<CT>(core_func);

  // for non-adaptive method and floats, 1e-6 fails immediately, 1e-5 fails quickly, 3e-5 seems to work
  // for doubles, can use 1e-6, will increase accuracy for slight performance hit (see vrm3d)
  static const CT nnls_eps = 1.e-6;
  // default to 1e-6, but drop to 1e-4 for adaptive with high overlap?
  static const CT nnls_thresh = 1.e-6;
#ifdef PLUGIN_SIMPLEX
  // default to 1e-6, but drop to 1e-4 for adaptive with high overlap?
  static const CT simplex_thresh = 1.e-6;
#endif

  // reset the arrays
  //std::cout << "\nSetting up Ax=b least-squares problem" << std::endl;
  assert(pts.size() <= static_cast<size_t>(max_near) && "Too many neighbors in VRM");
  b.setZero();
  A.resize(num_rows, pts.size());
  A.setZero();
  fractions.resize(pts.size());
  fractions.setZero();

  // fill it in
  for (size_t j=0; j<pts.size(); ++j) {
    // on a subsequent call, we've already done this for j=0..size-2
    // all distances are normalized to h_nu
    CT dx = (xc-pts[j].x) * oohnu;
    CT dy = (yc-pts[j].y) * oohnu;
    A(0,j) = 1.0;
    if (num_moments > 0) {
      A(1,j) = dx;
      A(2,j) = dy;
    }
    if (num_moments > 1) {
      A(3,j) = dx*dx;
      A(4,j) = dx*dy;
      A(5,j) = dy*dy;
      if (adapt_radii) {
        const CT toadd = core_second_mom * std::pow(pts[j].newr * oohnu, 2);
        A(3,j) += toadd;
        A(5,j) += toadd;
      }
    }
    if (num_moments > 2) {
      A(6,j) = dx*dx*dx;
      A(7,j) = dx*dx*dy;
      A(8,j) = dx*dy*dy;
      A(9,j) = dy*dy*dy;
    }
    // fourth moments should be 16?
    if (num_moments > 3) {
      A(10,j) = dx*dx*dx*dx;
      A(11,j) = dx*dx*dx*dy;
      A(12,j) = dx*dx*dy*dy;
      A(13,j) = dx*dy*dy*dy;
      A(14,j) = dy*dy*dy*dy;
      if (adapt_radii) {
        const CT toadd = core_fourth_mom * std::pow(pts[j].newr * oohnu, 4);
        A(10,j) += toadd;
        A(12,j) += toadd / 3.0;
        A(14,j) += toadd;
      }
    }
  }

  b(0) = 1.f;
  if (num_moments > 1) {
    b(3) = second_moment;
    b(5) = second_moment;
    if (adapt_radii) {
      const CT toadd = core_second_mom * std::pow(rc*oohnu, 2);
      b(3) += toadd;
      b(5) += toadd;
    }
  }
  if (num_moments > 3) {
    b(10) = fourth_moment;
    b(12) = fourth_moment / 3.0;
    b(14) = fourth_moment;
    if (adapt_radii) {
      const CT toadd = core_fourth_mom * std::pow(rc*oohnu, 4);
      b(10) += toadd;
      b(12) += toadd / 3.0;
      b(14) += toadd;
    }
  }
  //std::cout << "  Here is the matrix A^T:\n" << A.transpose() << std::endl;
  //std::cout << "  Here is the right hand side b:\n\t" << b.transpose() << std::endl;
  //std::cout << "  Here is the solution vector:\n\t" << fractions.transpose() << std::endl;


  if (use_solver == nnls) {
    //std::cout << "    using NNLS solver\n" << std::endl;

    // solve with non-negative least-squares
    Eigen::NNLS<Eigen::Matrix<CT,Eigen::Dynamic,Eigen::Dynamic> > nnls_solver(A, 100, nnls_eps);

    //std::cout << "A is" << std::endl << A << std::endl;
    //std::cout << "b is" << std::endl << b.transpose() << std::endl;
    if (nnls_solver.solve(b)) {
      fractions = nnls_solver.x();
      //std::cout << "  success! required " << nnls_solver.numLS() << " LS problems" << std::endl;
      //std::cout << "  check says " << nnls_solver.check(b) << std::endl;
    } else {
      for (size_t j=0; j<pts.size(); ++j) fractions(j) = 0.f;
      //std::cout << "  fail!" << std::endl;
    }

    //std::cout << "  fractions are:\n\t" << fractions.transpose() << std::endl;

    // measure the results
    Eigen::Matrix<CT,Eigen::Dynamic,1> err = A*fractions - b;
    //std::cout << "  error is:\n" << err.transpose() << std::endl;
    //std::cout << "  error magnitude is " << std::sqrt(err.dot(err)) << std::endl;

    // was this solution successful?
    if (err.dot(err) < nnls_thresh) {
      // this is good enough!
      haveSolution = true;
    }

  } else {
    //std::cout << "    using Simplex solver\n" << std::endl;

#ifdef PLUGIN_SIMPLEX
    // solve with simplex solver

    // first, adjust the RHS
    for (size_t j=0; j<num_rows; ++j) {
      b(j) -= 0.5 * A.row(j).sum();
    }
    //std::cout << "  New right hand side b:\n\t" << b.transpose() << std::endl;

    // finally call the solver
    CT LInfNorm = 1.0;
    int retval = undr_dtrmn_solvr<CT,num_rows,max_near>(A, b, fractions, LInfNorm);
    //std::cout << "  undr_dtrmn_solvr returned " << retval << " " << LInfNorm << std::endl;

    // if we used the simplex solver, adjust the fractions here
    fractions.array() += 0.5;
    //std::cout << "  final fractions are:\n\t" << fractions.transpose() << std::endl;

    // was this solution successful?
    if (retval == 0 and LInfNorm < 0.5 + simplex_thresh) {
      // this is good enough!
      haveSolution = true;
    }
#else
    // we should never get here
    throw "Simplex solver is not available.";
#endif
  }

  // set output fractions and result
  return haveSolution;
}

//
// read/write parameters to json
//

// read a json object and retrieve all diffusion parameters
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::from_json(const nlohmann::json simj) {

  if (simj.find("VRM") != simj.end()) {
    nlohmann::json j = simj["VRM"];

    if (j.find("ignoreBelow") != j.end()) {
      ignore_thresh = j["ignoreBelow"];
      std::cout << "  setting ignore_thresh= " << ignore_thresh << std::endl;
    }

    if (j.find("relativeThresholds") != j.end()) {
      thresholds_are_relative = j["relativeThresholds"];
      std::cout << "  setting thresholds_are_relative= " << thresholds_are_relative << std::endl;
    }

    if (j.find("solver") != j.end()) {
      std::string solverstr = j["solver"];
      if (solverstr == "simplex") {
        use_solver = simplex;
      } else {
        // default is nnls
        use_solver = nnls;
      }
    } else {
      // default is nnls
      use_solver = nnls;
    }
#ifdef PLUGIN_SIMPLEX
    std::cout << "  setting VRM solver= " << (use_solver ? "simplex" : "nnls") << std::endl;
#else
    if (use_solver == simplex) {
      std::cout << "  setting VRM solver= nnls because PLUGIN_SIMPLEX is not set" << std::endl;
      use_solver = nnls;
    } else {
      std::cout << "  setting VRM solver= nnls" << std::endl;
    }
#endif
  }

  if (simj.find("AMR") != simj.end()) {
    nlohmann::json j = simj["AMR"];

    if (j.find("radiusGradient") != j.end()) {
      radius_lapse = j["radiusGradient"];
      std::cout << "  setting radius_lapse= " << radius_lapse << std::endl;
    }

    if (j.find("byVorticity") != j.end()) {
      adapt_by_vort = j["byVorticity"];
      std::cout << "  setting adapt_by_vort= " << adapt_by_vort << std::endl;
    }

    if (j.find("spreadBelow") != j.end()) {
      adapt_thresh = j["spreadBelow"];
      std::cout << "  setting adapt_thresh= " << adapt_thresh << std::endl;
    }

    if (j.find("byDistance") != j.end()) {
      adapt_by_distance = j["byDistance"];
      std::cout << "  setting adapt_by_distance= " << adapt_by_distance << std::endl;
    }

    if (j.find("distanceRatioMin") != j.end()) {
      adapt_dist_ratio_min = j["distanceRatioMin"];
      std::cout << "  setting adapt_dist_ratio_min= " << adapt_dist_ratio_min << std::endl;
    }

    if (j.find("distanceRatioMax") != j.end()) {
      adapt_dist_ratio_max = j["distanceRatioMax"];
      std::cout << "  setting adapt_dist_ratio_max= " << adapt_dist_ratio_max << std::endl;
    }
  }
}

// create and write a json object for all diffusion parameters
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::add_to_json(nlohmann::json& simj) const {

  // set vrm-specific parameters
  nlohmann::json j;
  j["ignoreBelow"] = ignore_thresh;
  j["relativeThresholds"] = thresholds_are_relative;
  j["solver"] = use_solver ? "simplex" : "nnls";
  simj["VRM"] = j;

  // set A-VRM params
  nlohmann::json aj;
  aj["radiusGradient"] = radius_lapse;
  aj["byVorticity"] = adapt_by_vort;
  aj["spreadBelow"] = adapt_thresh;
  aj["byDistance"] = adapt_by_distance;
  aj["distanceRatioMin"] = adapt_dist_ratio_min;
  aj["distanceRatioMax"] = adapt_dist_ratio_max;
  simj["AMR"] = aj;
}

#ifdef USE_IMGUI
//
// draw advanced options parts of the GUI
//
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::draw_advanced() {

  ImGui::Checkbox("Adapt by strength/vorticity", &adapt_by_vort);
  ImGui::SameLine();
  ShowHelpMarker("Adapt particle sizes based on their strengths.");

  if (adapt_by_vort) {
    ImGui::PushItemWidth(-270);
    float thisradgrad = radius_lapse;
    ImGui::SliderFloat("Radius gradient", &thisradgrad, 0.01, 0.5f, "%.2f");
    ImGui::SameLine();
    ShowHelpMarker("During adaptive diffusion, enforce a maximum spatial gradient for particle radii.");
    radius_lapse = thisradgrad;

    float log_adapt_thresh = std::log10(adapt_thresh);
    ImGui::SliderFloat("Threshold to adapt", &log_adapt_thresh, -12, 0, "%.1f");
    ImGui::SameLine();
    ShowHelpMarker("During diffusion, allow any particles with strength less than this power-of-ten threshold to grow in size.");
    adapt_thresh = std::pow(10.f,log_adapt_thresh);
    ImGui::PopItemWidth();
  }

  ImGui::Checkbox("Adapt by distance from origin", &adapt_by_distance);
  ImGui::SameLine();
  ShowHelpMarker("Adapt particle sizes based on their distance to the origin, where larger distances allow larger particles.");

  if (adapt_by_distance) {
    ImGui::PushItemWidth(-270);

    float this_adrmax = adapt_dist_ratio_max;
    ImGui::SliderFloat("Maximum radius per distance", &this_adrmax, 0.001, 0.1f, "%.3f");
    ImGui::SameLine();
    ShowHelpMarker("During adaptive diffusion, prevent particles from growing larger than this factor times the distance to the origin.");
    adapt_dist_ratio_max = this_adrmax;

    float this_adrmin = adapt_dist_ratio_min;
    ImGui::SliderFloat("Minimum radius per distance", &this_adrmin, 0.001, 0.1f, "%.3f");
    ImGui::SameLine();
    ShowHelpMarker("During adaptive diffusion, grow particles that are smaller than this factor times the distance to the origin.");
    adapt_dist_ratio_min = this_adrmin;

    ImGui::PopItemWidth();
  }
}
#endif
