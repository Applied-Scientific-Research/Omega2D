/*
 * VRM.h - the Vorticity Redistribution Method for 2D vortex particles
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Core.h"
#include "VectorHelper.h"
#include "nanoflann.hpp"
#ifdef PLUGIN_SIMPLEX
#include "simplex.h"
#endif
#include "nnls.h"

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

enum SolverType { nnls, simplex };

//
// Class to hold VRM parameters and temporaries
//
// templatized on storage type ST, compute type CT, and max number of moments to solve for
//
template <class ST, class CT, uint8_t MAXMOM>
class VRM {
public:
  VRM();
  VRM(const CT);

  void set_hnu(const ST);
  void set_adaptive_radii(const bool);
  void set_relative(const bool _in) { thresholds_are_relative = _in; }
  void set_ignore(const float _in) { ignore_thresh = _in; }
  void set_simplex(const bool _in) { use_solver = (_in ? simplex : nnls); }

  ST get_hnu();
  const bool get_adaptive_radii() const { return adapt_radii; }
  const bool get_relative() const { return thresholds_are_relative; }
  const float get_ignore() const { return ignore_thresh; }
  const bool get_simplex() const { return (use_solver==simplex); }

  // all-to-all diffuse; can change array sizes
  void diffuse_all(std::array<Vector<ST>,2>&,
                   Vector<ST>&,
                   Vector<ST>&,
                   const CoreType,
                   const ST);

  // other functions to eventually support:
  // two-to-one merge (when particles are close to each other)
  // one-to-many elongate (re-sphericalize a stretched particle)

  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

protected:
  // search for new target location
  std::pair<ST,ST> fill_neighborhood_search(const int32_t,
                                            const Vector<ST>&,
                                            const Vector<ST>&,
                                            const std::vector<int32_t>&,
                                            const ST);

  // set up and call the solver
  bool attempt_solution(const int32_t,
                        std::vector<int32_t>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        Vector<ST>&,
                        const CoreType,
                        Eigen::Matrix<CT, Eigen::Dynamic, 1>&);

private:
  // solve VRM to how many moments?
  static const int32_t num_moments = MAXMOM;
  static constexpr int32_t num_rows = (num_moments+1) * (num_moments+2) / 2;
  // we needed 16 here for static solutions, 32 for dynamic, and 64 for dynamic with adaptivity
  static constexpr int32_t max_near = 32 * num_moments;

  // h_nu is sqrt(dt*nu) or sqrt(dt/Re)
  CT h_nu;

  // new point insertion sites (normalized to h_nu and centered around origin)
  static constexpr size_t num_sites = 30 * ((MAXMOM>2) ? 2 : 1);
  std::array<ST,num_sites> xsite,ysite;
  void initialize_sites();

  // for adaptive particle size VRM
  bool adapt_radii = false;

  // do not perform VRM if source particle strength is less than
  //   this fraction of max particle strength
  ST ignore_thresh = 1.e-5;
  // are thresholds absolute or relative to strongest particle?
  bool thresholds_are_relative = true;

  // use nanoflann for nearest-neighbor searching? false uses direct search
  const bool use_tree = true;

  SolverType use_solver = nnls;
  //SolverType use_solver = simplex;
};

// delegating ctor
template <class ST, class CT, uint8_t MAXMOM>
VRM<ST,CT,MAXMOM>::VRM()
  : VRM(1.0)
  {}

// primary constructor
template <class ST, class CT, uint8_t MAXMOM>
VRM<ST,CT,MAXMOM>::VRM(const CT _hnu)
  : h_nu(_hnu) {

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
void VRM<ST,CT,MAXMOM>::set_hnu(const ST _newhnu) {
  h_nu = _newhnu;
}

template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::set_adaptive_radii(const bool _doamr) {
  //if (!adapt_radii and _doamr) std::cout << "Particle radii will adapt to solution" << std::endl;
  //if (adapt_radii and !_doamr) std::cout << "Particle radii will not adapt to solution" << std::endl;
  adapt_radii = _doamr;
}

template <class ST, class CT, uint8_t MAXMOM>
ST VRM<ST,CT,MAXMOM>::get_hnu() {
  return (ST)h_nu;
}

//
// use a ring of sites to determine the location of a new particle
//
template <class ST, class CT, uint8_t MAXMOM>
std::pair<ST,ST> VRM<ST,CT,MAXMOM>::fill_neighborhood_search(const int32_t idx,
                                                             const Vector<ST>& x,
                                                             const Vector<ST>& y,
                                                             const std::vector<int32_t>& inear,
                                                             const ST nom_sep) {

  // create array of potential sites
  std::array<ST,num_sites> tx,ty,nearest;
  for (size_t i=0; i<num_sites; ++i) {
    tx[i] = x[idx] + nom_sep * xsite[i];
    ty[i] = y[idx] + nom_sep * ysite[i];
  }

  // test all points vs. all sites
  size_t iopen = 0;
  ST maxmindist = 0.0;
  for (size_t i=0; i<num_sites; ++i) {
    // find the nearest particle to this site
    ST mindistsq = nom_sep * nom_sep;
    for (size_t j=0; j<inear.size(); ++j) {
      ST distsq = std::pow(x[inear[j]]-tx[i], 2) + std::pow(y[inear[j]]-ty[i], 2);
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
                                    const CoreType core_func,
                                    const ST particle_overlap) {

  // make sure all vector sizes are identical
  assert(pos[0].size()==pos[1].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==str.size() && "Input arrays are not uniform size");
  assert(str.size()==rad.size() && "Input arrays are not uniform size");
  size_t n = rad.size();

  std::cout << "  Running VRM with n " << n << std::endl;

  // start timer
  auto start = std::chrono::system_clock::now();

  // reference or generate the local set of vectors
  Vector<ST>& x = pos[0];
  Vector<ST>& y = pos[1];
  Vector<ST>& r = rad;
  Vector<ST>& s = str;

  Vector<ST> newr, ds;
  newr.resize(n);
  ds.resize(n);

  // zero out delta vector
  std::fill(ds.begin(), ds.end(), 0.0);

  // and copy the new radius
  std::copy(r.begin(), r.end(), newr.begin());

  // create the matrix elements, reusable, and dynamically allocated
  Eigen::Matrix<CT, Eigen::Dynamic, 1> fractions;

  const int32_t minNearby = 7;
  const int32_t maxNewParts = num_moments*8 - 4;

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
  my_kd_tree_t mat_index(xp, 20);
  if (use_tree) mat_index.index->buildIndex();
  std::vector<std::pair<EigenIndexType,ST> > ret_matches;
  ret_matches.reserve(max_near);
  nanoflann::SearchParams params;
  params.sorted = true;

  // do not adapt particle radii -- copy current to new
  newr = r;

  // for each particle (can parallelize this part)
  const size_t initial_n = n;
  size_t nsolved = 0;
  size_t nneibs = 0;
  //size_t ntooclose = 0;
  size_t minneibs = 999999;
  size_t maxneibs = 0;
  // note that an OpenMP loop here will need to use int32_t as the counter variable type
  for (size_t i=0; i<initial_n; ++i) {

    // find the nearest neighbor particles
    //std::cout << "\nDiffusing particle " << i << " with strength " << s[i] << std::endl;
    //std::cout << "  part " << i << " with strength " << s[i] << std::endl;

    // if current particle strength is very small, skip out
    //   (this particle could still core-spread if adaptive particle size is on)
    if ((thresholds_are_relative && (std::abs(s[i]) < maxAbsStr * ignore_thresh)) or
        (!thresholds_are_relative && (std::abs(s[i]) < ignore_thresh))) continue;

    nsolved++;

    // nominal separation for this particle (insertion distance)
    const ST nom_sep = r[i] / particle_overlap;

    // what is search radius?
    const ST search_rad = nom_sep * ((num_moments > 2) ? 2.5 : 1.6);

    // initialize vector of indexes of nearest particles
    std::vector<int32_t> inear;

    // switch on search method
    if (use_tree) {
      // tree-based search with nanoflann
      const ST distsq_thresh = std::pow(search_rad, 2);
      const ST query_pt[2] = { x[i], y[i] };
      (void) mat_index.index->radiusSearch(query_pt, distsq_thresh, ret_matches, params);
      //if (ret_matches.size() > 20) std::cout << "part " << i << " at " << x[i] << " " << y[i] << " " << z[i] << " has " << ret_matches.size() << " matches" << std::endl;

      // copy the indexes into my vector
      for (size_t j=0; j<ret_matches.size(); ++j) inear.push_back((int32_t)ret_matches[j].first);

      // now direct search over all newer particles
      for (size_t j=initial_n; j<n; ++j) {
        ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2);
        if (distsq < distsq_thresh) inear.push_back((int32_t)j);
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
      for (size_t j=0; j<n; ++j) {
        ST distsq = std::pow(x[i]-x[j], 2) + std::pow(y[i]-y[j], 2);
        if (distsq < distsq_thresh) inear.push_back((int32_t)j);
      }
      //std::cout << "radiusSearch(): radius " << search_rad << " found " << inear.size() << " matches";

    }
    //std::cout << "  found " << inear.size() << " particles close to particle " << i << std::endl;
    //std::cout << " :";
    //for (size_t j=0; j<inear.size(); ++j) std::cout << " " << inear[j];
    //std::cout << std::endl;

    // if there are less than, say, 6, we should just add some now
    while (inear.size() < minNearby) {
      auto newpt = fill_neighborhood_search(i, x, y, inear, nom_sep);
      // add the index to the near list
      inear.push_back(n);
      // add it to the master list
      x.push_back(newpt.first);
      y.push_back(newpt.second);
      const ST thisnewr = newr[i];
      r.push_back(thisnewr);
      newr.push_back(thisnewr);
      s.push_back(0.0);
      ds.push_back(0.0);
      n++;
      //std::cout << "  inear is";
      //for (int32_t j=0; j<inear.size(); ++j) std::cout << " " << inear[j];
      //std::cout << std::endl;
      //std::cout << "    fill neib with new part at " << newpt[0] << " " << newpt[1] << std::endl;
    }

    bool haveSolution = false;
    size_t numNewParts = 0;

    // assemble the underdetermined system
    while (not haveSolution and ++numNewParts < maxNewParts) {
      //std::cout << "  attempt solution with " << inear.size() << " close particles" << std::endl;

      // this does the heavy lifting - assemble and solve the VRM equations for the 
      //   diffusion from particle i to particles in inear
      haveSolution = attempt_solution(i, inear, x, y, r, newr, core_func, fractions);

      // if that didn't work, add a particle and try again
      if (not haveSolution) {
        // solution is bad, add a particle and try again
        auto newpt = fill_neighborhood_search(i, x, y, inear, nom_sep);
        // add the index to the near list
        inear.push_back(n);
        // add it to the master list
        x.push_back(newpt.first);
        y.push_back(newpt.second);
        const ST thisnewr = newr[i];
        r.push_back(thisnewr);
        newr.push_back(thisnewr);
        s.push_back(0.0);
        ds.push_back(0.0);
        n++;
        //std::cout << "    no solution, added part at " << newpt[0] << " " << newpt[1] << " " << newpt[2] << std::endl;
      }
    }

    // did we eventually reach a solution?
    if (numNewParts >= maxNewParts) {
      std::cout << "Something went wrong" << std::endl;
      std::cout << "  needed numNewParts= " << numNewParts << std::endl;
      exit(0);
    }

    nneibs += inear.size();
    if (inear.size() < minneibs) minneibs = inear.size();
    if (inear.size() > maxneibs) maxneibs = inear.size();

    // apply those fractions to the delta vector
    //std::cout << "Added strengths" << std::endl;
    for (size_t j=0; j<inear.size(); ++j) {
      size_t idx = inear[j];
      if (idx == i) {
        // self-influence
        ds[idx] += s[i] * (fractions(j) - 1.0);
      } else {
        ds[idx] += s[i] * fractions(j);
      }
      //std::cout << "  " << (s[i]*fractions(j)) << " to particle " << idx << std::endl;
    }

  } // end loop over all current particles

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
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    vrm.diffuse_all:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}

//
// Set up and solve the VRM equations
//
template <class ST, class CT, uint8_t MAXMOM>
bool VRM<ST,CT,MAXMOM>::attempt_solution(const int32_t idiff,
                                         std::vector<int32_t>& inear,
                                         Vector<ST>& x,
                                         Vector<ST>& y,
                                         Vector<ST>& r,
                                         Vector<ST>& newr,
                                         const CoreType core_func,
                                         Eigen::Matrix<CT, Eigen::Dynamic, 1>& fracout) {

  bool haveSolution = false;

  // the matricies that we will repeatedly work on
  static Eigen::Matrix<CT, num_rows, Eigen::Dynamic, 0, num_rows, max_near> A;
  static Eigen::Matrix<CT, num_rows, 1> b;
  static Eigen::Matrix<CT, Eigen::Dynamic, 1, 0, max_near, 1> fractions;

  // second moment in each direction
  // one dt should generate 4 hnu^2 of second moment, or when distances
  //   are normalized by hnu, 2.0 in each direction
  // core moments are the coefficients on the moments based on core type
  static const CT second_moment = 2.0;
  static const CT fourth_moment = 12.0;

  // the Ixx and Iyy moments of these core functions is half of the 2nd radial moment
  //static const CT core_second_mom = get_core_second_mom<CT>(core_func);
  // the Ixxxx and Iyyyy moments of these core functions is 3/8th of the 4th radial moment
  //static const CT core_fourth_mom = get_core_fourth_mom<CT>(core_func);

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
  b.setZero();
  A.resize(num_rows, inear.size());
  A.setZero();
  fractions.resize(inear.size());
  fractions.setZero();

  // fill it in
  for (size_t j=0; j<inear.size(); ++j) {
    const int32_t jdx = inear[j];
    // all distances are normalized to h_nu
    CT dx = (x[idiff]-x[jdx]) / h_nu;
    CT dy = (y[idiff]-y[jdx]) / h_nu;
    A(0,j) = 1.0;
    if (num_moments > 0) {
      A(1,j) = dx;
      A(2,j) = dy;
    }
    if (num_moments > 1) {
      A(3,j) = dx*dx;
      A(4,j) = dx*dy;
      A(5,j) = dy*dy;
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
    }
  }

  b(0) = 1.f;
  if (num_moments > 1) {
    b(3) = second_moment;
    b(5) = second_moment;
  }
  if (num_moments > 3) {
    b(10) = fourth_moment;
    b(12) = fourth_moment / 3.0;
    b(14) = fourth_moment;
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
      for (size_t j=0; j<inear.size(); ++j) fractions(j) = 0.f;
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
  fracout = fractions;
  return haveSolution;
}

//
// read/write parameters to json
//

// create and write a json object for all diffusion parameters
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
  }
}

// create and write a json object for all diffusion parameters
template <class ST, class CT, uint8_t MAXMOM>
void VRM<ST,CT,MAXMOM>::add_to_json(nlohmann::json& simj) const {

  // set vrm-specific parameters
  nlohmann::json j;
  j["ignoreBelow"] = ignore_thresh;
  j["relativeThresholds"] = thresholds_are_relative;
  simj["VRM"] = j;
}

