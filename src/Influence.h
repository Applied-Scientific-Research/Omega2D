/*
 * Influence.h - Non-class influence calculations
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "Kernels.h"
#include "Points.h"
#include "Surfaces.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>



//
// Vc and x86 versions of Points/Particles affecting Points/Particles
//
template <class S, class A>
void points_affect_points (Points<S> const& src, Points<S>& targ) {
  auto start = std::chrono::system_clock::now();

  // is this where we dispatch the OpenGL compute shader?

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const Vector<S>&                        sr = src.get_rad();
  const Vector<S>&                        ss = src.get_str();

  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

#ifdef USE_VC
  // create float_v versions of the source vectors
  const Vc::Memory<Vc::Vector<S>> sxv = stdvec_to_vcvec<S>(sx[0], 0.0);
  const Vc::Memory<Vc::Vector<S>> syv = stdvec_to_vcvec<S>(sx[1], 0.0);
  const Vc::Memory<Vc::Vector<S>> srv = stdvec_to_vcvec<S>(sr,    1.0);
  const Vc::Memory<Vc::Vector<S>> ssv = stdvec_to_vcvec<S>(ss,    0.0);
#endif

  float flops = (float)targ.get_n();

  // here is where we can dispatch on solver type, grads-or-not, core function, etc.?
  if (targ.is_inert()) {
    std::cout << "    0v_0p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // targets are field points

    #pragma omp parallel for
    for (size_t i=0; i<targ.get_n(); ++i) {
#ifdef USE_VC
      const Vc::Vector<S> txv = tx[0][i];
      const Vc::Vector<S> tyv = tx[1][i];
      // care must be taken if S != A, because these vectors must have the same length
      Vc::Vector<A> accumu = 0.0;
      Vc::Vector<A> accumv = 0.0;
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        kernel_0v_0p<Vc::Vector<S>,Vc::Vector<A>>(
                          sxv.vector(j), syv.vector(j), srv.vector(j), ssv.vector(j),
                          txv, tyv,
                          &accumu, &accumv);
      }
      tu[0][i] += accumu.sum();
      tu[1][i] += accumv.sum();
#else  // no Vc
      A accumu = 0.0;
      A accumv = 0.0;
      for (size_t j=0; j<src.get_n(); ++j) {
        kernel_0v_0p<S,A>(sx[0][j], sx[1][j], sr[j], ss[j], 
                          tx[0][i], tx[1][i],
                          &accumu, &accumv);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
#endif // no Vc
    }
    flops *= 2.0 + 13.0*(float)src.get_n();

  } else {
    std::cout << "    0v_0v compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // targets are particles
    const Vector<S>&				tr = targ.get_rad();

    #pragma omp parallel for
    for (size_t i=0; i<targ.get_n(); ++i) {
#ifdef USE_VC
      const Vc::Vector<S> txv = tx[0][i];
      const Vc::Vector<S> tyv = tx[1][i];
      const Vc::Vector<S> trv = tr[i];
      // care must be taken if S != A, because these vectors must have the same length
      Vc::Vector<A> accumu = 0.0;
      Vc::Vector<A> accumv = 0.0;
      for (size_t j=0; j<sxv.vectorsCount(); ++j) {
        kernel_0v_0v<Vc::Vector<S>,Vc::Vector<A>>(
                          sxv.vector(j), syv.vector(j), srv.vector(j), ssv.vector(j),
                          txv, tyv, trv,
                          &accumu, &accumv);
        /*
        if (false) {
          // this is how to print
          Vc::Vector<S> temp = sxv.vector(j,0);
          std::cout << "src " << j << " has sxv " << temp << std::endl;
        }
        */
      }
      tu[0][i] += accumu.sum();
      tu[1][i] += accumv.sum();
#else  // no Vc
      A accumu = 0.0;
      A accumv = 0.0;
      for (size_t j=0; j<src.get_n(); ++j) {
        kernel_0v_0v<S,A>(sx[0][j], sx[1][j], sr[j], ss[j], 
                          tx[0][i], tx[1][i], tr[i],
                          &accumu, &accumv);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
#endif // no Vc
    }
    flops *= 2.0 + 15.0*(float)src.get_n();
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_affect_points: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// Vc and x86 versions of Panels/Surfaces affecting Points/Particles
//
template <class S, class A>
void panels_affect_points (Surfaces<S> const& src, Points<S>& targ) {
  std::cout << "    1_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::vector<Int>&                 si = src.get_idx();
  const Vector<S>&                        vs = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  // and get the source strengths, if they exist
  const bool           have_source_strengths = src.have_src_str();
  const Vector<S>&                        ss = src.get_src_str();

  float flops = (float)targ.get_n();

#ifdef USE_VC
  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

  // sources are panels, assemble de-interleaved vectors
  Vc::Memory<StoreVec> vsx0(src.get_npanels());
  Vc::Memory<StoreVec> vsy0(src.get_npanels());
  Vc::Memory<StoreVec> vsx1(src.get_npanels());
  Vc::Memory<StoreVec> vsy1(src.get_npanels());
  Vc::Memory<StoreVec> vsvs(src.get_npanels());
  for (size_t j=0; j<src.get_npanels(); ++j) {
    const size_t id0 = si[2*j];
    const size_t id1 = si[2*j+1];
    vsx0[j]     = sx[0][id0];
    vsy0[j]     = sx[1][id0];
    vsx1[j]     = sx[0][id1];
    vsy1[j]     = sx[1][id1];
    vsvs[j]     = vs[j];
    //std::cout << "  src panel " << j << " has " << vsx1[j] << " " << vsy1[j] << " and str " << vsvs[j] << std::endl;
  }
  for (size_t j=src.get_npanels(); j<vsvs.vectorsCount()*StoreVec::size(); ++j) {
    vsx0[j] = 0.0;
    vsy0[j] = 0.0;
    vsx1[j] = 1.0;
    vsy1[j] = 0.0;
    vsvs[j] = 0.0;
    //std::cout << "  src panel " << j << " has " << vsx1[j] << " " << vsy1[j] << " and str " << vsvs[j] << std::endl;
  }

  // set vectors for source source strength
  Vc::Memory<StoreVec> vsss(src.get_npanels());
  if (have_source_strengths) {
    for (size_t j=0; j<src.get_npanels(); ++j) {
      vsss[j] = ss[j];
    }
    for (size_t j=src.get_npanels(); j<vsss.vectorsCount()*StoreVec::size(); ++j) {
      vsss[j] = 0.0;
    }
  } else {
    for (size_t j=0; j<vsss.vectorsCount()*StoreVec::size(); ++j) {
      vsss[j] = 0.0;
    }
  }

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_n(); ++i) {

    // spread the target points out over a vector
    const StoreVec vtx = tx[0][i];
    const StoreVec vty = tx[1][i];

    // generate accumulator
    AccumVec accumu(0.0);
    AccumVec accumv(0.0);
    AccumVec resultu(0.0);
    AccumVec resultv(0.0);

    if (have_source_strengths) {
      for (size_t j=0; j<vsvs.vectorsCount(); ++j) {
        // note that this is the same kernel as panels_affect_points!
        kernel_1_0vs<StoreVec,AccumVec>(vsx0.vector(j), vsy0.vector(j),
                                        vsx1.vector(j), vsy1.vector(j),
                                        vsvs.vector(j), vsss.vector(j),
                                        vtx, vty,
                                        &resultu, &resultv);
        accumu += resultu;
        accumv += resultv;
      }
    } else {
      // only vortex strengths
      for (size_t j=0; j<vsvs.vectorsCount(); ++j) {
        // note that this is the same kernel as panels_affect_points!
        kernel_1_0v<StoreVec,AccumVec>(vsx0.vector(j), vsy0.vector(j),
                                       vsx1.vector(j), vsy1.vector(j),
                                       vsvs.vector(j),
                                       vtx, vty,
                                       &resultu, &resultv);
        accumu += resultu;
        accumv += resultv;
      }
    }

    // use this as normal
    tu[0][i] += accumu.sum();
    tu[1][i] += accumv.sum();
    //std::cout << "    new 1_0 vel on " << i << " is " << accumu.sum() << " " << accumv.sum() << std::endl;
  }
#else

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_n(); ++i) {

    A accumu  = 0.0;
    A accumv  = 0.0;
    A resultu = 0.0;
    A resultv = 0.0;

    if (have_source_strengths) {
      // source and vortex strengths
      for (size_t j=0; j<src.get_npanels(); ++j) {
        const size_t jp0 = si[2*j];
        const size_t jp1 = si[2*j+1];

        // note that this is the same kernel as points_affect_panels
        kernel_1_0vs<S,A>(sx[0][jp0], sx[1][jp0], 
                          sx[0][jp1], sx[1][jp1],
                          vs[j],      ss[j],
                          tx[0][i],   tx[1][i],
                          &resultu, &resultv);

        accumu += resultu;
        accumv += resultv;
      }
    } else {
      // only vortex strengths
      for (size_t j=0; j<src.get_npanels(); ++j) {
        const size_t jp0 = si[2*j];
        const size_t jp1 = si[2*j+1];

        // note that this is the same kernel as points_affect_panels
        kernel_1_0v<S,A>(sx[0][jp0], sx[1][jp0], 
                         sx[0][jp1], sx[1][jp1],
                         vs[j],
                         tx[0][i],   tx[1][i],
                         &resultu, &resultv);

        accumu += resultu;
        accumv += resultv;
      }
    }

    // use this as normal
    tu[0][i] += accumu;
    tu[1][i] += accumv;
    //std::cout << "    new 1_0 vel on " << i << " is " << accumu << " " << accumv << std::endl;
  }
#endif

  if (have_source_strengths) {
    flops *= 2.0 + 49.0*(float)src.get_npanels();
  } else {
    flops *= 2.0 + 37.0*(float)src.get_npanels();
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    panels_affect_points: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// Vc and x86 versions of Points/Particles affecting Panels/Surfaces
//
template <class S, class A>
void points_affect_panels (Points<S> const& src, Surfaces<S>& targ) {
  std::cout << "    0_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const Vector<S>&                        vs = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<Int>&                 ti = targ.get_idx();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

#ifdef USE_VC
  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

  const Vc::Memory<StoreVec> sxv = stdvec_to_vcvec<S>(sx[0], 0.0);
  const Vc::Memory<StoreVec> syv = stdvec_to_vcvec<S>(sx[1], 0.0);
  const Vc::Memory<StoreVec> vsv = stdvec_to_vcvec<S>(vs,    0.0);
#endif

  float flops = (float)targ.get_npanels();

  #pragma omp parallel for
  for (size_t i=0; i<targ.get_npanels(); ++i) {

    const size_t ip0 = ti[2*i];
    const size_t ip1 = ti[2*i+1];

    // scale by the panel size
    const A plen = 1.0 / std::sqrt(std::pow(tx[0][ip1]-tx[0][ip0],2) + std::pow(tx[1][ip1]-tx[1][ip0],2));

#ifdef USE_VC
    // spread the target out over a vector
    const StoreVec vtx0 = tx[0][ip0];
    const StoreVec vty0 = tx[1][ip0];
    const StoreVec vtx1 = tx[0][ip1];
    const StoreVec vty1 = tx[1][ip1];

    // generate accumulator
    AccumVec accumu(0.0);
    AccumVec accumv(0.0);
    AccumVec resultu(0.0);
    AccumVec resultv(0.0);

    for (size_t j=0; j<vsv.vectorsCount(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      kernel_1_0v<StoreVec,AccumVec>(vtx0, vty0,
                                     vtx1, vty1,
                                     vsv.vector(j),
                                     sxv.vector(j), syv.vector(j),
                                     &resultu, &resultv);
      accumu += resultu;
      accumv += resultv;
    }

    // but we use it backwards, so the resulting velocities are negative
    tu[0][i] -= plen*accumu.sum();
    tu[1][i] -= plen*accumv.sum();
    //std::cout << "    new 0_1 vel on " << i << " is " << (-plen*accumu.sum()) << " " << (-plen*accumv.sum()) << std::endl;

#else
    A accumu  = 0.0;
    A accumv  = 0.0;
    A resultu = 0.0;
    A resultv = 0.0;

    for (size_t j=0; j<src.get_n(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      kernel_1_0v<S,A>(tx[0][ip0], tx[1][ip0],
                       tx[0][ip1], tx[1][ip1],
                       vs[j],
                       sx[0][j],   sx[1][j],
                       &resultu, &resultv);
      accumu += resultu;
      accumv += resultv;
    }

    // but we use it backwards, so the resulting velocities are negative
    tu[0][i] -= plen*accumu;
    tu[1][i] -= plen*accumv;
    //std::cout << "  vel on " << i << " is " << tu[0][i] << " " << tu[1][i] << std::endl;
    //std::cout << "    new 0_1 vel on " << i << " is " << (-plen*accumu) << " " << (-plen*accumv) << std::endl;
#endif
  }
  flops *= 11.0 + 37.0*(float)src.get_n();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_affect_panels: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


template <class S, class A>
void panels_affect_panels (Surfaces<S> const& src, Surfaces<S>& targ) {
  std::cout << "    1_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // run panels_affect_points instead

  // generate temporary colocation points as Points - is this inefficient?
  std::vector<S> xysr = targ.represent_as_particles(0.0001, 0.0001);
  Points<float> temppts(xysr, active, lagrangian, nullptr);

  // run the calculation
  panels_affect_points<S,A>(src, temppts);

  // and copy the velocities to the real target
  std::array<Vector<S>,Dimensions>& fromvel = temppts.get_vel();
  std::array<Vector<S>,Dimensions>& tovel   = targ.get_vel();
  for (size_t i=0; i<Dimensions; ++i) {
    std::copy(fromvel[i].begin(), fromvel[i].end(), tovel[i].begin());
  }

  // get references to use locally
/*
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::vector<Int>&                 si = src.get_idx();
  const Vector<S>&                        ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<Int>&                 ti = targ.get_idx();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();
*/

  // run panels_affect_points on those
  // or do it on the fly here
  // return results
}


//
// helper struct for dispatching through a variant
//
template <class A>
struct InfluenceVisitor {
  // source collection, target collection
  void operator()(Points<float> const& src,   Points<float>& targ)   { points_affect_points<float,A>(src, targ); } 
  void operator()(Surfaces<float> const& src, Points<float>& targ)   { panels_affect_points<float,A>(src, targ); } 
  void operator()(Points<float> const& src,   Surfaces<float>& targ) { points_affect_panels<float,A>(src, targ); } 
  void operator()(Surfaces<float> const& src, Surfaces<float>& targ) { panels_affect_panels<float,A>(src, targ); } 
};

