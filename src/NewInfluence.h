/*
 * Influence.h - Non-class influence calculations
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "NewKernels.h"
#include "Points.h"
//#include "Panels.h"

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>


template <class S, class A>
void points_affect_points (Points<S> const& src, Points<S>& targ) {
  auto start = std::chrono::system_clock::now();

  // is this where we dispatch the OpenGL compute shader?

  // get references to use locally
  const std::array<Vector<S>,Dimensions>&	sx = src.get_pos();
  const Vector<S>&				sr = src.get_rad();
  const Vector<S>&				ss = src.get_str();

  const std::array<Vector<S>,Dimensions>&	tx = targ.get_pos();
  std::array<Vector<S>,Dimensions>&		tu = targ.get_vel();

#ifdef USE_VC
  // create float_v versions of the source vectors
  const Vc::Memory<Vc::Vector<S>> sxv = stdvec_to_vcvec<S>(sx[0], 0.0);
  const Vc::Memory<Vc::Vector<S>> syv = stdvec_to_vcvec<S>(sx[1], 0.0);
  const Vc::Memory<Vc::Vector<S>> srv = stdvec_to_vcvec<S>(sr,    1.0);
  const Vc::Memory<Vc::Vector<S>> ssv = stdvec_to_vcvec<S>(ss,    0.0);
#endif

  float flops = (float)targ.getn();

  // here is where we can dispatch on solver type, grads-or-not, core function, etc.?
  if (targ.is_inert()) {
    std::cout << "    0v_0p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // targets are field points

    #pragma omp parallel for
    for (size_t i=0; i<targ.getn(); ++i) {
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
      for (size_t j=0; j<src.getn(); ++j) {
        kernel_0v_0p<S,A>(sx[0][j], sx[1][j], sr[j], ss[j], 
                          tx[0][i], tx[1][i],
                          &accumu, &accumv);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
#endif // no Vc
    }
    flops *= 2.0 + 13.0*(float)src.getn();

  } else {
    std::cout << "    0v_0v compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // targets are particles
    const Vector<S>&				tr = targ.get_rad();

    #pragma omp parallel for
    for (size_t i=0; i<targ.getn(); ++i) {
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
      for (size_t j=0; j<src.getn(); ++j) {
        kernel_0v_0v<S,A>(sx[0][j], sx[1][j], sr[j], ss[j], 
                          tx[0][i], tx[1][i], tr[i],
                          &accumu, &accumv);
      }
      tu[0][i] += accumu;
      tu[1][i] += accumv;
#endif // no Vc
    }
    flops *= 2.0 + 15.0*(float)src.getn();

  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_affect_points: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}

/*
template <class S, class A>
void panels_affect_points (Panels<S> const& src, Points<S>& targ) {
  std::cout << "    1_0 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  //const Vector<S>&                      sr = src.get_rad();
  const std::vector<uint16_t>&            si = src.get_idx();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  //const Vector<S>&                      tr = targ.get_rad();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  #pragma omp parallel for
  for (size_t i=0; i<targ.getn(); ++i) {
    //std::array<A,3> accum = {0.0};
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    for (size_t j=0; j<src.getn(); ++j) {
      const size_t jp0 = si[2*j];
      const size_t jp1 = si[2*j+1];
      //kernel_1_0v<S,A>(&sx[2*si[2*j]], &sx[2*si[2*j+1]], ss[j],
      //                &tx[2*i], accum.data());
      kernel_1_0s<S,A>(sx[0][jp0], sx[1][jp0], sx[2][jp0],
                       sx[0][jp1], sx[1][jp1], sx[2][jp1],
                       ss[0][j], ss[1][j], ss[2][j],
                       tx[0][i], tx[1][i], tx[2][i],
                       //accum.data());
                       &accumu, &accumv, &accumw);
    }
    //tu[0][i] += accum[0];
    //tu[1][i] += accum[1];
    //tu[2][i] += accum[2];
    tu[0][i] += accumu;
    tu[1][i] += accumv;
    tu[2][i] += accumw;
  }
}


template <class S, class A>
void points_affect_panels (Points<S> const& src, Panels<S>& targ) {
  std::cout << "    0_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();
  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  const std::vector<uint16_t>&            ti = targ.get_idx();
  std::array<Vector<S>,Dimensions>&       tu = targ.get_vel();

  #pragma omp parallel for
  for (size_t i=0; i<targ.getn(); ++i) {
    //std::array<A,3> accum = {0.0};
    A accumu = 0.0;
    A accumv = 0.0;
    A accumw = 0.0;
    const size_t ip0 = ti[2*i];
    const size_t ip1 = ti[2*i+1];
    for (size_t j=0; j<src.getn(); ++j) {
      // note that this is the same kernel as panels_affect_points!
      //kernel_1_0v<S,A>(&tx[2*ti[2*i]], &tx[2*ti[2*i+1]], ss[j],
      //                 &sx[2*j], accum.data());
      kernel_1_0s<S,A>(tx[0][ip0], tx[1][ip0], tx[2][ip0],
                       tx[0][ip1], tx[1][ip1], tx[2][ip1],
                       ss[0][j], ss[1][j], ss[2][j],
                       sx[0][j], sx[1][j], sx[2][j],
                       //accum.data());
                       &accumu, &accumv, &accumw);
    }
    // we use it backwards, so the resulting velocities are negative
    //tu[0][i] -= accum[0];
    //tu[1][i] -= accum[1];
    //tu[2][i] -= accum[2];
    tu[0][i] -= accumu;
    tu[1][i] -= accumv;
    tu[2][i] -= accumw;
  }
}

template <class S, class A>
void panels_affect_panels (Panels<S> const& src, Panels<S>& targ) {
  std::cout << "    1_1 compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
  // not sure how to do this - find field points of one and apply a function above?
}
*/


// helper struct for dispatching through a variant
template <class A>
struct InfluenceVisitor {
  // source collection, target collection
  void operator()(Points<float> const& src, Points<float>& targ) { points_affect_points<float,A>(src, targ); } 
  //void operator()(Panels<float> const& src, Points<float>& targ) { panels_affect_points<float,A>(src, targ); } 
  //void operator()(Points<float> const& src, Panels<float>& targ) { points_affect_panels<float,A>(src, targ); } 
  //void operator()(Panels<float> const& src, Panels<float>& targ) { panels_affect_panels<float,A>(src, targ); } 
};

