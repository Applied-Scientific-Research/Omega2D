/*
 * Influence.h - Non-class influence calculations
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Boundaries.h"
#include "Kernels.h"
#include "Particles.h"
#include "Vorticity.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <array>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <vector>
#include <cmath>

//
// Particle stuff
//

//
// naive caller for the O(N^2) particle-particle kernel
//
template <class S, class A, uint8_t SRCSTRIDE, uint8_t TARGSTRIDE>
float add_influence(const std::vector<S>& srcpos,
                    const std::vector<S>& targpos,
                    std::vector<S>&       targvel) {

  //std::cout << "      inside add_influence(vec, vec, vec)" << std::endl;

  // make sure we have enough room for the results
  assert(targpos.size() == targvel.size());
  assert(SRCSTRIDE > 2);
  assert(TARGSTRIDE > 1);

  // accumulate results into targvel
  #pragma omp parallel for
  for (size_t i=0; i<targpos.size(); i+=TARGSTRIDE) {

    // generate accumulator
    std::array<A,TARGSTRIDE> accum = {{0.0}};

    // iterate and accumulate
    if (SRCSTRIDE == 4 and TARGSTRIDE == 4) {
      for (size_t j=0; j<srcpos.size(); j+=SRCSTRIDE) {
        nbody_kernel_44<S,A>(&srcpos[j], &targpos[i], accum.data());
      }
    } else if (SRCSTRIDE == 4 and TARGSTRIDE == 2) {
      for (size_t j=0; j<srcpos.size(); j+=SRCSTRIDE) {
        nbody_kernel_42<S,A>(&srcpos[j], &targpos[i], accum.data());
      }
    }

    // add to running sum
    for (size_t j=0; j<TARGSTRIDE; ++j) {
      targvel[i+j] += accum[j];
    }
    //std::cout << "  targ pt " << i/TARGSTRIDE << " has new vel " << targvel[i] << " " << targvel[i+1] << std::endl;
  }

  float flops = (float)(targpos.size()/TARGSTRIDE) * (float)(srcpos.size()/SRCSTRIDE);
  if (TARGSTRIDE == 4) { flops *= 15.0; }
  if (TARGSTRIDE == 2) { flops *= 13.0; }

  return flops;
}

#ifdef USE_VC

// convert a std::vector into native Vc array
template <class S>
inline const Vc::Memory<Vc::Vector<S>> deinterleave (const std::vector<S>& in,
                                                     const size_t idx, const size_t nper,
                                                     const float defaultval) {
  Vc::Memory<Vc::Vector<S>> out(in.size()/nper);
  for (size_t i=0; i<in.size()/nper; ++i) out[i] = in[nper*i+idx];
  for (size_t i=in.size()/nper; i<out.vectorsCount()*Vc::Vector<S>::size(); ++i) out[i] = defaultval;
  return out;
}

//
// naive caller for the O(N^2) particle-particle kernel - using Vc
//
template <class S, class A, uint8_t SRCSTRIDE, uint8_t TARGSTRIDE>
float add_influence_Vc(const std::vector<S>& srcpos,
                       const std::vector<S>& targpos,
                       std::vector<S>&       targvel) {

  //std::cout << "      inside add_influence_Vc(vec, vec, vec)" << std::endl;

  // make sure we have enough room for the results
  assert(targpos.size() == targvel.size());
  assert(SRCSTRIDE > 2);
  assert(TARGSTRIDE > 1);

  typedef Vc::Vector<S> StoreVec;
  // make a type that has just as many entries as Vector<S>
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;
  //std::cout << "      inside add_influence_Vc, source simd have " << Vc::Vector<S>::size() << " entries" << std::endl;
  //std::cout << "        and " << sizeof(Vc::Vector<S>) << " bytes per simd vector" << std::endl;
  //std::cout << "      inside add_influence_Vc, accumu simd have " << AccumVec::size() << " entries" << std::endl;
  //std::cout << "        and " << sizeof(AccumVec) << " bytes per simd vector" << std::endl;

  // create the simd-ized versions of the source points
  //std::vector<StoreVec> sx, sy, ss, sr;
  const Vc::Memory<StoreVec> sx = deinterleave(srcpos, 0, 4, 0.0);
  const Vc::Memory<StoreVec> sy = deinterleave(srcpos, 1, 4, 0.0);
  const Vc::Memory<StoreVec> ss = deinterleave(srcpos, 2, 4, 0.0);
  const Vc::Memory<StoreVec> sr = deinterleave(srcpos, 3, 4, 1.0);

  // accumulate results into targvel
  #pragma omp parallel for
  for (size_t i=0; i<targpos.size(); i+=TARGSTRIDE) {

    // spread the target out over a vector
    const StoreVec vtx = targpos[i+0];
    const StoreVec vty = targpos[i+1];

    // generate accumulator (only for vels)
    AccumVec accumu(0.0);
    AccumVec accumv(0.0);

    // iterate and accumulate
    if (SRCSTRIDE == 4 and TARGSTRIDE == 4) {
      const StoreVec vtr = targpos[i+3];
      for (size_t j=0; j<sx.vectorsCount(); ++j) {
        nbody_kernel_Vc_44<StoreVec,AccumVec>(sx.vector(j), sy.vector(j), ss.vector(j), sr.vector(j),
                                              vtx, vty, vtr, &accumu, &accumv);
      }
    } else if (SRCSTRIDE == 4 and TARGSTRIDE == 2) {
      for (size_t j=0; j<sx.vectorsCount(); ++j) {
        nbody_kernel_Vc_42<StoreVec,AccumVec>(sx.vector(j), sy.vector(j), ss.vector(j), sr.vector(j),
                                              vtx, vty, &accumu, &accumv);
      }
    }

    // add to running sum
    targvel[i+0] += accumu.sum();
    targvel[i+1] += accumv.sum();
    //std::cout << "  targ pt " << i/TARGSTRIDE << " has new vel " << targvel[i] << " " << targvel[i+1] << std::endl;
  }

  float flops = (float)(targpos.size()/TARGSTRIDE) * (float)(srcpos.size()/SRCSTRIDE);
  if (TARGSTRIDE == 4) { flops *= 15.0; }
  if (TARGSTRIDE == 2) { flops *= 13.0; }

  return flops;
}
#endif

//
// Particles affect Particles
// HACK - assumes Elements objects are actually Particles! see architecture for a resolution for this
//
template <class S, class A, class I>
float add_influence_pp (Elements<S,I>& _src, Elements<S,I>& _targ) {
//void add_influence_pp (Particles<S,I>& _src, Particles<S,I>& _targ) {
  //std::cout << "    inside add_influence(Particles, Particles)" << std::endl;
#ifdef USE_VC
  return add_influence_Vc<S,A,4,4>(_src.get_x(), _targ.get_x(), _targ.get_u());
#else
  return add_influence<S,A,4,4>(_src.get_x(), _targ.get_x(), _targ.get_u());
#endif
}

//
// High-level driver for all-affects-all
//
template <class S, class A, class I>
void add_influence (Vorticity<S,I>& _src, Vorticity<S,I>& _targ) {
  std::cout << "  inside add_influence(Vorticity, Vorticity)" << std::endl;

  // start timers
  std::chrono::system_clock::time_point start, end;
  std::chrono::duration<double> elapsed_seconds;

  for (auto & targ_elem: _targ.get_collections()) {
    for (auto & src_elem: _src.get_collections()) {
      start = std::chrono::system_clock::now();

      std::cout << "    computing influence of " << src_elem->get_n() << " particles on " << targ_elem->get_n() << " particles" << std::endl;
      // here's the problem: the routine here doesn't know what type each of these is!
      // HACK - let's assume it's always Particles
      const float flops = add_influence_pp<S,A,I>(*src_elem, *targ_elem);

      end = std::chrono::system_clock::now();
      elapsed_seconds = end-start;
      const float gflops = (flops / 1.e+9) / (float)elapsed_seconds.count();
      printf("    add_influence_pp:\t[%.6f] cpu seconds at [%.3f] GFlop/s\n", (float)elapsed_seconds.count(), gflops);
    }
  }
}


//
// Panel stuff
//

#ifdef USE_VC
//
// naive caller for the O(N^2) particle-panel kernel
//
template <class S, class A, class I>
float add_influence_ppan_Vc (Elements<S,I> const & _src, Panels<S,I>& _targ) {

  // get handles for the vectors
  std::vector<S> const& sx = _src.get_x();	// contains 0=x, 1=y, 2=str, 3=rad
  std::vector<S> const& tx = _targ.get_x();
  std::vector<I> const& ti = _targ.get_idx();
  std::vector<S>&       tu = _targ.get_vel();

  // make sure we have enough room for the results
  assert(tx.size() == tu.size());

  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

  // sources are points, so de-interleave
  const Vc::Memory<StoreVec> vsx = deinterleave(sx, 0, 4, 0.0);
  const Vc::Memory<StoreVec> vsy = deinterleave(sx, 1, 4, 0.0);
  const Vc::Memory<StoreVec> vss = deinterleave(sx, 2, 4, 0.0);
  //const Vc::Memory<StoreVec> sr = deinterleave(sx, 3, 4, 1.0);

  // accumulate results into targvel
  #pragma omp parallel for
  for (size_t i=0; i<_targ.get_n(); ++i) {

    // helper constants
    const I id1 = ti[2*i];
    const I id2 = ti[2*i+1];

    // scale by the panel size
    const AccumVec plen = 1.0 / std::sqrt(std::pow(tx[2*id2]-tx[2*id1],2) + std::pow(tx[2*id2+1]-tx[2*id1+1],2));

    // spread the target out over a vector
    const StoreVec vtx1 = tx[2*id1];
    const StoreVec vty1 = tx[2*id1+1];
    const StoreVec vtx2 = tx[2*id2];
    const StoreVec vty2 = tx[2*id2+1];

    // generate accumulator
    AccumVec accumu(0.0);
    AccumVec accumv(0.0);
    AccumVec resultu(0.0);
    AccumVec resultv(0.0);

    // iterate and accumulate
    for (size_t j=0; j<vsx.vectorsCount(); ++j) {
      vortex_panel_affects_point_Vc<StoreVec,AccumVec>(vtx1, vty1, vtx2, vty2,
                                                       vss.vector(j), vsx.vector(j), vsy.vector(j),
                                                       &resultu, &resultv);
      // and add to accumulator
      accumu += plen*resultu;
      accumv += plen*resultv;
    }

    // subtract from running sum because of way that kernel is written
    tu[2*i+0] -= accumu.sum();
    tu[2*i+1] -= accumv.sum();
    //if (std::isnan(tu[2*i+0]) or std::isnan(tu[2*i+1])) exit(0);
  }

  return (float)_targ.get_n() * (9.0 + 42.0 * (float)_src.get_n());
}
#endif

//
// naive caller for the O(N^2) particle-panel kernel
//
template <class S, class A, class I>
float add_influence_ppan (Elements<S,I> const & _src, Panels<S,I>& _targ) {

  // get handles for the vectors
  std::vector<S> const& sx = _src.get_x();	// contains 0=x, 1=y, 2=str, 3=rad
  std::vector<S> const& tx = _targ.get_x();
  std::vector<I> const& ti = _targ.get_idx();
  std::vector<S>&       tu = _targ.get_vel();

  // make sure we have enough room for the results
  assert(tx.size() == tu.size());

  // accumulate results into targvel
  #pragma omp parallel for
  for (size_t i=0; i<_targ.get_n(); ++i) {

    // helper constants
    const I id1 = ti[2*i];
    const I id2 = ti[2*i+1];

    // scale by the panel size
    const A plen = 1.0 / std::sqrt(std::pow(tx[2*id2]-tx[2*id1],2) + std::pow(tx[2*id2+1]-tx[2*id1+1],2));

    // generate accumulator
    std::array<A,2> accum = {{0.0}};
    std::array<A,2> result = {{0.0}};

    // iterate and accumulate
    for (size_t j=0; j<_src.get_n(); ++j) {
      result = vortex_panel_affects_point<S,A>(tx[2*id1], tx[2*id1+1],
                                               tx[2*id2], tx[2*id2+1],
                                               sx[4*j+2], sx[4*j+0], sx[4*j+1]);
      // and add to accumulator
      for (size_t j=0; j<2; ++j) {
        accum[j] += plen*result[j];
      }
    }

    // subtract from running sum because of way that kernel is written
    for (size_t j=0; j<2; ++j) {
      tu[2*i+j] -= accum[j];
    }
  }

  return (float)_targ.get_n() * (9.0 + 42.0 * (float)_src.get_n());
}

#ifdef USE_VC
//
// naive caller for the O(N^2) panel-particle kernel for Vc
//
template <class S, class A, class I>
float add_influence_panp_Vc (Panels<S,I> const & _src, Elements<S,I>& _targ) {

  // get handles for the vectors
  std::vector<S> const& sx = _src.get_x();
  std::vector<I> const& si = _src.get_idx();
  std::vector<S> const& ss = _src.get_strengths();
  std::vector<S> const& tx = _targ.get_x();		// contains 0=x, 1=y, 2=str, 3=rad
  std::vector<S>&       tu = _targ.get_u();		// contains derivs 0=u, 1=v, 2=dstrdt, 3=draddt

  // make sure we have enough room for the results
  assert(tx.size() == tu.size());

  // define vector types for Vc (still only S==A supported here)
  typedef Vc::Vector<S> StoreVec;
  typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

  // sources are panels, assemble de-interleaved vectors
  Vc::Memory<StoreVec> vsx1(_src.get_n());
  Vc::Memory<StoreVec> vsy1(_src.get_n());
  Vc::Memory<StoreVec> vsx2(_src.get_n());
  Vc::Memory<StoreVec> vsy2(_src.get_n());
  Vc::Memory<StoreVec> vss(_src.get_n());
  for (size_t j=0; j<_src.get_n(); ++j) {
    const I id1 = 2*si[2*j];
    const I id2 = 2*si[2*j+1];
    vsx1[j]     = sx[id1];
    vsy1[j]     = sx[id1+1];
    vsx2[j]     = sx[id2];
    vsy2[j]     = sx[id2+1];
    vss[j]      = ss[j];
    //std::cout << "  src panel " << j << " has " << vsx1[j] << " " << vsy1[j] << " and str " << vss[j] << std::endl;
  }
  for (size_t j=_src.get_n(); j<vss.vectorsCount()*StoreVec::size(); ++j) {
    vsx1[j] = 0.0;
    vsy1[j] = 0.0;
    vsx2[j] = 1.0;
    vsy2[j] = 0.0;
    vss[j]  = 0.0;
    //std::cout << "  src panel " << j << " has " << vsx1[j] << " " << vsy1[j] << " and str " << vss[j] << std::endl;
  }

  // accumulate results into targvel
  #pragma omp parallel for
  for (size_t i=0; i<_targ.get_n(); ++i) {

    // spread the target points out over a vector
    const StoreVec vtx = tx[4*i+0];
    const StoreVec vty = tx[4*i+1];

    // generate accumulator
    AccumVec accumu(0.0);
    AccumVec accumv(0.0);
    AccumVec resultu(0.0);
    AccumVec resultv(0.0);

    // iterate and accumulate
    for (size_t j=0; j<vss.vectorsCount(); ++j) {
      vortex_panel_affects_point_Vc<StoreVec,AccumVec>(vsx1.vector(j), vsy1.vector(j),
                                                       vsx2.vector(j), vsy2.vector(j),
                                                       vss.vector(j), vtx, vty,
                                                       &resultu, &resultv);
      // and add to accumulator
      accumu += resultu;
      accumv += resultv;
    }
    //std::cout << "  targ pt " << i << " has old vel " << tu[4*i] << " " << tu[4*i+1] << std::endl;

    // add to running sum because of way that kernel is written
    tu[4*i+0] += accumu.sum();
    tu[4*i+1] += accumv.sum();
    //std::cout << "              and new vel " << tu[4*i] << " " << tu[4*i+1] << std::endl;
    if (std::isnan(tu[4*i])) exit(0);
    //std::cout << "  targ pt " << i << " has new vel " << tu[4*i] << " " << tu[4*i+1] << std::endl;
  }

  return (float)_targ.get_n() * (2.0 + 40.0 * (float)_src.get_n());
}
#endif

//
// naive caller for the O(N^2) panel-particle kernel
//
template <class S, class A, class I>
float add_influence_panp (Panels<S,I> const & _src, Elements<S,I>& _targ) {

  // get handles for the vectors
  std::vector<S> const& sx = _src.get_x();
  std::vector<I> const& si = _src.get_idx();
  std::vector<S> const& ss = _src.get_strengths();
  std::vector<S> const& tx = _targ.get_x();		// contains 0=x, 1=y, 2=str, 3=rad
  std::vector<S>&       tu = _targ.get_u();		// contains derivs 0=u, 1=v, 2=dstrdt, 3=draddt

  // make sure we have enough room for the results
  assert(tx.size() == tu.size());

  // accumulate results into targvel
  #pragma omp parallel for
  for (size_t i=0; i<_targ.get_n(); ++i) {

    // generate accumulator
    std::array<A,2> accum = {{0.0}};
    std::array<A,2> result = {{0.0}};

    // iterate and accumulate
    for (size_t j=0; j<_src.get_n(); ++j) {
      result = vortex_panel_affects_point<S,A>(sx[2*si[2*j]],   sx[2*si[2*j]+1],
                                               sx[2*si[2*j+1]], sx[2*si[2*j+1]+1],
                                               ss[j], tx[4*i+0], tx[4*i+1]);
      for (size_t j=0; j<2; ++j) {
        accum[j] += result[j];
      }
    }
    //std::cout << "  targ pt " << i << " has old vel " << tu[4*i] << " " << tu[4*i+1] << std::endl;

    // add to running sum because of way that kernel is written
    for (size_t j=0; j<2; ++j) {
      tu[4*i+j] += accum[j];
    }
    //std::cout << "             and new vel " << tu[4*i] << " " << tu[4*i+1] << std::endl;
    //std::cout << "  targ pt " << i << " has new vel " << tu[4*i] << " " << tu[4*i+1] << std::endl;
  }

  return (float)_targ.get_n() * (2.0 + 40.0 * (float)_src.get_n());
}

//
// High-level driver for all-affects-all
//
template <class S, class A, class I>
void add_influence (Vorticity<S,I> const & _src, Boundaries<S,I>& _targ) {
  std::cout << "  inside add_influence(Vorticity, Boundaries)" << std::endl;

  // start timers
  std::chrono::system_clock::time_point start, end;
  std::chrono::duration<double> elapsed_seconds;

  // just one set of panels for all boundaries now
  Panels<S,I>& targ_elem = _targ.get_panels();
  // iterate over all sets of source vorticities (currently only one Particles object)
  for (auto& src_elem : _src.get_collections()) {
    std::cout << "    computing influence of " << src_elem->get_n() << " particles on " << targ_elem.get_n() << " panels" << std::endl;
    start = std::chrono::system_clock::now();

    // here's the problem: the routine here doesn't know what type each of these is!
    // HACK - let's assume it's always Particles and Panels
#ifdef USE_VC
    const float flops = add_influence_ppan_Vc<S,A,I>(*src_elem, targ_elem);
#else
    const float flops = add_influence_ppan<S,A,I>(*src_elem, targ_elem);
#endif

    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    const float gflops = (flops / 1.e+9) / (float)elapsed_seconds.count();
    printf("    add_influence_ppan:\t[%.6f] cpu seconds at [%.3f] GFlop/s\n", (float)elapsed_seconds.count(), gflops);
  }
}

template <class S, class A, class I>
void add_influence (Boundaries<S,I> const & _src, Vorticity<S,I>& _targ) {
  std::cout << "  inside add_influence(Boundaries, Vorticity)" << std::endl;

  // start timers
  std::chrono::system_clock::time_point start, end;
  std::chrono::duration<double> elapsed_seconds;

  // just one set of panels for all boundaries now
  Panels<S,I> const& src_elem = _src.get_panels();
  // iterate over all sets of target vorticities (currently only one Particles object)
  for (auto& targ_elem : _targ.get_collections()) {
    std::cout << "    computing influence of " << src_elem.get_n() << " panels on " << targ_elem->get_n() << " particles" << std::endl;
    start = std::chrono::system_clock::now();

    // here's the problem: the routine here doesn't know what type each of these is!
    // HACK - let's assume it's always Particles and Panels
#ifdef USE_VC
    const float flops = add_influence_panp_Vc<S,A,I>(src_elem, *targ_elem);
#else
    const float flops = add_influence_panp<S,A,I>(src_elem, *targ_elem);
#endif

    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    const float gflops = (flops / 1.e+9) / (float)elapsed_seconds.count();
    printf("    add_influence_panp:\t[%.6f] cpu seconds at [%.3f] GFlop/s\n", (float)elapsed_seconds.count(), gflops);
  }
}

