/*
 * Reflect.h - Non-class particle-panel reflecting operation
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "Points.h"
#include "Surfaces.h"

#include <cstdlib>
#include <limits>
#include <vector>
#include <cmath>

enum ClosestType { panel, node };

//
// on output of this routine, jidx is:
//   -1 is a panel is the closest
//    0 if node 0 is closest
//    1 if node 1 is closest
// but then jidx will be changed to the actual panel or node index
//
template <class S> 
struct ClosestReturn { 
  int jidx; 
  S distsq; 
  S cpx, cpy;
  ClosestType disttype; 
}; 

//
// find closest distance from a point to a line segment
// logic taken from pointElemDistance2d
//
// minimum 22 flops, maximum 26, average is probably 23
//
template <class S>
ClosestReturn<S> panel_point_distance(const S sx0, const S sy0,
                                      const S sx1, const S sy1,
                                      const S tx,  const S ty) {

  ClosestReturn<S> retval;
  retval.jidx = -1;

  // segment vector
  const S bx    = sx1-sx0;
  const S by    = sy1-sy0;
  const S blensq  = bx*bx + by*by;
  //std::cout << "point is " << tx << " " << ty << std::endl;
  //std::cout << "  panel is " << sx0 << " " << sy0 << " to " << sx1 << " " << sy1 << std::endl;

  // leg vector
  const S ax    = tx-sx0;
  const S ay    = ty-sy0;

  // t is a parametric value between 0 and 1 along the segment
  const S t     = (ax*bx + ay*by) / blensq;
  //std::cout << "  t is " << t << std::endl;

  if (t > 0.0 and t < 1.0) {
    // point is closest to the segment, not the nodes
    retval.disttype = panel;
    // note that the normal of the panel is (-by, bx)/sqrt(blensq)
    //std::cout << "  returning " << std::sqrt((std::pow(ay*bx - ax*by, 2) / blensq)) << std::endl;
    retval.distsq = std::pow(ay*bx - ax*by, 2) / blensq;
    retval.cpx = (1.0-t)*sx0 + t*sx1;
    retval.cpy = (1.0-t)*sy0 + t*sy1;
    return retval;
  }

  // check distance to each end node
  retval.disttype = node;

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(ax,2) + std::pow(ay,2);
  const S rij12 = std::pow(tx-sx1,2) + std::pow(ty-sy1,2);
  //std::cout << "  rij2 is " << rij2 << " and rij12 is " << rij12 << std::endl;

  //std::cout << "  returning " << std::sqrt(std::min(rij2, rij12)) << std::endl;
  if (rij2 < rij12) {
    retval.jidx = 0;
    retval.distsq = rij2;
    retval.cpx = sx0;
    retval.cpy = sy0;
  } else {
    retval.jidx = 1;
    retval.distsq = rij12;
    retval.cpx = sx1;
    retval.cpy = sy1;
  }

  return retval;
}


//
// naive caller for the O(N^2) panel-particle reflection kernel
//
template <class S>
void reflect_panp2 (Surfaces<S> const& _src, Points<S>& _targ) {
  //std::cout << "  inside reflect(Surfaces, Points)" << std::endl;
  std::cout << "  Reflecting" << _targ.to_string() << " from near" << _src.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();

  // pre-compute the node normals
  std::array<Vector<S>,Dimensions> sn;
  for (size_t i=0; i<Dimensions; ++i) {
    sn[i].resize(_src.get_n());
    std::fill(sn[i].begin(), sn[i].end(), 0.0);
  }
  for (size_t j=0; j<_src.get_npanels(); ++j) {
    // nodes in this panel
    const Int ip0 = si[2*j];
    const Int ip1 = si[2*j+1];
    // find the normal of this panel
    const S bx = sx[0][ip1] - sx[0][ip0];
    const S by = sx[1][ip1] - sx[1][ip0];
    const S blen = 1.0 / std::sqrt(bx*bx + by*by);
    const S normx = -by * blen;
    const S normy =  bx * blen;
    // add the normal to each of the nodes
    sn[0][ip0] += normx;
    sn[1][ip0] += normy;
    sn[0][ip1] += normx;
    sn[1][ip1] += normy;
  }

  size_t num_reflected = 0;
  const S eps = 10.0*std::numeric_limits<S>::epsilon();

  #pragma omp parallel for reduction(+:num_reflected)
  for (size_t i=0; i<_targ.get_n(); ++i) {
    S mindist = std::numeric_limits<S>::max();
    std::vector<ClosestReturn<S>> hits;

    // iterate and search for closest panel
    for (size_t j=0; j<_src.get_npanels(); ++j) {
      ClosestReturn<S> result = panel_point_distance<S>(sx[0][si[2*j]],   sx[1][si[2*j]],
                                                        sx[0][si[2*j+1]], sx[1][si[2*j+1]],
                                                        tx[0][i],         tx[1][i]);

      if (result.distsq < mindist - eps) {
        // we blew the old one away
        mindist = result.distsq;
        if (result.disttype == node) {
          result.jidx = si[2*j+result.jidx];
        } else {
          result.jidx = j;
        }
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      } else if (result.distsq < mindist + eps) {
        // we are effectively the same as the old closest
        if (result.disttype == node) {
          result.jidx = si[2*j+result.jidx];
        } else {
          result.jidx = j;
        }
        hits.push_back(result);
        //std::cout << "  THIS TIES THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;
      }
    }

    // dump out the hits
    if (false) {
      std::cout << "point " << i << " is " << tx[0][i] << " " << tx[1][i] << std::endl;
      for (auto & ahit: hits) {
        if (ahit.disttype == node) {
          std::cout << "  node " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        } else {
          std::cout << "  panel " << ahit.jidx << " is " << std::sqrt(ahit.distsq) << std::endl;
        }
      }
    }

    // if no hits, then something is wrong
    assert(hits.size() > 0 && "No nearest neighbors");

    //std::cout << "  REFLECTING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // no matter how many hits, find the mean norm and mean contact point
    S normx = 0.0;
    S normy = 0.0;
    S cpx = 0.0;
    S cpy = 0.0;

    // accumulate mean normal and the mean contact point
    for (size_t k=0; k<hits.size(); ++k) {
      //std::cout << "    cp at " << hits[k].cpx << " " << hits[k].cpy << std::endl;
      const size_t j = hits[k].jidx;
      if (hits[k].disttype == panel) {
        // hit a panel, use the norm
        const S bx = sx[0][si[2*j+1]] - sx[0][si[2*j]];
        const S by = sx[1][si[2*j+1]] - sx[1][si[2*j]];
        const S blen = 1.0 / std::sqrt(bx*bx + by*by);
        normx += -by * blen;
        normy +=  bx * blen;
      } else {
        // hit a node, use the vector to the node (this could be backwards!)
        //const S bx = tx[0][i] - hits[k].cpx;
        //const S by = tx[1][i] - hits[k].cpy;
        //const S blen = 1.0 / std::sqrt(bx*bx + by*by);
        //normx += bx * blen;
        //normy += by * blen;
        // hit a node, use the cached node norm
        normx += sn[0][hits[k].jidx];
        normy += sn[1][hits[k].jidx];
      }
      cpx += hits[k].cpx;
      cpy += hits[k].cpy;
    }

    // finish computing the mean norm and mean cp
    const S normilen = 1.0 / std::sqrt(normx*normx + normy*normy);
    normx *= normilen;
    normy *= normilen;
    cpx /= (S)hits.size();
    cpy /= (S)hits.size();

    // compare this mean norm to the vector from the contact point to the particle
    const S dotp = normx*(tx[0][i]-cpx) + normy*(tx[1][i]-cpy);
    if (dotp < 0.0) {
      // this point is under the panel - reflect it off entry 0
      // this is reasonable for most cases, except very sharp angles between adjacent panels
      const S dist = std::sqrt(hits[0].distsq);
      //std::cout << "  REFLECTING " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]);
      tx[0][i] = cpx + dist*normx;
      tx[1][i] = cpy + dist*normy;
      //std::cout << " to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) << std::endl;
      //std::cout << "    TYPE 2 TO " << tx[0][i] << " " << tx[1][i] << std::endl;
      //if (std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) < 0.49) {
      //  std::cout << "    cp is " << cpx << " " << cpy << std::endl;
      //  std::cout << "    norm is " << normx << " " << normy << std::endl;
      //  assert(false && "Die");
      //}
      num_reflected++;
    }
  }

  std::cout << "    reflected " << num_reflected << " particles" << std::endl;
  const S flops = _targ.get_n() * (62.0 + 27.0*_src.get_npanels());

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    reflect_panp2:\t[%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// reflect interior particles to exterior because VRM only works in free space
//
template <class S>
void reflect_interior(std::vector<Collection> const & _bdry,
                      std::vector<Collection>&        _vort) {

  // may need to do this multiple times to clear out concave zones!
  // this should only function when _vort is Points and _bdry is Surfaces
  for (auto &targ : _vort) {
    if (std::holds_alternative<Points<S>>(targ)) {
      Points<S>& pts = std::get<Points<S>>(targ);

      for (auto &src : _bdry) {
        if (std::holds_alternative<Surfaces<S>>(src)) {
          Surfaces<S> const & surf = std::get<Surfaces<S>>(src);

          // call the specific panels-affect-points routine
          (void) reflect_panp2<S>(surf, pts);
        }
      }
    }
  }
}


//
// generate the cut tables
//
template <class S>
std::vector<std::tuple<S,S,S>> init_cut_tables (const S _dx) {
  // get the vector ready
  std::vector<std::tuple<S,S,S>> ct;

  // all distances are normalized to default vdelta (not nominal separation)
  // going out to 3 vdeltas catches all but 1e-4 of the strength
  const S max_rad = 3.0;
  // how many entries?
  const int nx = (int)(0.5 + max_rad / _dx);
  const S dx = max_rad / (float)nx;
  //std::cout << "Making cut tables with nx " << nx << " and dx " << dx << std::endl;

  // add the first entry (remove all strength, set dshift later)
  ct.push_back(std::make_tuple((S)(-nx-0.5)*dx, 0.0, 0.0));

  // generate the entries
  S twgt = 0.0;
  S tmom = 0.0;
  for (int i=-nx; i<nx+1; ++i) {
    const S distxs = std::pow((S)i * dx, 2);

    // compute the weight of this row
    S rwgt = 0.0;
    for (int j=-nx; j<nx+1; ++j) {
      const S distxys = distxs + std::pow((S)j * dx, 2);
      // this is a pure 2D Gaussian
      rwgt += std::exp(-distxys);
      // this is for a 2D compact Gaussian
      //rwgt += std::exp(-std::pow(distxys, 1.5));
    }
    // scale by the cell size
    rwgt *= std::pow(dx, 2);
    // this is a pure 2D Gaussian
    rwgt *= 1.0/M_PI;
    // this is for a 2D compact Gaussian
    //rwgt *= 0.3526021;

    //std::cout << "  row " << i << " at " << ((S)i*dx) << " has weight " << rwgt << std::endl;

    // total weight accumulated is the fraction of strength to keep
    twgt += rwgt;

    // total moment is the amount to shift
    tmom += ((S)i*dx) * rwgt;

    // add an entry
    ct.push_back(std::make_tuple(((S)i+0.5)*dx, twgt, -tmom/twgt));
  }
  //std::cout << "  total weight " << twgt << std::endl;

  // add the last entry (keep all strength, set dshift to zero)
  ct.push_back(std::make_tuple((S)(nx+1.5)*dx, 1.0, 0.0));

  //std::cout << "Cut table is" << std::endl;
  //for (auto &entry : ct) {
    //std::cout << "  " << std::get<0>(entry) << " " << std::get<1>(entry) << " " << std::get<2>(entry) << std::endl;
  //}

  return ct;
}

//
// use the cut tables - assume _pos is normalized by vdelta
//
template <class S>
std::pair<S,S> get_cut_entry (std::vector<std::tuple<S,S,S>>& ct, const S _pos) {
  // set defaults (change nothing)
  S smult = 1.0;
  S dshift = 0.0;

  const size_t inum = ct.size();
  assert(inum > 0 && "Cut tables not initialized");

  // check vs. low and high bounds
  if (_pos < std::get<0>(ct[0])) {
    // particle is too far inside to survive
    smult = 0.0;
    dshift = 1.0 - _pos;
  } else if (_pos > std::get<0>(ct[inum-1])) {
    // particle is too far outside to be affected
    smult = 1.0;
    dshift = 0.0;
  } else {
    // linearly interpolate between two values
    for (size_t i=1; i<inum; ++i) {
      if (_pos < std::get<0>(ct[i])) {
        // point lies between this and the previous entry
        const S frac = (_pos - std::get<0>(ct[i-1])) / (std::get<0>(ct[i]) - std::get<0>(ct[i-1]));
        smult = frac*std::get<1>(ct[i]) + (1.0-frac)*std::get<1>(ct[i-1]);
        dshift = frac*std::get<2>(ct[i]) + (1.0-frac)*std::get<2>(ct[i-1]);
        //std::cout << "new dist " << _pos << " found smult " << smult << " and dshift " << dshift << std::endl;
        break;
      }
    }
  }

  return std::pair<S,S>(smult,dshift);
}


//
// naive caller for the O(N^2) panel-particle clear-inner-layer kernel
//
template <class S>
void clear_inner_panp2 (Surfaces<S> const & _src,
                        Points<S>& _targ,
                        const S _cutoff_mult,
                        const S _ips) {

  //std::cout << "  inside clear_inner_layer(Surfaces, Points)" << std::endl;
  std::cout << "  Clearing" << _targ.to_string() << " from near" << _src.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  static bool made_cut_tables = false;
  static std::vector<std::tuple<S,S,S>> ct;
  if (not made_cut_tables) {
    ct = init_cut_tables<S>(0.1);
    made_cut_tables = true;
  }

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();
  Vector<S>&                              ts = _targ.get_str();

  // if called on field points, there is no tr
  Vector<S>&                              tr = _targ.get_rad();
  const bool are_fldpts = tr.empty();

  // pre-compute the node normals
  std::array<Vector<S>,Dimensions> sn;
  for (size_t i=0; i<Dimensions; ++i) {
    sn[i].resize(_src.get_n());
    std::fill(sn[i].begin(), sn[i].end(), 0.0);
  }
  for (size_t j=0; j<_src.get_npanels(); ++j) {
    // nodes in this panel
    const Int ip0 = si[2*j];
    const Int ip1 = si[2*j+1];
    // find the normal of this panel
    const S bx = sx[0][ip1] - sx[0][ip0];
    const S by = sx[1][ip1] - sx[1][ip0];
    const S blen = 1.0 / std::sqrt(bx*bx + by*by);
    const S normx = -by * blen;
    const S normy =  bx * blen;
    // add the normal to each of the nodes
    sn[0][ip0] += normx;
    sn[1][ip0] += normy;
    sn[0][ip1] += normx;
    sn[1][ip1] += normy;
  }

  size_t num_cropped = 0;
  const S eps = 10.0*std::numeric_limits<S>::epsilon();

  #pragma omp parallel for reduction(+:num_cropped)
  for (size_t i=0; i<_targ.get_n(); ++i) {
    S mindist = std::numeric_limits<S>::max();
    std::vector<ClosestReturn<S>> hits;

    // iterate and search for closest panel/node
    for (size_t j=0; j<_src.get_npanels(); ++j) {
      ClosestReturn<S> result = panel_point_distance<S>(sx[0][si[2*j]],   sx[1][si[2*j]],
                                                        sx[0][si[2*j+1]], sx[1][si[2*j+1]],
                                                        tx[0][i],         tx[1][i]);

      if (result.distsq < mindist - eps) {
        // we blew the old one away
        mindist = result.distsq;
        if (result.disttype == node) {
          result.jidx = si[2*j+result.jidx];
        } else {
          result.jidx = j;
        }
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      } else if (result.distsq < mindist + eps) {
        // we are effectively the same as the old closest
        if (result.disttype == node) {
          result.jidx = si[2*j+result.jidx];
        } else {
          result.jidx = j;
        }
        hits.push_back(result);
        //std::cout << "  THIS TIES THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;
      }
    } // end of loop over panels

    // dump out the hits
    if (false) {
      std::cout << "point " << i << " is " << tx[0][i] << " " << tx[1][i] << std::endl;
      for (auto & ahit: hits) {
        if (ahit.disttype == node) {
          std::cout << "  node " << ahit.jidx << " dist " << std::sqrt(ahit.distsq) << std::endl;
          std::cout << "    cp is " << ahit.cpx << " " << ahit.cpy << std::endl;
        } else {
          std::cout << "  panel " << ahit.jidx << " dist " << std::sqrt(ahit.distsq) << std::endl;
          std::cout << "    cp is " << ahit.cpx << " " << ahit.cpy << std::endl;
        }
      }
    }

    // if no hits, then something is wrong
    assert(hits.size() > 0 && "No nearest neighbors");

    //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // no matter how many hits, find the mean norm and mean contact point
    S normx = 0.0;
    S normy = 0.0;
    S cpx = 0.0;
    S cpy = 0.0;

    // accumulate mean normal and the mean contact point
    for (size_t k=0; k<hits.size(); ++k) {
      //std::cout << "    cp at " << hits[k].cpx << " " << hits[k].cpy << std::endl;
      const size_t j = hits[k].jidx;
      if (hits[k].disttype == panel) {
        // hit a panel, use the norm
        const S bx = sx[0][si[2*j+1]] - sx[0][si[2*j]];
        const S by = sx[1][si[2*j+1]] - sx[1][si[2*j]];
        const S blen = 1.0 / std::sqrt(bx*bx + by*by);
        normx += -by * blen;
        normy +=  bx * blen;
        //std::cout << "    panel norm is " << (-by * blen) << " " << (bx * blen) << std::endl;
      } else {
        normx += sn[0][hits[k].jidx];
        normy += sn[1][hits[k].jidx];
        //std::cout << "    REAL norm is " << sn[0][hits[k].jidx] << " " << sn[1][hits[k].jidx] << std::endl;
      }
      cpx += hits[k].cpx;
      cpy += hits[k].cpy;
    }

    // finish computing the mean norm and mean cp
    const S normilen = 1.0 / std::sqrt(normx*normx + normy*normy);
    normx *= normilen;
    normy *= normilen;
    cpx /= (S)hits.size();
    cpy /= (S)hits.size();
    //std::cout << "  mean norm is " << normx << " " << normy << std::endl;
    //std::cout << "  mean cp is " << cpx << " " << cpy << std::endl;

    // compare this mean norm to the vector from the contact point to the particle
    const S dotp = normx*(tx[0][i]-cpx) + normy*(tx[1][i]-cpy) - _cutoff_mult*_ips;
    // now dotp is how far this point is above the cutoff layer
    // if dotp == 0.0 then the point is exactly on the cutoff layer, and it loses half of its strength
    // if dotp < -vdelta then the point loses all of its strength
    // if dotp > vdelta then the point passes unmodified

    const S this_radius = are_fldpts ? _ips : tr[i];

    if (dotp < this_radius) {
      //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << " because dotp " << dotp << " and norm " << norm[0] << " " << norm[1] << std::endl;

      // use precomputed table lookups for new position and remaining strength
      std::pair<S,S> entry = get_cut_entry(ct, dotp/this_radius);
      if (not are_fldpts) ts[i] *= std::get<0>(entry);
      tx[0][i] += std::get<1>(entry) * this_radius * normx;
      tx[1][i] += std::get<1>(entry) * this_radius * normy;

      //std::cout << "  SHIFTING dotp/tr " << (dotp/this_radius) << " str " << sfac << " and shift " << (shiftd/this_radius) << std::endl;
      //assert(shiftd > 0.0 && "Shift in clear is less than zero");
      //std::cout << "  PUSHING " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]);

      //std::cout << "    to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) << " and scale str by " << sfac << std::endl;
      // do not change radius yet
      //std::cout << "    TO " << tx[0][i] << " " << tx[1][i] << " and weaken by " << sfac << std::endl;
      //if (std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) < 0.51) {
      //  std::cout << "    dotp is " << dotp << " and tr[i] is " << this_radius << std::endl;
      //  std::cout << "    pt is " << tx[0][i] << " " << tx[1][i] << std::endl;
      //  //std::cout << "    cp is " << cpx << " " << cpy << std::endl;
      //  std::cout << "    norm is " << normx << " " << normy << std::endl;
      //  assert(false && "Die");
      //}

      num_cropped++;
    }

  } // end loop over particles

  // we did not resize the x array, so we don't need to touch the u array

  std::cout << "    cropped " << num_cropped << " particles" << std::endl;
  // flops count here is taken from reflect - might be different here
  const S flops = _targ.get_n() * (62.0 + 27.0*_src.get_npanels());

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    clear_inner_panp2:\t[%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// clean up by removing the innermost layer - the one that will be represented by boundary strengths
//
template <class S>
void clear_inner_layer(std::vector<Collection> const & _bdry,
                       std::vector<Collection>&        _vort,
                       const S                         _cutoff_factor,
                       const S                         _ips) {

  // may need to do this multiple times to clear out concave zones!
  // this should only function when _vort is Points and _bdry is Surfaces
  for (auto &targ : _vort) {
    if (std::holds_alternative<Points<S>>(targ)) {
      Points<S>& pts = std::get<Points<S>>(targ);

      for (auto &src : _bdry) {
        if (std::holds_alternative<Surfaces<S>>(src)) {
          Surfaces<S> const & surf = std::get<Surfaces<S>>(src);

          // call the specific panels-affect-points routine
          (void) clear_inner_panp2<S>(surf, pts, _cutoff_factor, _ips);
        }
      }
    }
  }
}

