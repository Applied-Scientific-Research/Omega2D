/*
 * Reflect.h - Non-class particle-panel reflecting operation
 *
 * (c)2017-21 Applied Scientific Research, Inc.
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

#include "Omega2D.h"
#include "Points.h"
#include "Surfaces.h"

#include <cstdlib>
#include <limits>
#include <vector>
#include <cmath>
#include <chrono>

enum ClosestType { panel, node };

//
// on output of this routine, jidx is:
//   -1 if a panel is the closest
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
/*
template <class S>
S panel_point_distance_simd(const S sx0, const S sy0,
                            const S sx1, const S sy1,
                            const S tx,  const S ty) {

  S retval;

  // segment vector
  const S bx    = sx1-sx0;
  const S by    = sy1-sy0;
  const S ooblensq  = S(1.0) / (bx*bx + by*by);
  //std::cout << "point is " << tx << " " << ty << std::endl;
  //std::cout << "  panel is " << sx0 << " " << sy0 << " to " << sx1 << " " << sy1 << std::endl;

  // leg vector
  const S ax    = tx-sx0;
  const S ay    = ty-sy0;

  // t is a parametric value between 0 and 1 along the segment
  const S t     = (ax*bx + ay*by) * ooblensq;
  //std::cout << "  t is " << t << std::endl;

  if (t > 0.0 and t < 1.0) {
    // point is closest to the segment, not the nodes
    // note that the normal of the panel is (-by, bx)/sqrt(blensq)
    //std::cout << "  returning " << std::sqrt((std::pow(ay*bx - ax*by, 2) / blensq)) << std::endl;
    retval.distsq = std::pow(ay*bx - ax*by, 2) * ooblensq;
    return retval;
  }

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(ax,2) + std::pow(ay,2);
  const S rij12 = std::pow(tx-sx1,2) + std::pow(ty-sy1,2);
  //std::cout << "  rij2 is " << rij2 << " and rij12 is " << rij12 << std::endl;

  //std::cout << "  returning " << std::sqrt(std::min(rij2, rij12)) << std::endl;
  if (rij2 < rij12) {
    retval.distsq = rij2;
  } else {
    retval.distsq = rij12;
  }

  return retval;
}
*/


//
// find closest distance from a point to a line segment and return some related data
// logic taken from pointElemDistance2d
// does not work on simd data structures
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
  assert(blensq != 0); // Can't divide by 0
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
  auto start = std::chrono::steady_clock::now();

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,Dimensions> const& sn = _src.get_norm();
  const size_t nper = si.size() / _src.get_npanels();   // we can have >2 nodes per elem now
  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();

  // pre-compute the *node* normals
  std::array<Vector<S>,Dimensions> nn;
  for (size_t i=0; i<Dimensions; ++i) {
    nn[i].resize(_src.get_n());
    std::fill(nn[i].begin(), nn[i].end(), 0.0);
  }
  for (size_t j=0; j<_src.get_npanels(); ++j) {
    // nodes in this panel
    const Int ip0 = si[nper*j];
    const Int ip1 = si[nper*j+1];
    // add the normal to each of the nodes
    nn[0][ip0] += sn[0][j];
    nn[1][ip0] += sn[1][j];
    nn[0][ip1] += sn[0][j];
    nn[1][ip1] += sn[1][j];
  }
  //for (size_t j=0; j<_src.get_n(); ++j) {
  //  std::cout << "Node " << j << " at " << sx[0][j] << " " << sx[1][j] << " has norm " << nn[0][j] << " " << nn[1][j] << std::endl;
  //}

  size_t num_reflected = 0;
  //const S eps = 10.0*std::numeric_limits<S>::epsilon();

  #pragma omp parallel for reduction(+:num_reflected)
  for (int32_t i=0; i<(int32_t)_targ.get_n(); ++i) {
    S mindist = std::numeric_limits<S>::max();
    std::vector<ClosestReturn<S>> hits;

    // iterate and search for closest panel
    for (size_t j=0; j<_src.get_npanels(); ++j) {
      const Int jp0 = si[nper*j+0];
      const Int jp1 = si[nper*j+1];
      ClosestReturn<S> result = panel_point_distance<S>(sx[0][jp0], sx[1][jp0],
                                                        sx[0][jp1], sx[1][jp1],
                                                        tx[0][i],   tx[1][i]);

      //if (result.distsq < mindist - eps) {
      if (result.distsq < std::nextafter(mindist,0.0)) {
        // we blew the old one away
        mindist = result.distsq;
        if (result.disttype == node) {
          result.jidx = si[nper*j+result.jidx];
        } else {
          result.jidx = j;
        }
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      //} else if (result.distsq < mindist + eps) {
      } else if (result.distsq < std::nextafter(mindist, std::numeric_limits<S>::max())) {
        // we are effectively the same as the old closest
        if (result.disttype == node) {
          result.jidx = si[nper*j+result.jidx];
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
        normx += sn[0][j];
        normy += sn[1][j];
      } else {
        // hit a node, use the vector to the node (this could be backwards!)
        //const S bx = tx[0][i] - hits[k].cpx;
        //const S by = tx[1][i] - hits[k].cpy;
        //const S blen = 1.0 / std::sqrt(bx*bx + by*by);
        //normx += bx * blen;
        //normy += by * blen;
        // hit a node, use the cached node norm
        normx += nn[0][j];
        normy += nn[1][j];
      }
      cpx += hits[k].cpx;
      cpy += hits[k].cpy;
    }

    // if by some fluke the sum of the normals is zero (like when we initialize between two walls)
    //   just yank out one of the points, or choose the last
    if (normx*normx + normy*normy == 0.0) {
      const size_t j = hits[0].jidx;
      if (hits[0].disttype == panel) {
        // hit a panel, use the norm
        normx += sn[0][j];
        normy += sn[1][j];
      } else {
        normx += nn[0][j];
        normy += nn[1][j];
      }
      cpx = hits[0].cpx;
      cpy = hits[0].cpy;
    }
    // make sure we fixed it
    assert(std::sqrt(normx*normx + normy*normy) != 0.0); // Can't divide by 0

    // finish computing the mean norm and mean cp
    const S normilen = 1.0 / std::sqrt(normx*normx + normy*normy);
    normx *= normilen;
    normy *= normilen;
    assert((S)hits.size() != 0); // Can't divide by 0
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

  auto end = std::chrono::steady_clock::now();
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
  ct.push_back(std::make_tuple((S)(-nx-0.5)*dx, (S)0.0, (S)0.0));

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
    ct.push_back(std::make_tuple((S)(i+0.5)*dx, twgt, -tmom/twgt));
  }
  //std::cout << "  total weight " << twgt << std::endl;

  // add the last entry (keep all strength, set dshift to zero)
  ct.push_back(std::make_tuple((S)(nx+1.5)*dx, (S)1.0, (S)0.0));

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
        assert(std::get<0>(ct[i]) - std::get<0>(ct[i-1]) != 0); // Can't divide by zero
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
// return value is the amount of circulation removed
//
template <class S>
S clear_inner_panp2 (const int _method,
                     Surfaces<S> const & _src,
                     Points<S>& _targ,
                     const S _cutoff_mult,
                     const S _ips) {

  std::cout << "  Clearing" << _targ.to_string() << " from near" << _src.to_string() << std::endl;
  auto start = std::chrono::steady_clock::now();

  static bool made_cut_tables = false;
  static std::vector<std::tuple<S,S,S>> ct;
  if (not made_cut_tables) {
    ct = init_cut_tables<S>((S)0.1);
    made_cut_tables = true;
  }

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,Dimensions> const& sn = _src.get_norm();
  const size_t nper = si.size() / _src.get_npanels();   // we can have >2 nodes per elem now

  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();
  Vector<S>&                              ts = _targ.get_str();
  // if called on field points, there is no tr
  Vector<S>&                              tr = _targ.get_rad();
  const bool are_fldpts = tr.empty();

  // pre-compute the node normals
  std::array<Vector<S>,Dimensions> nn;
  for (size_t i=0; i<Dimensions; ++i) {
    nn[i].resize(_src.get_n());
    std::fill(nn[i].begin(), nn[i].end(), 0.0);
  }
  for (size_t j=0; j<_src.get_npanels(); ++j) {
    // nodes in this panel
    const Int ip0 = si[nper*j];
    const Int ip1 = si[nper*j+1];
    // add the normal to each of the nodes
    nn[0][ip0] += sn[0][j];
    nn[1][ip0] += sn[1][j];
    nn[0][ip1] += sn[0][j];
    nn[1][ip1] += sn[1][j];
  }

  if (_method==0 and not are_fldpts) {
    S this_circ = 0.0;
    for (int32_t i=0; i<(int32_t)_targ.get_n(); ++i) this_circ += ts[i];
    std::cout << "    circulation before: " << this_circ << std::endl;
  }

  size_t num_cropped = 0;
  S circ_removed = 0.0;
  //const S eps = 10.0*std::numeric_limits<S>::epsilon();

  // create array of flags - any moved particle will be tested again
  std::vector<bool> untested;
  untested.assign(_targ.get_n(), true);

  // iterate more than once to make sure particles get cleared from corners
  while (std::any_of(untested.begin(), untested.end(), [](bool x){return x;})) {

    #pragma omp parallel for reduction(+:num_cropped,circ_removed)
    for (int32_t i=0; i<(int32_t)_targ.get_n(); ++i) {

      if (untested[i]) {

        //S mindist = std::numeric_limits<S>::max();
        S mindist = 0.1*std::numeric_limits<S>::max();
        std::vector<ClosestReturn<S>> hits;

        // iterate and search for closest panel/node
        for (size_t j=0; j<_src.get_npanels(); ++j) {
          const Int jp0 = si[nper*j+0];
          const Int jp1 = si[nper*j+1];
          ClosestReturn<S> result = panel_point_distance<S>(sx[0][jp0], sx[1][jp0],
                                                            sx[0][jp1], sx[1][jp1],
                                                            tx[0][i],   tx[1][i]);

          //if (result.distsq < mindist - eps) {
          if (result.distsq < std::nextafter(mindist,0.0)) {
            // we blew the old one away
            mindist = result.distsq;
            if (result.disttype == node) {
              result.jidx = si[nper*j+result.jidx];
            } else {
              result.jidx = j;
            }
            hits.clear();
            hits.push_back(result);
            //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

          } else if (result.distsq < std::nextafter(mindist, std::numeric_limits<S>::max())) {
            // we are effectively the same as the old closest
            if (result.disttype == node) {
              result.jidx = si[nper*j+result.jidx];
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
        //assert(hits.size() > 0 && "No nearest neighbors");
        // this messes up long runs, ignore it (break out of this loop iteration)
        if (hits.size() == 0) continue;

        //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

        // init the mean normal and the mean contact point
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
            normx += sn[0][j];
            normy += sn[1][j];
            //std::cout << "    panel norm is " << (-by * blen) << " " << (bx * blen) << std::endl;
          } else {
            normx += nn[0][j];
            normy += nn[1][j];
            //std::cout << "    REAL norm is " << nn[0][hits[k].jidx] << " " << nn[1][hits[k].jidx] << std::endl;
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

        if (_method == 0) {
          if (dotp < this_radius) {
            //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << " because dotp " << dotp << " and norm " << norm[0] << " " << norm[1] << std::endl;

            // use precomputed table lookups for new position and remaining strength
            if (not are_fldpts) {
              assert(this_radius != 0); // Can't divide by 0
              const std::pair<S,S> entry = get_cut_entry(ct, dotp/this_radius);

              // ensure that this "reabsorbed" circulation is accounted for in BEM
              circ_removed += ts[i] * (1.0-std::get<0>(entry));

              // modify the particle in question
              ts[i] *= std::get<0>(entry);
              tx[0][i] += std::get<1>(entry) * this_radius * normx;
              tx[1][i] += std::get<1>(entry) * this_radius * normy;
            }

            num_cropped++;

          } else {
            // don't test this point again
            untested[i] = false;
          }

        } else if (_method == 1) {
          // simply push the particles until it is a certain distance from the body, keeping all strength
          //std::cout << "  particle " << i << " at " << tx[0][i] << " " << tx[1][i] << " has dotp " << dotp << std::endl;
          if (dotp < -0.001*_ips) {
            // modify the particle in question
            //std::cout << "  pushing " << tx[0][i] << " " << tx[1][i];
            tx[0][i] -= dotp * normx;
            tx[1][i] -= dotp * normy;
            //std::cout << " to " << tx[0][i] << " " << tx[1][i] << std::endl;
            num_cropped++;

          } else {
            // don't test this point again
            untested[i] = false;
          }
        }

      } // end if (untested)
    } // end loop over particles
  } // end loop over iterations

  // we did not resize the x array, so we don't need to touch the u array

  if (_method == 0) {
    std::cout << "    cropped " << num_cropped << " particles" << std::endl;
  } else if (_method == 1) {
    std::cout << "    pushed " << num_cropped << " particles" << std::endl;
  }

  if (_method==0 and not are_fldpts) {
    S this_circ = 0.0;
    for (int32_t i=0; i<(int32_t)_targ.get_n(); ++i) this_circ += ts[i];
    std::cout << "    circulation after: " << this_circ << std::endl;
    std::cout << "    removed: " << circ_removed << std::endl;
  }

  // flops count here is taken from reflect - might be different here
  const S flops = (_targ.get_n()+num_cropped) * (62.0 + 27.0*_src.get_npanels());

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    clear_inner_panp2:\t[%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return circ_removed;
}


//
// clean up by removing the innermost layer - the one that will be represented by boundary strengths
//
template <class S>
void clear_inner_layer(const int                _method,
                       std::vector<Collection>& _bdry,
                       std::vector<Collection>& _vort,
                       const S                  _cutoff_factor,
                       const S                  _ips) {

  // may need to do this multiple times to clear out concave zones!
  // this should only function when _vort is Points and _bdry is Surfaces
  for (auto &targ : _vort) {
    if (std::holds_alternative<Points<S>>(targ)) {
      Points<S>& pts = std::get<Points<S>>(targ);

      // don't push anything if they don't move themselves
      if (pts.get_movet() == lagrangian) {

        for (auto &src : _bdry) {
          if (std::holds_alternative<Surfaces<S>>(src)) {
            Surfaces<S>& surf = std::get<Surfaces<S>>(src);

            // call the specific panels-affect-points routine
            const S lost_circ = clear_inner_panp2<S>(_method, surf, pts, _cutoff_factor, _ips);

            // and tell the boundary collection that it reabsorbed that much
            surf.add_to_reabsorbed(lost_circ);
          }
        }
      }
    }

    // if lagrangian Surfaces (to be supported someday), we also need to push the nodes out
  }
}


//
// calculate the distance from every point in Points to the nearest part of Surfaces
//
template <class S>
Vector<S> get_nearest_distances(const Surfaces<S>& _src,
                                const Points<S>&   _targ) {

  std::cout << "  Finding distances to " << _targ.to_string() << " from " << _src.to_string() << std::endl;
  auto start = std::chrono::steady_clock::now();

  int32_t itest = -1;
  //int32_t itest = 46035;

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,Dimensions> const& tx = _targ.get_pos();
  const size_t nper = si.size() / _src.get_npanels();   // we can have >2 nodes per elem now

  // prepare the output vector
  Vector<S> distances;
  distances.resize(_targ.get_n());

  #pragma omp parallel for
  for (int32_t i=0; i<(int32_t)_targ.get_n(); ++i) {
    S mindist = std::numeric_limits<S>::max();
    if (itest==i) std::cout << "  testing pt " << itest << " at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // iterate and search for closest panel
    for (size_t j=0; j<_src.get_npanels(); ++j) {
      const Int jp0 = si[nper*j+0];
      const Int jp1 = si[nper*j+1];
      ClosestReturn<S> result = panel_point_distance<S>(sx[0][jp0], sx[1][jp0],
                                                        sx[0][jp1], sx[1][jp1],
                                                        tx[0][i],   tx[1][i]);
      if (itest==i) std::cout << "    dist to panel " << j << " at " << sx[0][jp0] << " " << sx[1][jp0] << " is " << std::sqrt(result.distsq) << std::endl;

      if (result.distsq < std::nextafter(mindist,S(0.0))) {
        // we blew the old one away
        mindist = result.distsq;
        if (itest==i) std::cout << "    new close panel " << j << " at " << std::sqrt(mindist) << " next lower float is " << std::sqrt(std::nextafter(mindist,S(0.0))) << std::endl;
      }
    }

    distances[i] = std::sqrt(mindist);
  }

  const S flops = _targ.get_n() * (1.0 + 27.0*_src.get_npanels());

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    get_nearest_distances: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);

  return distances;
}

