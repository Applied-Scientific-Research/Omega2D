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
#define _USE_MATH_DEFINES // Required by MSVC to define M_PI,etc. in <cmath>
#include <cmath>

enum ClosestType { panel, node };

template <class S> 
struct ClosestReturn { 
  size_t jpanel; 
  S distsq; 
  S cpx, cpy;
  ClosestType disttype; 
}; 

//
// find closest distance from a point to a line segment
// logic taken from pointElemDistance2d
//
template <class S>
ClosestReturn<S> panel_point_distance(const S sx0, const S sy0,
                                      const S sx1, const S sy1,
                                      const S tx,  const S ty) {

  ClosestReturn<S> retval;
  retval.jpanel = 0;

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
    retval.distsq = rij2;
    retval.cpx = sx0;
    retval.cpy = sy0;
  } else {
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

  size_t num_reflected = 0;
  const S eps = 10.0*std::numeric_limits<S>::epsilon();

  // accumulate results into targvel
  #pragma omp parallel for
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
        result.jpanel = j;
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      } else if (result.distsq < mindist + eps) {
        // we are effectively the same as the old closest
        result.jpanel = j;
        hits.push_back(result);
        //std::cout << "  THIS TIES THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;
      }
    }

    // dump out the hits
    if (false) {
    std::cout << "point " << i << " is " << tx[0][i] << " " << tx[1][i] << std::endl;
    for (auto & ahit: hits) {
      if (ahit.disttype == node) {
        std::cout << "  node " << ahit.jpanel << " is " << std::sqrt(ahit.distsq) << std::endl;
      } else {
        std::cout << "  panel " << ahit.jpanel << " is " << std::sqrt(ahit.distsq) << std::endl;
      }
    }
    }

    //std::cout << "  REFLECTING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // now look at the vector of return values and decide if we're under or above the panel!
    if (hits.size() == 1) {
      // this is easy if the closest is a panel!
      if (hits[0].disttype == panel) {
        // compare vector to normal vector
        const size_t j = hits[0].jpanel;
        const S bx = sx[0][si[2*j+1]] - sx[0][si[2*j]];
        const S by = sx[1][si[2*j+1]] - sx[1][si[2*j]];
        const S blen  = 1.0 / std::sqrt(bx*bx + by*by);
        const S normx = -by * blen;
        const S normy =  bx * blen;
        const S dist  = normx*(tx[0][i]-hits[0].cpx) + normy*(tx[1][i]-hits[0].cpy);
        //std::cout << "  dist is actually " << dist << std::endl;
        if (dist < 0.0) {
          // this point is under the panel - reflect it
          tx[0][i] -= 2.0*dist*normx;
          tx[1][i] -= 2.0*dist*normy;
          //std::cout << "    TYPE 1 TO " << tx[0][i] << " " << tx[1][i] << std::endl;
          num_reflected++;
        }
      } else {
        // this should never happen - any single hit must be a panel because every
        //   node effectively gets checked twice
        const S dist = std::sqrt(hits[0].distsq);
        // if the distance is large enough, we don't need to care
        if (dist < 0.1) {
          std::cout << "WARNING: point at " << tx[0][i] << " " << tx[1][i] << std::endl;
          std::cout << "  only hits one node on panel " << hits[0].jpanel << " with dist " << dist << std::endl;
        }
        // we cannot define a distance!
      }
    } else if (hits.size() == 2) {
      // if two hits, they should both be nodes, but in rare cases can be anything
      S normx = 0.0;
      S normy = 0.0;
      S cpx = 0.0;
      S cpy = 0.0;
      // find the mean normal and the mean contact point
      for (size_t k=0; k<hits.size(); ++k) {
        //std::cout << "    cp at " << hits[k].cpx << " " << hits[k].cpy << std::endl;
        const size_t j = hits[k].jpanel;
        const S bx = sx[0][si[2*j+1]] - sx[0][si[2*j]];
        const S by = sx[1][si[2*j+1]] - sx[1][si[2*j]];
        const S blen = 1.0 / std::sqrt(bx*bx + by*by);
        normx += -by * blen;
        normy +=  bx * blen;
        cpx += hits[k].cpx;
        cpy += hits[k].cpy;
      }
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
        tx[0][i] = hits[0].cpx + dist*normx;
        tx[1][i] = hits[0].cpy + dist*normy;
        //std::cout << " to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) << std::endl;
        //std::cout << "    TYPE 2 TO " << tx[0][i] << " " << tx[1][i] << std::endl;
        num_reflected++;
      }
    } else {
      // this should never happen!
      std::cout << "WARNING: point at " << tx[0][i] << " " << tx[1][i] << " hits " << hits.size() << " nodes/panels!" << std::endl;
    }
  }

  std::cout << "    reflected " << num_reflected << " particles" << std::endl;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    reflect_panp2:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
  //const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  //printf("    panels_affect_points: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// naive caller for the O(N^2) panel-particle reflection kernel
//
template <class S>
void clear_inner_panp2 (Surfaces<S> const & _src, Points<S>& _targ, const S _cutoff) {
  //std::cout << "  inside clear_inner_layer(Surfaces, Points)" << std::endl;
  std::cout << "  Clearing" << _targ.to_string() << " from near" << _src.to_string() << std::endl;
  auto start = std::chrono::system_clock::now();

  // get handles for the vectors
  std::array<Vector<S>,Dimensions> const& sx = _src.get_pos();
  std::vector<Int> const&                 si = _src.get_idx();
  std::array<Vector<S>,Dimensions>&       tx = _targ.get_pos();
  Vector<S>&                              ts = _targ.get_str();
  Vector<S>&                              tr = _targ.get_rad();

  size_t num_cropped = 0;
  const S eps = 10.0*std::numeric_limits<S>::epsilon();

  // accumulate results into targvel
  #pragma omp parallel for reduction(+:num_cropped)
  for (size_t i=0; i<_targ.get_n(); ++i) {

    // reserve space for oft-used temporaries
    S along[2], norm[2], thisnorm[2];
    S oopanlen, dotp, dist_sqrd, x0, y0, x1, y1;

    S near_dist_sqrd = std::numeric_limits<S>::max();
    // inear is the NODE id that the particle is closest to
    Int inear = std::numeric_limits<Int>::max();
    norm[0] = 0.0; norm[1] = 0.0;

    // iterate and search for closest panel/node
    for (size_t j=0; j<_src.get_n(); ++j) {

      x0 = sx[0][si[2*j]];
      y0 = sx[1][si[2*j]];
      x1 = sx[0][si[2*j+1]];
      y1 = sx[1][si[2*j+1]];
      // is the point closest to the panel line?

      // find vector along panel segment
      along[0] = x1 - x0;
      along[1] = y1 - y0;
      // one over the panel length is useful
      oopanlen = 1.0 / std::sqrt(along[0]*along[0] + along[1]*along[1]);
      // normalize along
      along[0] *= oopanlen;
      along[1] *= oopanlen;
      // and normal
      thisnorm[0] = -along[1];
      thisnorm[1] =  along[0];
      // find projection
      dotp = (tx[0][i]-x0)*along[0] + (tx[1][i]-y0)*along[1];
      dotp *= oopanlen;
      // dotp is now 0.0 to 1.0 for points ON the panel

      if (dotp < 0.0 + eps) {
        // particle is closer to first node
        dist_sqrd = std::pow(tx[0][i]-x0, 2) + pow(tx[1][i]-y0, 2);
        if (dist_sqrd < near_dist_sqrd - eps) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          inear = si[2*j];
          // replace the running normal
          norm[0] = thisnorm[0];
          norm[1] = thisnorm[1];
        } else if (dist_sqrd < near_dist_sqrd + eps) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += thisnorm[0];
          norm[1] += thisnorm[1];
        }

      } else if (dotp > 1.0 - eps) {
        // particle is closer to second node
        dist_sqrd = std::pow(tx[0][i]-x1, 2) + std::pow(tx[1][i]-y1, 2);
        if (dist_sqrd < near_dist_sqrd - eps) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          inear = si[2*j+1];
          // replace the running normal
          norm[0] = thisnorm[0];
          norm[1] = thisnorm[1];
        } else if (dist_sqrd < near_dist_sqrd + eps) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += thisnorm[0];
          norm[1] += thisnorm[1];
        }

      } else {
        // particle is closest to segment, not ends
        // find actual distance
        // first find closest point along segment
        const S cpx = x0 + dotp*along[0]/oopanlen;
        const S cpy = y0 + dotp*along[1]/oopanlen;
        // then compute distance
        dist_sqrd = std::pow(tx[0][i]-cpx, 2) + std::pow(tx[1][i]-cpy, 2);
        // and compare to saved
        if (dist_sqrd < near_dist_sqrd - eps) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          // and it doesn't matter which node we pick, the normal calculation is always the same
          inear = si[2*j];
          // replace the normal with this normal
          norm[0] = thisnorm[0];
          norm[1] = thisnorm[1];
        } else if (dist_sqrd < near_dist_sqrd + eps) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += thisnorm[0];
          norm[1] += thisnorm[1];
        }
      }

    } // end of loop over panels

    //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << std::endl;

    // compare the particle to the normal to see which side of the object it is on
    // then, find out which side the point is on (dot with normal)
    oopanlen = 1.0 / std::sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0] *= oopanlen;
    norm[1] *= oopanlen;
    dotp = (tx[0][i]-sx[0][inear])*norm[0] + (tx[1][i]-sx[1][inear])*norm[1];
    if (dotp > 0.0) {
      // particle is above panel
    } else {
      // particle is beneath panel
    }
    //std::cout << "    part is " << dotp << " away from node " << inear << std::endl;

    // HACK, weaken and push it out
    // eventually precompute table lookups for new position and remaining strength
    dotp -= _cutoff;
    if (dotp < tr[i]) {
      //std::cout << "  CLEARING pt at " << tx[0][i] << " " << tx[1][i] << " because dotp " << dotp << " and norm " << norm[0] << " " << norm[1] << std::endl;

      // smoothly vary the strength scaling factor from 0..1 over range dotp/vdelta = -1..1
      const S sfac = 0.5 * (1.0 + sin(0.5*M_PI*std::max((S)-1.0, dotp/tr[i])));
      ts[i] *= sfac;

      // move particle up to above the cutoff
      const S shiftd = tr[i] * 0.5 * (1.0 - dotp/tr[i]);
      //std::cout << "  PUSHING " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]);
      tx[0][i] += shiftd * norm[0];
      tx[1][i] += shiftd * norm[1];
      //std::cout << "    to " << std::sqrt(tx[0][i]*tx[0][i]+tx[1][i]*tx[1][i]) << " and weaken by " << sfac << std::endl;
      // do not change radius yet
      //std::cout << "    TO " << tx[0][i] << " " << tx[1][i] << " and weaken by " << sfac << std::endl;

      num_cropped++;
    }

  } // end loop over particles

  // we did not resize the x array, so we don't need to touch the u array

  std::cout << "    cropped " << num_cropped << " particles" << std::endl;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    clear_inner_panp2:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}


// helper struct for dispatching through a variant - don't need it
//struct ReflectVisitor {
  // source collection, target collection
//  void operator()(Points<float> const& src,   Points<float>& targ)   { points_affect_points<float>(src, targ); }
//  void operator()(Surfaces<float> const& src, Points<float>& targ)   { panels_affect_points<float>(src, targ); }
//  void operator()(Points<float> const& src,   Surfaces<float>& targ) { reflect_panp2<float>(src, targ); }
//  void operator()(Surfaces<float> const& src, Surfaces<float>& targ) { panels_affect_panels<float>(src, targ); }
//};
