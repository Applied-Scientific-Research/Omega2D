/*
 * Reflect.h - Non-class particle-panel reflecting operation
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Boundaries.h"
#include "Particles.h"
#include "Vorticity.h"

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
template <class S, class I>
void reflect_panp (Panels<S,I> const& _src, Elements<S,I>& _targ) {

  // get handles for the vectors
  std::vector<S> const& sx = _src.get_x();
  std::vector<I> const& si = _src.get_idx();
  std::vector<S>&       tx = _targ.get_x();		// contains 0=x, 1=y, 2=str, 3=rad

  size_t num_reflected = 0;

  // accumulate results into targvel
  #pragma omp parallel for
  for (size_t i=0; i<_targ.get_n(); ++i) {
    //S result = std::numeric_limits<S>::max();
    S mindist = std::numeric_limits<S>::max();
    std::vector<ClosestReturn<S>> hits;

    // iterate and search for closest panel
    for (size_t j=0; j<_src.get_n(); ++j) {
      ClosestReturn<S> result = panel_point_distance<S>(sx[2*si[2*j]],   sx[2*si[2*j]+1],
                                                        sx[2*si[2*j+1]], sx[2*si[2*j+1]+1],
                                                        tx[4*i+0], tx[4*i+1]);

      if (result.distsq < mindist - std::numeric_limits<S>::epsilon()) {
        // we blew the old one away
        mindist = result.distsq;
        result.jpanel = j;
        hits.clear();
        hits.push_back(result);
        //std::cout << "  THIS BLOWS AWAY THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;

      } else if (result.distsq < mindist + std::numeric_limits<S>::epsilon()) {
        // we are effectively the same as the old closest
        result.jpanel = j;
        hits.push_back(result);
        //std::cout << "  THIS TIES THE CLOSEST, AT " << std::sqrt(mindist) << std::endl;
      }
    }

    // dump out the hits
    //std::cout << "point is " << tx[4*i+0] << " " << tx[4*i+1] << std::endl;
    //for (auto & ahit: hits) {
    //  if (ahit.disttype == node) {
    //    std::cout << "  node " << ahit.jpanel << " is " << std::sqrt(ahit.distsq) << std::endl;
    //  } else {
    //    std::cout << "  panel " << ahit.jpanel << " is " << std::sqrt(ahit.distsq) << std::endl;
    //  }
    //}

    //std::cout << "  REFLECTING pt at " << tx[4*i+0] << " " << tx[4*i+1] << std::endl;

    // now look at the vector of return values and decide if we're under or above the panel!
    if (hits.size() == 1) {
      // this is easy if the closest is a panel!
      if (hits[0].disttype == panel) {
        // compare vector to normal vector
        const size_t j = hits[0].jpanel;
        const S bx = sx[2*si[2*j+1]] - sx[2*si[2*j]];
        const S by = sx[2*si[2*j+1]+1] - sx[2*si[2*j]+1];
        const S blen  = 1.0 / std::sqrt(bx*bx + by*by);
        const S normx = -by * blen;
        const S normy =  bx * blen;
        const S dist  = normx*(tx[4*i+0]-hits[0].cpx) + normy*(tx[4*i+1]-hits[0].cpy);
        //std::cout << "  dist is actually " << dist << std::endl;
        if (dist < 0.0) {
          // this point is under the panel - reflect it
          tx[4*i+0] -= 2.0*dist*normx;
          tx[4*i+1] -= 2.0*dist*normy;
          //std::cout << "    TYPE 1 TO " << tx[4*i+0] << " " << tx[4*i+1] << std::endl;
          num_reflected++;
        }
      } else {
        // this should never happen - any single hit must be a panel because every
        //   node effectively gets checked twice
        const S dist = std::sqrt(hits[0].distsq);
        // if the distance is large enough, we don't need to care
        if (dist < 2.0 * tx[4*i+3]) {
          std::cout << "WARNING: point at " << tx[4*i+0] << " " << tx[4*i+1] << std::endl;
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
        const S bx = sx[2*si[2*j+1]] - sx[2*si[2*j]];
        const S by = sx[2*si[2*j+1]+1] - sx[2*si[2*j]+1];
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
      const S dotp = normx*(tx[4*i+0]-cpx) + normy*(tx[4*i+1]-cpy);
      if (dotp < 0.0) {
        // this point is under the panel - reflect it off entry 0
        // this is reasonable for most cases, except very sharp angles between adjacent panels
        const S dist = std::sqrt(hits[0].distsq);
        //std::cout << "  REFLECTING " << std::sqrt(tx[4*i+0]*tx[4*i+0]+tx[4*i+1]*tx[4*i+1]);
        tx[4*i+0] = hits[0].cpx + dist*normx;
        tx[4*i+1] = hits[0].cpy + dist*normy;
        //std::cout << " to " << std::sqrt(tx[4*i+0]*tx[4*i+0]+tx[4*i+1]*tx[4*i+1]) << std::endl;
        //std::cout << "    TYPE 2 TO " << tx[4*i+0] << " " << tx[4*i+1] << std::endl;
        num_reflected++;
      }
    } else {
      // this should never happen!
      std::cout << "WARNING: point at " << tx[4*i+0] << " " << tx[4*i+1] << " hits " << hits.size() << " nodes/panels!" << std::endl;
    }
  }

  std::cout << "    reflected " << num_reflected << " particles" << std::endl;
}

//
// High-level driver for all-affects-all
//
// boundaries are the source (const) and particles are the targets (modified)
//
template <class S, class I>
void reflect (Boundaries<S,I> const& _src, Vorticity<S,I>& _targ) {
  std::cout << "  inside reflect(Boundaries, Vorticity)" << std::endl;

  // just one set of panels for all boundaries now
  Panels<S,I> const& src_elem = _src.get_panels();

  // iterate over all sets of target vorticities (currently only one Particles object)
  for (auto& targ_elem : _targ.get_collections()) {
    std::cout << "    computing reflection of " << targ_elem->get_n() << " particles off " << src_elem.get_n() << " panels" << std::endl;

    // here's the problem: the routine here doesn't know what type each of these is!
    // HACK - let's assume it's always Particles and Panels
    reflect_panp<S,I>(src_elem, *targ_elem);
  }
}


//
// naive caller for the O(N^2) panel-particle reflection kernel
//
template <class S, class I>
void clear_inner_panp (Panels<S,I> const & _src, Elements<S,I>& _targ, const S _cutoff) {

  // get handles for the vectors
  std::vector<S> const&   x = _src.get_x();
  std::vector<I> const& idx = _src.get_idx();
  std::vector<S>&      pmod = _targ.get_x();	// contains 0=x, 1=y, 2=str, 3=rad

  size_t num_cropped = 0;

  // accumulate results into targvel
  #pragma omp parallel for reduction(+:num_cropped)
  for (size_t i=0; i<_targ.get_n(); ++i) {

    // reserve space for oft-used temporaries
    S along[2], norm[2], thisnorm[2];
    S oopanlen, dotp, dist_sqrd, x0, y0, x1, y1;

    S near_dist_sqrd = std::numeric_limits<S>::max();
    // inear is the NODE id that the particle is closest to
    I inear = std::numeric_limits<I>::max();
    norm[0] = 0.0; norm[1] = 0.0;

    // iterate and search for closest panel/node
    for (size_t j=0; j<_src.get_n(); ++j) {

      x0 = x[2*idx[2*j]];
      y0 = x[2*idx[2*j]+1];
      x1 = x[2*idx[2*j+1]];
      y1 = x[2*idx[2*j+1]+1];
      // is the point closest to the panel line?

      // find vector along panel segment
      along[0] = x1 - x0;
      along[1] = y1 - y0;
      // one over the panel length is useful
      oopanlen = 1.0 / std::sqrt(along[0]*along[0] + along[1]*along[1]);
      // and normal
      thisnorm[0] = -along[1] * oopanlen;
      thisnorm[1] =  along[0] * oopanlen;
      // find projection
      dotp = (pmod[4*i]-x0)*along[0] + (pmod[4*i+1]-y0)*along[1];
      dotp *= oopanlen*oopanlen;

      if (dotp < 0.0) {
        // particle is closer to first node
        dist_sqrd = std::pow(pmod[4*i]-x0, 2) + pow(pmod[4*i+1]-y0, 2);
        if (dist_sqrd < near_dist_sqrd - std::numeric_limits<S>::epsilon()) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          inear = idx[2*j];
          // replace the running normal
          norm[0] = thisnorm[0];
          norm[1] = thisnorm[1];
        } else if (dist_sqrd < near_dist_sqrd + std::numeric_limits<S>::epsilon()) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += thisnorm[0];
          norm[1] += thisnorm[1];
        }

      } else if (dotp > 1.0) {
        // particle is closer to second node
        dist_sqrd = std::pow(pmod[4*i]-x1, 2) + std::pow(pmod[4*i+1]-y1, 2);
        if (dist_sqrd < near_dist_sqrd - std::numeric_limits<S>::epsilon()) {
          // point is clearly the closest
          near_dist_sqrd = dist_sqrd;
          inear = idx[2*j+1];
          // replace the running normal
          norm[0] = thisnorm[0];
          norm[1] = thisnorm[1];
        } else if (dist_sqrd < near_dist_sqrd + std::numeric_limits<S>::epsilon()) {
          // point is just as close as another point
          // add the normal to the running sum
          norm[0] += thisnorm[0];
          norm[1] += thisnorm[1];
        }

      } else {
        // particle is closest to segment, not ends
        // find actual distance
        // first find closest point along segment
        along[0] = x0 + dotp*along[0];
        along[1] = y0 + dotp*along[1];
        // then compute distance
        dist_sqrd = std::pow(pmod[4*i]-along[0], 2) + std::pow(pmod[4*i+1]-along[1], 2);
        // and compare to saved
        if (dist_sqrd < near_dist_sqrd + std::numeric_limits<S>::epsilon()) {
          // even if its a tie, this one wins
          near_dist_sqrd = dist_sqrd;
          // and it doesn't matter which node we pick, the normal calculation is always the same
          inear = idx[2*j];
          // replace the normal with this normal
          norm[0] = thisnorm[0];
          norm[1] = thisnorm[1];
        }
      }

    } // end of loop over panels

    //std::cout << "  CLEARING pt at " << pmod[4*i+0] << " " << pmod[4*i+1] << std::endl;

    // compare the particle to the normal to see which side of the object it is on
    // then, find out which side the point is on (dot with normal)
    oopanlen = 1.0 / std::sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0] *= oopanlen;
    norm[1] *= oopanlen;
    dotp = (pmod[4*i]-x[2*inear])*norm[0] + (pmod[4*i+1]-x[2*inear+1])*norm[1];
    if (dotp > 0.0) {
      // particle is above panel
    } else {
      // particle is beneath panel
    }
    //std::cout << "  part " << i/4 << " is " << dotp << " away from node " << inear << std::endl;

    // HACK, weaken and push it out
    // evantually precompute table lookups for new position and remaining strength
    dotp -= _cutoff;
    if (dotp < pmod[4*i+3]) {
      // smoothly vary the strength scaling factor from 0..1 over range dotp/vdelta = -1..1
      const S sfac = 0.5 * (1.0 + sin(0.5*M_PI*std::max((S)-1.0, dotp/pmod[4*i+3])));
      pmod[4*i+2] *= sfac;

      // move particle up to above the cutoff
      const S shiftd = pmod[4*i+3] * 0.5 * (1.0 - dotp/pmod[4*i+3]);
      //std::cout << "  PUSHING " << std::sqrt(pmod[4*i+0]*pmod[4*i+0]+pmod[4*i+1]*pmod[4*i+1]);
      pmod[4*i+0] += shiftd * norm[0];
      pmod[4*i+1] += shiftd * norm[1];
      //std::cout << " to " << std::sqrt(pmod[4*i+0]*pmod[4*i+0]+pmod[4*i+1]*pmod[4*i+1]) << std::endl;
      // do not change radius yet
      //pmod[i+2] = 0.0;
      //std::cout << "    TO " << pmod[4*i+0] << " " << pmod[4*i+1] << " and weaken by " << sfac << std::endl;

      num_cropped++;
    }

  } // end loop over particles

  // we did not resize the x array, so we don't need to touch the u array

  std::cout << "    cropped " << num_cropped << " particles" << std::endl;
}


//
// Crop all near-body particles and replace them with their remainders
//
template <class S, class I>
void clear_inner_layer (Boundaries<S,I> const& _src,
                        Vorticity<S,I>&        _targ,
                        const S                _cutoff) {
  std::cout << "  inside clear_inner_layer(Boundaries, Vorticity)" << std::endl;

  // just one set of panels for all boundaries now
  Panels<S,I> const& src_elem = _src.get_panels();

  // iterate over all sets of target vorticities (currently only one Particles object)
  for (auto& targ_elem : _targ.get_collections()) {
    std::cout << "    clearing innermost parts of " << targ_elem->get_n() << " particles from " << src_elem.get_n() << " panels" << std::endl;

    // here's the problem: the routine here doesn't know what type each of these is!
    // HACK - let's assume it's always Particles and Panels
    clear_inner_panp<S,I>(src_elem, *targ_elem, _cutoff);
  }
}

