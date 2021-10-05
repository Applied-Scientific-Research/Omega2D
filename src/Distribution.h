/*
 * Distribution.h - Node distributions on a line
 *
 * (c)2020-21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *            Blake B Hillier <blakehillier@mac.com>
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

// this will eventually hold several simple patterns for node distributions
//   along a line: uniform, Chebyshev (1- and 2-sided), exponential (1- and 2-sided), etc.

#include <cmath>


// node distribution with dense points at each end
template <class T>
T chebyshev_node(T a, T b, T k, T n) {
  return (a+b)*0.5+(b-a)*0.5*std::cos((2*(n-k)-1)*M_PI*0.5/n);
}

// node distribution with tunable density at each end
// first, a measure of how many panels are needed given densities at each end and relative panel size
template <class T>
size_t chebyshev_node_2_count(const T dens_left, const T dens_right, const T rel_pan_size) {

  // first, compute theta bounds given node densities
  assert(dens_left > 0.0 && dens_left < 1.0 && dens_right > 0.0 && dens_right < 1.0 && "Panel densities out of range");
  const T theta_left = M_PI - std::asin(dens_left);
  const T theta_right = std::asin(dens_right);
  // theta range
  const T theta_range = theta_left - theta_right;
  // x range from these
  const T x_range = std::cos(theta_right) - std::cos(theta_left);
  // now, size of the middle panel
  const T x_panel = x_range * rel_pan_size;
  // convert that to an angle
  const T theta_panel = 2.0 * std::asin(0.5*x_panel);
  // finally, find panel count
  return std::max((int)3, (int)(0.5 + theta_range/theta_panel));
}

// then the mathematics itself to generate the node positions
template <class T>
T chebyshev_node_2(const T dens_left, const T dens_right, const size_t k, const size_t n) {

  // first, compute theta bounds given node densities
  assert(dens_left > 0.0 && dens_left < 1.0 && dens_right > 0.0 && dens_right < 1.0 && "Panel densities out of range");
  const T theta_left = M_PI - std::asin(dens_left);
  const T theta_right = std::asin(dens_right);
  // theta range
  const T theta_range = theta_left - theta_right;
  // x range from these
  const T x_range = std::cos(theta_right) - std::cos(theta_left);

  // now, find test theta from input indices
  const T theta = theta_right + theta_range*k/(T)n;

  // finally, compute and scale the return value (from 0..1)
  return (std::cos(theta_right) - std::cos(theta)) / x_range;
}

