/*
 * Core.h - Useful values for diffusion for various core functions
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
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

enum CoreType { gaussian, compactg };

//
// non-class templated functions (on real type)
//

// this is the per-axis second moment of a unit-volume core
template <class RT>
RT get_core_second_mom(const CoreType thiscore) {
  return (thiscore == CoreType::gaussian) ? (0.5) : (0.329727376);
}

template <class RT>
RT get_core_fourth_mom(const CoreType thiscore) {
  return (thiscore == CoreType::gaussian) ? (0.75) : (0.276933042);
}

template <class RT>
RT get_core_diffusion_time(const CoreType thiscore) {
  return (thiscore == CoreType::gaussian) ? (1.0) : (0.659454753);
}
