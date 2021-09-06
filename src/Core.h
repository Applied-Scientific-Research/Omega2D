/*
 * Core.h - Useful values for diffusion for various core functions
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
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
