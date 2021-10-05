/*
 * VectorHelper.h - Allows the same code to use std::vector of floats or Vc's SIMD float types
 *
 * (c)2018 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
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

#include <vector>

#ifdef USE_VC
#include <Vc/Vc>

// Must use Vc's allocator to ensure memory alignment
template <class S> using Vector = std::vector<S, Vc::Allocator<S>>;
// now we can use Vector<float> in code

//using VectorF = std::vector<float, Vc::Allocator<float>>;
//using VectorD = std::vector<double, Vc::Allocator<double>>;

// here's a way to make other vectors with the same number of entries as float_v
//using Vc::float_v;
//typedef Vc::SimdArray<double, float_v::size()> double_v;
//typedef Vc::SimdArray<std::uint16_t, float_v::size()> uint16_v;
//typedef Vc::SimdArray<std::uint32_t, float_v::size()> uint32_v;

// use this to safely convert a std::vector<S> to an array of Vc::Vector<S>
template <class S>
inline const Vc::Memory<Vc::Vector<S>> stdvec_to_vcvec (const Vector<S>& in, const S defaultval) {
    // create the new vector with memory to hold the given number of elements, plus buffer
    Vc::Memory<Vc::Vector<S>> out(in.size());
    //out.setZero();	// no need for this
    // because in.size() == out.entriesCount() always, we explicitly set the buffer region
    out.vector(out.vectorsCount()-1) = Vc::Vector<S>(defaultval);
    // now we copy the input vector
    for (size_t i=0; i<in.size(); ++i) out[i] = in[i];
    return out;
}

#else	// no Vc present, use stdlib instead

template <class S> using Vector = std::vector<S>;

//using VectorF = std::vector<float>;
//using VectorD = std::vector<double>;

#endif
