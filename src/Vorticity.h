/*
 * Vorticity.h - template class representing collections of 2d vortex particles
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Particles.h"

#include <array>
#include <iostream>
#include <memory>
#include <vector>

//
// Class containing all participating non-boundary-condition vorticity
//
// template parameters are data types for 'S'torage and 'I'ndicies
//
template <class S, class I>
class Vorticity {
public:
  Vorticity() {}
  Vorticity(Vorticity const&);
  ~Vorticity();

  size_t num_collections() const;
  size_t get_n() const;
  auto& get_collections() { return elems; }
  auto const& get_collections() const { return elems; }
  void reset();
  void add_new(std::vector<S>&);
  //void add_new(std::vector<S>&, std::vector<std::pair<I,I> >&);
  void reset_vels();
  void update_max_str();
  void scale_and_add_freestream(std::array<double,2> const&);
  void step_in_place(const S);
  void step_in_place(const S, Vorticity<S,I> const&);

  // all influence-finding is done outside of this class

  //virtual void debug(std::ostream& os) const = 0;
  //virtual std::string to_string() const = 0;

  // graphics pass-through calls
  void initGL(std::vector<float>&, float*, float*);
  void updateGL();
  void drawGL(std::vector<float>&, float*, float*);

protected:

private:
  // one entry for each collection of participating vorticity-carriers (like Particles)
  std::vector<std::unique_ptr<Elements<S,I>>> elems;
};


//
// custom copy constructor because we need copies of objects pointed to my std::unique_ptr
//
template <class S, class I>
Vorticity<S,I>::Vorticity(Vorticity const & src) {
  elems.reserve(src.elems.size());
  for (const auto & e : src.elems) {
    elems.push_back(e->clone());
  }
}

// the destructor needs to be deep, too
template <class S, class I>
Vorticity<S,I>::~Vorticity() {
  for (auto& e : elems) {
    e.reset();
  }
  elems.clear();
}

//
// will we ever need to know how many collections we have?
//
template <class S, class I>
size_t Vorticity<S,I>::num_collections() const {
  return elems.size();
}

//
// add up the total number of elements
//
template <class S, class I>
size_t Vorticity<S,I>::get_n() const {
  size_t n = 0;
  for (auto const& v : elems) {
    n += v->get_n();
  }
  return n;
}

//
// lose everything
//
template <class S, class I>
void Vorticity<S,I>::reset() {
  elems.clear();
}

//
// add the given vector of (x,y,str,rad) to the matching collection of Particles
//
template <class S, class I>
void Vorticity<S,I>::add_new(std::vector<S>& _newelems) {

  // each collection will add whichever are appropriate
  for (auto const& v : elems) {
    v->add_new(_newelems);
  }

  // if it was added to none, or no collections exist, or not all points got assigned
  if (_newelems.size() > 0) {
    // make a new collection
    elems.emplace_back(std::unique_ptr<Particles<S,I> >(new Particles<S,I>()));

    // and add all remaining points to it
    elems.back()->add_new(_newelems);
  }
}

/*
//
// add the given vector of (x,y,str,rad) to the matching collection of Particles
//
template <class S, class I>
void Vorticity<S,I>::add_new(std::vector<S>& _newnodes, std::vector<std::pair<I,I> >& _newelems) {
  // each collection will add whichever are appropriate
  for (auto const& v : elems) {
//    v.add_new(_newnodes,_newelems);
  }
}
*/

//
// in preparation for velocity-finding, zero the velocities on each element
//
template <class S, class I>
void Vorticity<S,I>::reset_vels() {
  for (auto const& v : elems) {
    v->reset_vels();
  }
}

//
// must reset known max str before drawing
//
template <class S, class I>
void Vorticity<S,I>::update_max_str() {
  for (auto const& v : elems) {
    v->update_max_str();
  }
}

//
// after velocity-finding, scale by 1/2pi and add the freestream
//
template <class S, class I>
void Vorticity<S,I>::scale_and_add_freestream(const std::array<double,2>& _fs) {
  for (auto const& v : elems) {
    v->scale_and_add_freestream(_fs);
  }
}

//
// apply the self-velocity to update the self-positions
//
template <class S, class I>
void Vorticity<S,I>::step_in_place(const S _dt) {
  for (auto const& v : elems) {
    v->step_in_place(_dt);
  }
}

//
// apply a combination of the self-velocity and the update velocity to update the self-positions
//
template <class S, class I>
void Vorticity<S,I>::step_in_place(const S _dt, Vorticity<S,I> const & _other) {
  auto& otherelems = _other.get_collections();
  assert(elems.size() == otherelems.size());
  // loop over both collections simultaneusly
  for (size_t i=0; i<elems.size(); ++i) {
    elems[i]->step_in_place(_dt, otherelems[i]->get_u());
  }
}

//
// OpenGL calls, pass on to children
//
template <class S, class I>
void Vorticity<S,I>::initGL(std::vector<float>& _projmat,
                            float* _poscolor, float* _negcolor) {
  //std::cout << "inside Vorticity.initGL" << std::endl;
  for (auto const& v : elems) {
    v->initGL(_projmat, _poscolor, _negcolor);
  }
}

template <class S, class I>
void Vorticity<S,I>::updateGL() {
  for (auto const& v : elems) {
    v->updateGL();
  }
}

template <class S, class I>
void Vorticity<S,I>::drawGL(std::vector<float>& _projmat,
                            float* _poscolor, float* _negcolor) {
  //std::cout << "inside Vorticity.drawGL" << std::endl;
  for (auto const& v : elems) {
    v->drawGL(_projmat, _poscolor, _negcolor);
  }
}

//std::ostream& operator<<(std::ostream& os, Elements const& elems);

