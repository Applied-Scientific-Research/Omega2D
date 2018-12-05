/*
 * Elements.h - abstract class for arrays of any computational elements in 2D
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <iostream>
#include <vector>
#include <memory>

//
// Abstract class for collection of elements: points/particles, lines/panels,
// sheets, volumes, etc.
//
// template parameters are data types for 'S'torage and 'I'ndicies
//
template <class S, class I>
class Elements {
public:
  explicit
  Elements() {}
  virtual ~Elements() {}

  virtual size_t get_n() = 0;
  virtual size_t const get_n() const = 0;
  virtual std::vector<S>& get_x() = 0;
  virtual std::vector<S> const& get_x() const = 0;
  //virtual std::vector<S>& get_x() = 0;
  virtual std::vector<S>& get_u() = 0;

  virtual std::unique_ptr<Elements<S,I>> clone() = 0;
  virtual void add_new(std::vector<S>&) = 0;
  virtual void reset_vels() = 0;
  virtual void update_max_str() = 0;
  virtual void scale_and_add_freestream(const std::array<S,2>&) = 0;
  virtual void step_in_place(const S) = 0;
  virtual void step_in_place(const S, std::vector<S> const&) = 0;
  virtual void increment_in_place() = 0;

  //virtual void debug(std::ostream& os) const = 0;
  //virtual std::string to_string() const = 0;
  //virtual std::vector<float> init_particles(float) const = 0;
  //virtual std::vector<float> step_particles(float) const = 0;

  // emit particles as vector of float4

  virtual void initGL(std::vector<float>&, float*, float*) = 0;
  virtual void updateGL() = 0;
  virtual void drawGL(std::vector<float>&, float*, float*) = 0;

protected:
private:
};

//std::ostream& operator<<(std::ostream &os, Elements const &elems);


/*
//
// Concrete class for a collection of 1-dimensional constant vortex sheets (segments)
//
class Segments : public Elements {
public:
  Segments(float _x, float _y, float _str, float _rad, float _soft)
    : Elements(_x, _y, _str),
      m_rad(_rad),
      m_softness(_soft)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

private:
  // the panel arrays - includes all panels from all bodies
  alignas(32) std::vector<S> x;
  alignas(32) std::vector< std::pair<I,I> > idx;
  // flow properties of panels
  alignas(32) std::vector<S> strengths;
};


//
// Concrete class for a collection of 2-dimensional constant vortex quads
//
class Quads : public Elements {
public:
  Quads(float _x, float _y, float _xsize, float _ysize, float _minstr, float _maxstr, int _num)
    : Elements(_x, _y),
      m_xsize(_xsize),
      m_num(_num)
    {}

  void debug(std::ostream& os) const override;
  std::string to_string() const override;
  std::vector<float> init_particles(float) const override;
  std::vector<float> step_particles(float) const override;

private:
  // the panel arrays - includes all panels from all bodies
  alignas(32) std::vector<S> x;
  //alignas(32) std::vector< std::array<I,4> > idx;
  alignas(32) std::vector<I> idx;
  // flow properties of panels
  alignas(32) std::vector<S> strengths;
};
*/

// possible subclassing of Quads to support higher-order quads (not constant strength)

