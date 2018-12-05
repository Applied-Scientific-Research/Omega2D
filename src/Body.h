/*
 * Body.h - class for an independent solid boundary
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#define _USE_MATH_DEFINES // Required by MSVC to define M_PI,etc. in <cmath>
#include <cmath>
#include <vector>
#include <array>

//------------------------------------------------------------------------
//
// A single rigid body
//
// Draws and moves itself, so needs to know its panels? no.
//
// movement subtypes DynamicBody, KinematicBody, StaticBody
// shape subtypes Square, Rectangle, Circle, Oval, Arbitrary, Plate, NACA4, Deforming, etc.
//
// how would we support inside-out bodies? say we're simulating flow inside a circle or box?
//

template <class S>
class Body {
public:
  Body();
  Body(const S, const S);
  virtual ~Body() {}

  // following are handled by the base class
  std::array<S,2> getPos();

  // following must be implemented by derived classes
  virtual std::vector<S> discretize(const S) = 0;

protected:
  // 2D position and velocity
  std::array<S,2> pos;
  std::array<S,2> vel;

  // angular position and velocity
  S apos;
  S avel;

private:
};

// delegating ctor
template <class S>
Body<S>::Body() :
  Body(0.0, 0.0)
  {}

// primary constructor
template <class S>
Body<S>::Body(const S _x, const S _y) :
  pos(std::array<S,2>({{_x, _y}})),
  vel(std::array<S,2>({{0.0, 0.0}})),
  apos(0.0),
  avel(0.0)
  {}

// getters/setters
template <class S>
std::array<S,2> Body<S>::getPos() {
  return pos;
}


//------------------------------------------------------------------------
//
// Subclass for true circles (really a specialization of oval?)
//
template <class S>
class Circle : public Body<S> {
public:
  Circle();
  Circle(const S);
  Circle(const S, const S, const S);

  S getRad();
  std::vector<S> discretize(const S);

protected:

private:
  // position and velocity are in Body
  S rad;
};

// delegating ctor
template <class S>
Circle<S>::Circle()
  : Circle(1.0)
  {}

// delegating ctor
template <class S>
Circle<S>::Circle(const S _diam)
  : Circle(0.0, 0.0, _diam)
  {}

// primary constructor
template <class S>
Circle<S>::Circle(const S _x, const S _y, const S _diam)
  : Body<S>(_x, _y),
    rad(0.5 * _diam)
  {}

// getters/setters
template <class S>
S Circle<S>::getRad() {
  return rad;
}

// discretize into an array of nodes
template <class S>
std::vector<S> Circle<S>::discretize(const S _ds) {

  // how many panels?
  size_t num_panels = std::min(10000, std::max(5, (int)(2.0*rad*M_PI / _ds)));

  printf("Adding circle with %ld panels\n", num_panels);

  // coordinates of nodes
  std::vector<S> retval(2*num_panels);

  // outside is to the left walking from one point to the next
  // so go CW around the circle starting at theta=0 (+x axis)
  for (size_t i=0; i<num_panels; i++) {
    retval[2*i]   = this->pos[0] + rad * cos(2.0 * M_PI * (S)i / (S)num_panels);
    retval[2*i+1] = this->pos[1] - rad * sin(2.0 * M_PI * (S)i / (S)num_panels);
  }

  //return std::array< std::vector<S>,2 >(x,y);
  return retval;
}


//------------------------------------------------------------------------
//
// Subclass for squares (really a specialization of rectangle?)
//
template <class S>
class Square : public Body<S> {
public:
  Square();
  Square(const S);
  Square(const S, const S, const S);

  S getSide();
  std::vector<S> discretize(const S);

protected:

private:
  // position and velocity are in Body
  S side;
};

// delegating ctor
template <class S>
Square<S>::Square()
  : Square(1.0)
  {}

// delegating ctor
template <class S>
Square<S>::Square(const S _side)
  : Square(0.0, 0.0, _side)
  {}

// primary constructor
template <class S>
Square<S>::Square(const S _x, const S _y, const S _side)
  : Body<S>(_x, _y),
    side(_side)
  {}

// getters/setters
template <class S>
S Square<S>::getSide() {
  return side;
}

//
// discretize into an array of nodes
// generate a number of particles in a square - use num to determine the square size
//
template <class S>
std::vector<S> Square<S>::discretize(const S _ds) {

  // how many panels?
  size_t num_panels = 4 * std::min(2500, std::max(1, (int)(side / _ds)));

  printf("Adding square with %ld panels\n", num_panels);

  // coordinates of nodes
  std::vector<S> x(2*num_panels);

  // outside is to the left walking from one point to the next
  // so go CW around the body
  size_t idx = 0;
  for (size_t i=0; i<num_panels/4; i++) {
    x[idx++] = this->pos[0] + side * -0.5;
    x[idx++] = this->pos[1] + side * (-0.5 + (S)i / (S)(num_panels/4));
  }
  for (size_t i=0; i<num_panels/4; i++) {
    x[idx++] = this->pos[0] + side * (-0.5 + (S)i / (S)(num_panels/4));
    x[idx++] = this->pos[1] + side * 0.5;
  }
  for (size_t i=0; i<num_panels/4; i++) {
    x[idx++] = this->pos[0] + side * 0.5;
    x[idx++] = this->pos[1] + side * (0.5 - (S)i / (S)(num_panels/4));
  }
  for (size_t i=0; i<num_panels/4; i++) {
    x[idx++] = this->pos[0] + side * (0.5 - (S)i / (S)(num_panels/4));
    x[idx++] = this->pos[1] + side * -0.5;
  }

  // return these nodes to the caller
  return x;
}

