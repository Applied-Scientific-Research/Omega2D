/*
 * Particles.h - a class for arrays of vortex particles in 2D
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Elements.h"
#include "OglHelper.h"
#include "ShaderHelper.h"
#include "glad.h"

#include <algorithm>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES // Required by MSVC to define M_PI,etc. in <cmath>
#include <cmath>


//
// Concrete class for a collection of 0-dimensional particles
//
template <class S, class I>
class Particles : public Elements<S,I> {
public:
  Particles()
    : Elements<S,I>(), max_strength(-1.0)
    {}

  size_t get_n() override { return x.size()/4; }
  size_t const get_n() const override { return x.size()/4; }
  std::vector<S>& get_x() override { return x; }
  std::vector<S> const& get_x() const override { return x; }
  std::vector<S>& get_u() override { return u; }
  //std::vector<S> const & get_u() const override { return u; }

  std::unique_ptr<Elements<S,I>> clone() override;
  void add_new(std::vector<S>&) override;
  void reset_vels() override;
  void update_max_str() override;
  void scale_and_add_freestream(const std::array<double,2>&) override;
  void step_in_place(const S) override;
  void step_in_place(const S, std::vector<S> const &) override;
  void increment_in_place() override;

  //void debug(std::ostream& os) const override;
  //std::string to_string() const override;
  //std::vector<float> init_particles(float) const override;
  //std::vector<float> step_particles(float) const override;

  // graphics calls
  void initGL(std::vector<float>&, float*, float*) override;
  void updateGL() override;
  void drawGL(std::vector<float>&, float*, float*) override;

protected:
  // state: x,y,str,rad
  alignas(8) std::vector<S> x;
  // differential: dx/dt, dy/dt, dstr/dt, drad/dt
  alignas(8) std::vector<S> u;
  //alignas(32) does not seem to work

private:
  // drawing data
  GLuint vao, vbo;
  GLuint blob_program;
  GLint projmat_attribute, projmat_attribute2, quad_attribute;
  GLint pos_color_attribute, neg_color_attribute, str_scale_attribute;
  float max_strength;
};


//
// needed for deep copy of Vorticity object, needed in multi-time-stepping
//
template <class S, class I>
std::unique_ptr<Elements<S,I>> Particles<S,I>::clone() {
  return std::make_unique<Particles>(*this);
}

//
// add the array of particles to the collection
//
// eventually this routine will determine which of the incoming particles are
//   most appropriate for this collection, maybe based on dt or size
//
template <class S, class I>
void Particles<S,I>::add_new(std::vector<S>& _in) {
  //std::cout << "  inside Particles::add_new with " << _in.size()/4 << std::endl;

  if (_in.size() == 0) return;
  assert(_in.size() % 4 == 0);

  std::cout << "  adding " << (_in.size()/4) << " particles to simulation...";

  // append to the end of the state vector
  x.insert(x.end(), _in.begin(), _in.end());

  // and resize the derivative vector with 0s
  size_t orig_size = u.size();
  u.resize(x.size());

  // copy the new data into the state array
  std::fill_n(&u[orig_size], _in.size(), 0.0);

  // update the maximum particle strength
  update_max_str();

  // and clean out the input vector
  _in.clear();

  std::cout << "there are " << x.size()/4 << " now" << std::endl;
}

//
// zero out the velocities
//
template <class S, class I>
void Particles<S,I>::reset_vels() {
  std::fill(u.begin(), u.end(), 0.0);
}

//
// reset saved max strength
//
template <class S, class I>
void Particles<S,I>::update_max_str() {
  float _thismax = 0.0;
  for (size_t i = 2; i < x.size(); i+=4) {
    const float _thisstr = std::abs(x[i]);
    if (_thisstr > _thismax) _thismax = _thisstr;
  }
  if (max_strength < 0.0) {
    max_strength = _thismax;
  } else {
    max_strength = 0.1*_thismax + 0.9*max_strength;
  }
  //std::cout << "  inside update_max_str with " << max_strength << std::endl;
}

//
// zero out the velocities
//
template <class S, class I>
void Particles<S,I>::scale_and_add_freestream(const std::array<double,2>& _fs) {
  assert(x.size() == u.size());

  for (size_t i=0; i<x.size(); i+=4) {
    u[i]   = _fs[0] + u[i]   / (2.0 * M_PI);
    u[i+1] = _fs[1] + u[i+1] / (2.0 * M_PI);
  }
}

//
// finish an Euler step
//
template <class S, class I>
void Particles<S,I>::step_in_place(const S _dt) {
  assert(x.size() == u.size());

  for (size_t i=0; i<x.size(); ++i) {
    x[i] += _dt * u[i];
  }
}

//
// finish a RK2 step
//
template <class S, class I>
void Particles<S,I>::step_in_place(const S _dt, std::vector<S> const & _u2) {
  assert(x.size() == u.size());
  assert(u.size() == _u2.size());

  for (size_t i=0; i<x.size(); ++i) {
    x[i] += 0.5 * _dt * (u[i] + _u2[i]);
    //std::cout << "  " << i/4 << " has two vels: " << u[i] << " and " << _u2[i] << std::endl;
  }
}

//
// increment radius and strength - NOT scaled by time, these are DISCRETE changes
//
template <class S, class I>
void Particles<S,I>::increment_in_place() {
  assert(x.size() == u.size());

  // here we add the strength change and replace the radius
  for (size_t i=0; i<x.size(); i+=4) {
    x[i+2] += u[i+2];	// add strength change
    x[i+3] = u[i+3];	// replace radius
  }
}


//
// OpenGL functions
//

// this gets done once - load the shaders, set up the vao
template <class S, class I>
void Particles<S, I>::initGL(std::vector<float>& _projmat,
                             float*              _poscolor,
                             float*              _negcolor) {

  // std::cout << "inside Particles.initGL" << std::endl;

  // Use a Vertex Array Object
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  // Create a Vector Buffer Object that will store the vertices on video memory
  glGenBuffers(1, &vbo);

  // Allocate space, but don't upload the data from CPU to GPU yet
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, 0, x.data(), GL_STATIC_DRAW);

  //
  // load and create the blob-drawing shader program
  //
  blob_program = create_particle_program();

  // Get the location of the attributes that enters in the vertex shader
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  GLint position_attribute = glGetAttribLocation(blob_program, "position");

  // Specify how the data for position can be accessed
  glVertexAttribPointer(position_attribute, 4, get_gl_type<S>, GL_FALSE, 0, 0);

  // Enable the attribute
  glEnableVertexAttribArray(position_attribute);

  // and tell it to advance two primitives per point
  glVertexAttribDivisor(position_attribute, 1);

  // upload the projection matrix
  projmat_attribute = glGetUniformLocation(blob_program, "Projection");
  glUniformMatrix4fv(projmat_attribute, 1, GL_FALSE, _projmat.data());

  // locate where the colors and color scales go
  pos_color_attribute = glGetUniformLocation(blob_program, "pos_color");
  neg_color_attribute = glGetUniformLocation(blob_program, "neg_color");
  str_scale_attribute = glGetUniformLocation(blob_program, "str_scale");

  // send the current values
  glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_poscolor);
  glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_negcolor);
  glUniform1f (str_scale_attribute, (const GLfloat)1.0);

  // and indicate the fragment color output
  glBindFragDataLocation(blob_program, 0, "frag_color");

  // Initialize the quad attributes
  std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
  GLuint qvbo;
  glGenBuffers(1, &qvbo);
  glBindBuffer(GL_ARRAY_BUFFER, qvbo);
  glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

  quad_attribute = glGetAttribLocation(blob_program, "quad_attr");
  glVertexAttribPointer(quad_attribute, 2, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(quad_attribute);

}


// this gets done every time we change the size of the positions array
template <class S, class I>
void Particles<S,I>::updateGL() {
  //std::cout << "inside Particles.updateGL" << std::endl;

  // has this been init'd yet?
  if (glIsVertexArray(vao) == GL_FALSE) return;

  if (x.size() > 0) {
    // Indicate and upload the data from CPU to GPU
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, x.size()*sizeof(S), x.data(), GL_DYNAMIC_DRAW);
  }
}

// OpenGL3 stuff to display points
template <class S, class I>
void Particles<S,I>::drawGL(std::vector<float>& _projmat,
                            float*              _poscolor,
                            float*              _negcolor) {

  //std::cout << "inside Particles.drawGL" << std::endl;
  //std::cout << "  inside display with n=" << n << std::endl;
  //std::cout << "  inside display with poscolor=" << (*_poscolor) << std::endl;
  //std::cout << "  inside display with maxstr=" << max_strength << std::endl;
  //std::cout << "  attrs " << projmat_attribute << " " << pos_color_attribute << " " << neg_color_attribute << std::endl;

  // has this been init'd yet?
  if (glIsVertexArray(vao) == GL_FALSE) {
    initGL(_projmat, _poscolor, _negcolor);
    updateGL();
  }

  if (x.size() > 0) {

    glBindVertexArray(vao);

    // get blending ready
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    // draw as colored clouds
    if (true) {
      glUseProgram(blob_program);

      glEnableVertexAttribArray(quad_attribute);

      // upload the current projection matrix
      glUniformMatrix4fv(projmat_attribute, 1, GL_FALSE, _projmat.data());

      // upload the current color values
      glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_poscolor);
      glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_negcolor);
      glUniform1f (str_scale_attribute, (const GLfloat)(0.4f/max_strength));

      // the one draw call here
      glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, get_n());
    }

    // return state
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
  }
}

