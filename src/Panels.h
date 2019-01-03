/*
 * Panels.h - a class for discretized surfaces in 2D
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "BEM.h"
#include "Kernels.h"
#include "OglHelper.h"
#include "ShaderHelper.h"
#include "glad.h"

#define _USE_MATH_DEFINES

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

//
// Class to hold BEM parameters and temporaries
//
// templatized on 'S'torage type and 'I'ndex type
//
template <class S, class I>
class Panels {
public:
  Panels();
  Panels(const S);

  void panels_are_made();
  bool are_panels_made();
  S get_panel_size();
  I get_n() const { return (I)idx.size()/2; }
  I get_npanels() const { return (I)idx.size()/2; }
  I get_npoints() const { return (I)x.size()/2; }
  void write_first(const I);

  auto const& get_x() const { return x; }
  auto const& get_idx() const { return idx; }
  auto& get_vel() { return vel; }
  auto& get_strengths() { return strengths; }
  auto const& get_strengths() const { return strengths; }

  std::vector<S> get_collocation_points();
  void add_panels_from_closed_body(const std::vector<S>);
  void reset();
  void reset_vels();
  void scale_and_add_freestream(const std::array<double,2>&);
  void solve_bem();
  void set_rhs(const std::vector<S>&);
  void find_strengths();
  void set_max_strength();
  void add_this_influence(const std::vector<S>&, std::vector<S>&, size_t);

  std::vector<S> diffuse_onto(const S, const S, const S);

  // graphics calls
  void initGL(std::vector<float>&, float*, float*);
  void updateGL();
  void drawGL(std::vector<float>&, float*, float*);

protected:
  // influence functions
  //std::array<S,2> all_panels_affect_point(const S, const S, const std::array<S,2>);

private:
  // the panel arrays - includes all panels from all bodies - 2 entries per element
  alignas(32) std::vector<S> x;
  alignas(32) std::vector<I> idx;

  // remember which rows/columns were from separate bodies
  std::vector<std::pair<I,I>> body_idx;

  // flow properties of panels
  alignas(32) std::vector<S> strengths;

  // velocity influence on panels
  alignas(32) std::vector<S> vel;

  // to control resolution, set this to approximately minimum particle separation
  S min_panel_size;

  // the BEM equation and solver
  BEM<S,I> bem;

  // flag to indicate that we've made all of the panels
  bool all_panels_made;

  // drawing data
  GLuint vao, vbo, vbos, eab;
  GLuint shaderProgram;
  GLint projmat_unif;
  GLint pos_color_unif, neg_color_unif, str_scale_unif;
  float max_strength;
};

// delegating constructor
template <class S, class I>
Panels<S,I>::Panels() : Panels(0.1) {}

// primary constructor
template <class S, class I>
Panels<S,I>::Panels(const S _panel_size)
  : x(), idx(), body_idx(), strengths(), vel(),
    min_panel_size(_panel_size), bem(),
    all_panels_made(false)
  {}


//
// set/get status of panel data
//
template <class S, class I>
void Panels<S,I>::panels_are_made() {
  all_panels_made = true;
}
template <class S, class I>
bool Panels<S,I>::are_panels_made() {
  return all_panels_made;
}


//
// helps a caller determine how to size incoming panels
//
template <class S, class I>
S Panels<S,I>::get_panel_size() {
  return min_panel_size;
}

//
// send the caller a vector of collocation points
//
template <class S, class I>
std::vector<S> Panels<S,I>::get_collocation_points() {

  std::vector<S> cp;
  cp.resize(idx.size());

  for (size_t i=0; i<idx.size(); i+=2) {
    // collocation point for panel i
    cp[2*i+0] = 0.5 * (x[2*idx[i+1]]   + x[2*idx[i]]);
    cp[2*i+1] = 0.5 * (x[2*idx[i+1]+1] + x[2*idx[i]+1]);
  }

  return cp;
}


//
// send the caller a vector of diffused particles
//
template <class S, class I>
std::vector<S> Panels<S,I>::diffuse_onto(const S _dt, const S _re, const S _vdelta) {

  // create the vector of new points
  std::vector<S> np;
  np.resize(4*get_n());
  S along[2];
  S oopanlen;

  // the amount of time to diffuse does not change the amount of vorticity that diffuses,
  //   only the distance to the center of the vorticity (1st moment), which scales as
  //   the square root of dt
  // this is local h_nu times a factor (0.4*vdelta in 3D)
  const S fac = 1.683 * std::sqrt(_dt / _re);

  for (size_t i=0; i<get_n(); i++) {
    I id0 = idx[2*i];
    I id1 = idx[2*i+1];
    // start at center of panel
    np[4*i+0] = 0.5 * (x[2*id1]   + x[2*id0]);
    np[4*i+1] = 0.5 * (x[2*id1+1] + x[2*id0+1]);
    // push out a fixed distance
    along[0] = x[2*id1]   - x[2*id0];
    along[1] = x[2*id1+1] - x[2*id0+1];
    // one over the panel length is useful
    oopanlen = 1.0 / std::sqrt(along[0]*along[0] + along[1]*along[1]);
    // this assumes properly resolved, vdelta and dt
    np[4*i+0] += fac * -along[1] * oopanlen;
    np[4*i+1] += fac *  along[0] * oopanlen;
    // complete the element with a strength and radius
    np[4*i+2] = strengths[i] / oopanlen;
    np[4*i+3] = _vdelta;
    //std::cout << "  new part is " << np[4*i+0] << " " << np[4*i+1] << " " << np[4*i+2] << " " << np[4*i+3] << std::endl;
  }

  return np;
}

//
// accept a list of coordinates (x0,y0,x1,y1,...), assume it represents a
//   single closed body, and add the panels to the internal list
//
template <class S, class I>
void Panels<S,I>::add_panels_from_closed_body(const std::vector<S> _xp) {

  // must be an even number of entries
  assert( _xp.size()%2 == 0 );
  bool echo = false;

  const I new_n = _xp.size()/2;

  // add these to the global list
  I orig_n = x.size()/2;
  if (echo) printf("  before adding points, length is %ld\n", x.size());
  for (size_t i=0; i<2*new_n; ++i) {
    x.push_back(_xp[i]);
  }
  if (echo) printf("  after adding points, length is %ld\n", x.size());

  // and set indexes properly
  for (size_t i=0; i<(size_t)new_n-1; ++i) {
    idx.push_back( (I)(orig_n + i) );
    idx.push_back( (I)(orig_n + i + 1) );
  }
  idx.push_back( (I)(orig_n + new_n - 1) );
  idx.push_back( (I)(orig_n) );

  // echo complete panel list
  if (echo) for (size_t i=0; i<get_n(); ++i) {
    printf("    %ld  %d %d  %g %g  %g %g\n", i, idx[2*i], idx[2*i+1],
           x[2*idx[2*i]], x[2*idx[2*i]+1], x[2*idx[2*i+1]], x[2*idx[2*i+1]+1]);
  }

  // save an indicator of which panels correspond to which separate body
  body_idx.push_back( std::pair<I,I>(orig_n, orig_n + new_n));
  printf("  body %ld contains panels %d %d\n", (body_idx.size() - 1), orig_n, (orig_n + new_n));

  // ensure that all other arrays are properly sized
  strengths.resize(orig_n+new_n, 0.0);
  vel.resize(2*(orig_n+new_n), 0.0);

  // tell BEM that our A matrix is invalid
  bem.panels_changed();
}


//
// remove all sim data in preparation for a restart
//
template <class S, class I>
void Panels<S,I>::reset() {
  x.clear();
  idx.clear();
  body_idx.clear();
  strengths.clear();
  vel.clear();
  all_panels_made = false;
}


//
// reset velocities in preparation for finding the RHS of the equations
//
template <class S, class I>
void Panels<S,I>::reset_vels() {
  std::fill(vel.begin(), vel.end(), 0.0);
}


//
// finish up velocity evaluation
//
template <class S, class I>
void Panels<S,I>::scale_and_add_freestream(const std::array<double,2>& _fs) {
  assert(get_n() == vel.size()/2);

  for (size_t i=0; i<vel.size(); i+=2) {
    vel[i]   = _fs[0] + vel[i]   / (2.0 * M_PI);
    vel[i+1] = _fs[1] + vel[i+1] / (2.0 * M_PI);
    //std::cout << "  elem " << i/2 << " final vel is " << vel[i] << " " << vel[i+1] << std::endl;
  }
}


//
// solve for BEM strengths, given a set of particles
//
template <class S, class I>
void Panels<S,I>::solve_bem() {

  //std::cout << "After scaling and fs, cp vels are:" << std::endl;
  //for (size_t i=0; i<vel.size(); i+=2) {
  //  std::cout << "  " << i/2 << "   " << vel[i] << " " << vel[i+1] << std::endl;
  //}

  // use those vels to set the RHS of the matrix equation
  set_rhs(vel);

  // assemble the influence matrix and solve
  find_strengths();

  // compute and set peak strength magnitude
  set_max_strength();

  // debug print?
  //write_first(25);
}


//
// given a vector of velocities on collocation points, set the rhs vector
//
template <class S, class I>
void Panels<S,I>::set_rhs(const std::vector<S>& _vel) {

  assert(vel.size() == idx.size());

  // prepare b for writing
  static std::vector<S> b;
  b.resize(get_n());

  // loop over panels
  for (size_t i=0; i<get_n(); ++i) {

    // panel vector
    S panelx = x[2*idx[2*i+1]]   - x[2*idx[2*i]];
    S panely = x[2*idx[2*i+1]+1] - x[2*idx[2*i]+1];
    S panell = std::sqrt(panelx*panelx + panely*panely);

    // dot product of tangent with local velocity
    b[i] = -(_vel[2*i]*panelx + _vel[2*i+1]*panely) / panell;
    //std::cout << "  elem " << i << " rhs is " << b[i] << std::endl;
  }

  // send this to the BEM solver
  bem.set_rhs(b);
}

//
// set up and solve the Boundary Element Method problem
//
template <class S, class I>
void Panels<S,I>::find_strengths() {

  assert(all_panels_made);
  std::cout << "  inside find_strengths with " << get_n() << " panels" << std::endl;

  // initialize timers
  std::chrono::system_clock::time_point start, end;
  std::chrono::duration<double> elapsed_seconds;

  // update all normals and panel lengths

  // find influence matrix
  if (not bem.is_A_current()) {
    start = std::chrono::system_clock::now();
    bem.assemble_influence_matrix(x, get_idx());
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    printf("    make A matrix:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());
  }

  // use existing rhs

  // solve BEM and get the results
  start = std::chrono::system_clock::now();
  bem.solve();
  strengths = bem.getStrengths();
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  printf("    total solver time:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());
}


//
// find the maximum strength magnitude - it gets used when drawing the panels
//
template <class S, class I>
void Panels<S,I>::set_max_strength() {

  assert(strengths.size() == get_n());

  S maxstr = -1.0;

  for (size_t i=0; i<get_n(); ++i) {
    maxstr = std::max(maxstr, std::abs(strengths[i]));
  }

  max_strength = maxstr;
}


//
// echo some properties of the first _num panels
//
template <class S, class I>
void Panels<S,I>::write_first(const I _num) {

  I nwrite = _num;
  // calling write_first(-1); will write all panels
  if (nwrite < 0) nwrite = get_n();

  std::vector<S> b = bem.getRhs();

  // write info for up to _num or n panels
  printf("Panels id, nodes, center, rhs, and strength\n");
  for (size_t i=0; i<std::min(idx.size(), (size_t)nwrite); ++i) {
    printf("    %ld  %d %d  %g %g  %g %g\n", i, idx[2*i], idx[2*i+1],
           0.5*(x[2*idx[2*i]]+x[2*idx[2*i+1]]), 0.5*(x[2*idx[2*i]+1]+x[2*idx[2*i+1]+1]),
           b[i], strengths[i]);
  }
}

//
// compute influence of all panels on a vector of points
//
template <class S, class I>
void Panels<S,I>::add_this_influence(const std::vector<S>& _x, std::vector<S>& _u, size_t targstride) {

  // make sure we have enough room for the results
  assert(_x.size() == _u.size());
  assert(targstride > 1);
  assert(targstride < 17);

  // accumulate results into targvel
  //#pragma omp parallel for
  for (size_t i=0; i<_x.size(); i+=targstride) {

    // generate accumulator
    std::array<S,16> accum = {0.0};

    // iterate and accumulate
    for (size_t j=0; j<get_n(); j++) {
      //nbody_kernel(&srcpos[j], &x[i], accum.data());
      // influence of vortex panel j with unit circulation on given point
      auto nv = vortex_panel_affects_point<S,S>(x[2*idx[2*j]],  x[2*idx[2*j]+1],
                                                x[2*idx[2*j+1]], x[2*idx[2*j+1]+1],
                                                strengths[j], _x[i], _x[i+1]);
      // add to running total
      accum[0] += nv[0];
      accum[1] += nv[1];
    }

    // add to running sum
    for (size_t j=0; j<targstride; ++j) {
      _u[i+j] += accum[j];
    }
  }
}

//
// OpenGL calls, handled here
//
template <class S, class I>
void Panels<S,I>::initGL(std::vector<float>& _projmat,
                         float*              _poscolor,
                         float*              _negcolor) {
  //std::cout << "inside Panels.initGL" << std::endl;

  // Use a Vertex Array Object
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  // load and create the shader program
  shaderProgram = create_panel_program();

  // first, fiddle with the attributes

  // Create a Vector Buffer Object that will store the vertices on video memory
  glGenBuffers(1, &vbo);
  // Allocate space and upload the data from CPU to GPU - nothing for now
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, 0, x.data(), GL_STATIC_DRAW);

  // position data
  GLint position_attribute = glGetAttribLocation(shaderProgram, "position");
  // Specify how the data for position can be accessed
  glVertexAttribPointer(position_attribute, 2, get_gl_type<S>, GL_FALSE, 0, 0);
  // Enable the attribute
  glEnableVertexAttribArray(position_attribute);

  // load up the strengths also
  glGenBuffers(1, &vbos);
  glBindBuffer(GL_ARRAY_BUFFER, vbos);
  glBufferData(GL_ARRAY_BUFFER, 0, strengths.data(), GL_STATIC_DRAW);

  // strength data
  GLint strength_attribute = glGetAttribLocation(shaderProgram, "rawstr");
  glVertexAttribPointer(strength_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(strength_attribute);

  // do the same for the indices
  glGenBuffers(1, &eab);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eab);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, idx.data(), GL_STATIC_DRAW);

  // now the uniforms

  // Get the location of the attributes that enters in the vertex shader
  projmat_unif = glGetUniformLocation(shaderProgram, "Projection");

  // upload the projection matrix
  glUniformMatrix4fv(projmat_unif, 1, GL_FALSE, _projmat.data());

  // locate where the colors and color scales go
  pos_color_unif = glGetUniformLocation(shaderProgram, "pos_color");
  neg_color_unif = glGetUniformLocation(shaderProgram, "neg_color");
  str_scale_unif = glGetUniformLocation(shaderProgram, "str_scale");

  // send the current values
  glUniform4fv(pos_color_unif, 1, (const GLfloat *)_poscolor);
  glUniform4fv(neg_color_unif, 1, (const GLfloat *)_negcolor);
  glUniform1f (str_scale_unif, (const GLfloat)1.0);

  // and indicate the fragment color output
  glBindFragDataLocation(shaderProgram, 0, "frag_color");
}

template <class S, class I>
void Panels<S,I>::updateGL() {
  //std::cout << "inside Panels.updateGL" << std::endl;

  // has this been init'd yet?
  if (glIsVertexArray(vao) == GL_FALSE) return;

  if (get_n() > 0) {
    // Why can't I just upload this once, back in initGL?
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(S)*x.size(), x.data(), GL_DYNAMIC_DRAW);
    //std::cout << "X " << (S)x[0] << " " << (S)x[1] << " " << (I)idx[0] << " " << (I)idx[1] << std::endl;

    glBindBuffer(GL_ARRAY_BUFFER, vbos);
    glBufferData(GL_ARRAY_BUFFER, sizeof(S)*strengths.size(), strengths.data(), GL_DYNAMIC_DRAW);
    //std::cout << "  str " << (S)strengths[0] << " " << (S)strengths[1] << " " << (S)strengths[2] << " " << (S)strengths[3] << std::endl;

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eab);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(I)*idx.size(), idx.data(), GL_DYNAMIC_DRAW);
  }
}

template <class S, class I>
void Panels<S,I>::drawGL(std::vector<float>& _projmat,
                         float* _poscolor, float* _negcolor) {
  //std::cout << "inside Panels.drawGL" << std::endl;

  // has this been init'd yet?
  if (glIsVertexArray(vao) == GL_FALSE) {
    initGL(_projmat, _poscolor, _negcolor);
    updateGL();
  }

  if (get_n() > 0) {
    glUseProgram(shaderProgram);

    // upload the current projection matrix
    glUniformMatrix4fv(projmat_unif, 1, GL_FALSE, _projmat.data());
  
    // upload the current color values
    glUniform4fv(pos_color_unif, 1, (const GLfloat *)_poscolor);
    glUniform4fv(neg_color_unif, 1, (const GLfloat *)_negcolor);
    glUniform1f (str_scale_unif, (const GLfloat)(max_strength));

    // get blending ready
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    glLineWidth(2.0);
    glBindVertexArray(vao);

    // Why can't I just upload this once, back in initGL?
    //glBindBuffer(GL_ARRAY_BUFFER, vbo);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(S)*x.size(), x.data(), GL_DYNAMIC_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, 0, x.data(), GL_DYNAMIC_DRAW);
    //std::cout << "X " << (S)x[0] << " " << (S)x[1] << " " << (I)idx[0] << " " << (I)idx[1] << std::endl;

    //glBindBuffer(GL_ARRAY_BUFFER, vbos);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(S)*strengths.size(), strengths.data(), GL_DYNAMIC_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, 0, strengths.data(), GL_DYNAMIC_DRAW);
    //std::cout << "  str " << (S)strengths[0] << " " << (S)strengths[1] << " " << (S)strengths[2] << " " << (S)strengths[3] << std::endl;

    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eab);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(I)*idx.size(), idx.data(), GL_STATIC_DRAW);

    //glDrawArrays(GL_POINTS, 0, 100);
    //glDrawArrays(GL_LINE_LOOP, 0, get_npoints());
    glDrawElements(GL_LINES, idx.size(), get_gl_type<I>, 0);
    
    // return state
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
  }
}

