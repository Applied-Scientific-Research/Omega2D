/*
 * Points.h - Specialized class for 2D points
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "ElementBase.h"

#ifdef USE_GL
#include "RenderParams.h"
#include "OglHelper.h"
#include "ShaderHelper.h"
#include "glad.h"
#endif

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <optional>
#include <random>
#include <cassert>
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>


// 0-D elements
template <class S>
class Points: public ElementBase<S> {
public:
  // flexible constructor - use input 4*n vector (x, y, s, r)
  //                         or input 2*n vector (x, y)
  Points(const std::vector<S>& _in,
         const elem_t _e,
         const move_t _m,
         std::shared_ptr<Body> _bp)
    : ElementBase<S>(0, _e, _m, _bp),
      max_strength(-1.0) {

    size_t nper = 4;
    if (_e == inert) {
      nper = 2;
      std::cout << "  new collection with " << (_in.size()/nper) << " tracers..." << std::endl;
    } else {
      nper = 4;
      std::cout << "  new collection with " << (_in.size()/nper) << " vortons..." << std::endl;
    }

    // need to reset the base class n
    this->n = _in.size()/nper;

    // make sure we have a complete input vector
    assert(_in.size() % nper == 0);

    // this initialization specific to Points
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        this->x[d][i] = _in[nper*i+d];
      }
    }

    // save untransformed positions if we are given a Body pointer
    if (_bp) {
      this->ux = this->x;
    }

    if (_e == inert) {
      // field points need no radius, but we must set one anyway so that vel evals work - not any more
      r.resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        r[i] = 0.0;
      }

    } else {
      // active vortons need a radius
      r.resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        r[i] = _in[4*i+3];
      }

      // optional strength in base class
      // need to assign it a vector first!
      Vector<S> new_s;
      new_s.resize(this->n);
      for (size_t i=0; i<this->n; ++i) {
        new_s[i] = _in[4*i+2];
      }
      this->s = std::move(new_s);
    }

    // velocity in base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(this->n);
    }
  }

  const Vector<S>& get_rad() const { return r; }
  Vector<S>&       get_rad()       { return r; }

  // find out the next row index in the BEM after this collection
  // once we start supporting BEM unknowns on points, we'll have to change these
  void set_first_row(const Int _i) { return; }
  const Int get_first_row() const { return 0; }
  const Int get_num_rows()  const { return 0; }
  const Int get_next_row()  const { return 0; }

  void add_new(std::vector<float>& _in) {
    // remember old size and incoming size
    const size_t nold = this->n;

    const size_t nper = (this->E == inert) ? 2 : 4;
    const size_t nnew = _in.size()/nper;
    std::cout << "  adding " << nnew << " particles to collection..." << std::endl;

    // must explicitly call the method in the base class first
    ElementBase<S>::add_new(_in);

    // then do local stuff
    r.resize(nold+nnew);
    if (this->E == inert) {
      for (size_t i=0; i<nnew; ++i) {
        r[nold+i] = 0.0;
      }
    } else {
      for (size_t i=0; i<nnew; ++i) {
        r[nold+i] = _in[4*i+3];
      }
    }

    // save the new untransformed positions if we have a Body pointer
    if (this->B) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*this->ux)[d].resize(nold+nnew);
        for (size_t i=nold; i<nold+nnew; ++i) {
          (*this->ux)[d][i] = this->x[d][i];
        }
      }
    }
  }

  // up-size all arrays to the new size, filling with sane values
  void resize(const size_t _nnew) {
    const size_t currn = this->n;
    //std::cout << "  inside Points::resize with " << currn << " " << _nnew << std::endl;

    // must explicitly call the method in the base class - this sets n
    ElementBase<S>::resize(_nnew);

    if (_nnew == currn) return;

    // radii here
    const size_t thisn = r.size();
    r.resize(_nnew);
    for (size_t i=thisn; i<_nnew; ++i) {
      r[i] = 1.0;
    }
  }

  void zero_vels() {
    // must explicitly call the method in the base class to zero the vels
    ElementBase<S>::zero_vels();
  }

  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    // must explicitly call the method in the base class, too
    ElementBase<S>::finalize_vels(_fs);
  }

  void transform(const double _time) {
    // must explicitly call the method in the base class
    ElementBase<S>::transform(_time);
    // no other specialization required here
  }

  //
  // 1st order Euler advection and stretch
  //
  void move(const double _time, const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt);

    // no specialization needed
    if (this->M == lagrangian and this->E != inert) {
      //std::cout << "  Stretching" << to_string() << " using 1st order" << std::endl;
      S thismax = 0.0;

      for (size_t i=0; i<this->n; ++i) {
        S this_s = (*this->s)[i];

        // compute stretch term
        std::array<S,2> wdu = {0.0};

        // add Cottet SFS

        // update strengths
        (*this->s)[i] = this_s + _dt * wdu[0];

        // check for max strength
        S thisstr = std::abs((*this->s)[i]);
        if (thisstr > thismax) thismax = thisstr;

      }
      if (max_strength < 0.0) {
        max_strength = thismax;
      } else {
        max_strength = 0.05*thismax + 0.95*max_strength;
      }
      //std::cout << "  New max_strength is " << max_strength << std::endl;
    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }

  //
  // 2nd order RK advection and stretch
  //
  void move(const double _time, const double _dt,
            const double _wt1, Points<S> const & _u1,
            const double _wt2, Points<S> const & _u2) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt, _wt1, _u1, _wt2, _u2);

    // must confirm that incoming time derivates include velocity

    // and specialize
    if (this->M == lagrangian and this->E != inert) {
      //std::cout << "  Stretching" << to_string() << " using 2nd order" << std::endl;
      S thismax = 0.0;

      for (size_t i=0; i<this->n; ++i) {

        // set up some convenient temporaries
        S this_s = (*this->s)[i];

        // compute stretch term
        std::array<S,2> wdu1 = {0.0};
        std::array<S,2> wdu2 = {0.0};
        std::array<S,2> wdu  = {0.0};

        wdu[0] = _wt1*wdu1[0] + _wt2*wdu2[0];
        wdu[1] = _wt1*wdu1[1] + _wt2*wdu2[1];

        // add Cottet SFS

        // update strengths
        (*this->s)[i] = this_s + _dt * wdu[0];

        // check for max strength
        S thisstr = std::abs((*this->s)[i]);
        if (thisstr > thismax) thismax = thisstr;
      }

      if (max_strength < 0.0) {
        max_strength = thismax;
      } else {
        max_strength = 0.05*thismax + 0.95*max_strength;
      }

    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }

  // find the new peak strength magnitude
  void update_max_str() {
    S thismax = ElementBase<S>::get_max_str();

    // and slowly update the current saved value
    if (max_strength < 0.0) {
      max_strength = thismax;
    } else {
      max_strength = 0.1*thismax + 0.9*max_strength;
    }
  }

  void add_body_motion(const S _factor, const double _time) {
    // no need to call base class now
    //ElementBase<S>::add_body_motion(_factor);
    // apply a factor times the body motion
    //for (size_t i=0; i<this->get_n(); ++i) {
    //  std::array<S,Dimensions> thisvel = B->get_vel(_time);
    //  // apply the velocity
    //  for (size_t d=0; d<Dimensions; ++d) {
    //    this->u[d][i] += _factor * thisvel[d];
    //  }
    //}
  }

#ifdef USE_GL
  //
  // OpenGL functions
  //

  // helper function to clean up initGL
  void prepare_opengl_buffer(GLuint _prog, GLuint _idx, const GLchar* _name) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo[_idx]);
    const GLint position_attribute = glGetAttribLocation(_prog, _name);
    // Specify how the data for position can be accessed
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    // Enable the attribute
    glEnableVertexAttribArray(position_attribute);
    // and tell it to advance two primitives per point (2 tris per quad)
    glVertexAttribDivisor(position_attribute, 1);
  }

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor,
              float*              _defcolor) {

    //std::cout << "inside Points.initGL" << std::endl;
    std::cout << "inside Points.initGL with E=" << this->E << " and M=" << this->M << std::endl;

    // Use a Vertex Array Object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Create four Vector Buffer Objects that will store the vertex arrays in gpu memory
    glGenBuffers(4, vbo);

    // Allocate space, but don't upload the data from CPU to GPU yet
    for (size_t i=0; i<Dimensions; ++i) {
      glBindBuffer(GL_ARRAY_BUFFER, vbo[i]);
      glBufferData(GL_ARRAY_BUFFER, 0, this->x[i].data(), GL_STATIC_DRAW);
    }

    // here is where we split on element type: active/reactive vs. inert
    if (this->E == inert) {

      // Load and create the blob-drawing shader program
      draw_point_program = create_draw_point_program();

      // Now do the four arrays
      prepare_opengl_buffer(draw_point_program, 0, "px");
      prepare_opengl_buffer(draw_point_program, 1, "py");

      // and for the compute shaders!

      // Get the location of the attributes that enters in the vertex shader
      projmat_attribute_pt = glGetUniformLocation(draw_point_program, "Projection");

      // upload the projection matrix
      glUniformMatrix4fv(projmat_attribute_pt, 1, GL_FALSE, _projmat.data());

      // locate where the colors and color scales go
      def_color_attribute = glGetUniformLocation(draw_point_program, "def_color");

      // send the current values
      glUniform4fv(def_color_attribute, 1, _defcolor);

      // locate where the point radius goes
      unif_rad_attribute = glGetUniformLocation(draw_point_program, "rad");

      // send the current values
      glUniform1f(unif_rad_attribute, (const GLfloat)0.01);

      // and indicate the fragment color output
      glBindFragDataLocation(draw_point_program, 0, "frag_color");

      // Initialize the quad attributes
      std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
      GLuint qvbo;
      glGenBuffers(1, &qvbo);
      glBindBuffer(GL_ARRAY_BUFFER, qvbo);
      glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

      quad_attribute_pt = glGetAttribLocation(draw_point_program, "quad_attr");
      glVertexAttribPointer(quad_attribute_pt, 2, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray(quad_attribute_pt);

    } else { // this->E is active or reactive

      glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
      glBufferData(GL_ARRAY_BUFFER, 0, r.data(), GL_STATIC_DRAW);
      if (this->s) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
        glBufferData(GL_ARRAY_BUFFER, 0, (*this->s).data(), GL_STATIC_DRAW);
      }

      // Load and create the blob-drawing shader program
      draw_blob_program = create_draw_blob_program();

      // Now do the four arrays
      prepare_opengl_buffer(draw_blob_program, 0, "px");
      prepare_opengl_buffer(draw_blob_program, 1, "py");
      prepare_opengl_buffer(draw_blob_program, 2, "r");
      prepare_opengl_buffer(draw_blob_program, 3, "sx");

      // and for the compute shaders!

      // Get the location of the attributes that enters in the vertex shader
      projmat_attribute_bl = glGetUniformLocation(draw_blob_program, "Projection");

      // upload the projection matrix
      glUniformMatrix4fv(projmat_attribute_bl, 1, GL_FALSE, _projmat.data());

      // locate where the colors and color scales go
      pos_color_attribute = glGetUniformLocation(draw_blob_program, "pos_color");
      neg_color_attribute = glGetUniformLocation(draw_blob_program, "neg_color");
      str_scale_attribute = glGetUniformLocation(draw_blob_program, "str_scale");

      // send the current values
      glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_poscolor);
      glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_negcolor);
      glUniform1f (str_scale_attribute, (const GLfloat)1.0);
      //std::cout << "init pos color as " << _poscolor[0] << " " << _poscolor[1] << " " << _poscolor[2] << " " << _poscolor[3] << std::endl;

      // and indicate the fragment color output
      glBindFragDataLocation(draw_blob_program, 0, "frag_color");

      // Initialize the quad attributes
      std::vector<float> quadverts = {-1,-1, 1,-1, 1,1, -1,1};
      GLuint qvbo;
      glGenBuffers(1, &qvbo);
      glBindBuffer(GL_ARRAY_BUFFER, qvbo);
      glBufferData(GL_ARRAY_BUFFER, quadverts.size()*sizeof(float), quadverts.data(), GL_STATIC_DRAW);

      quad_attribute_bl = glGetAttribLocation(draw_blob_program, "quad_attr");
      glVertexAttribPointer(quad_attribute_bl, 2, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray(quad_attribute_bl);

    } // end this->E is active or reactive

    glBindVertexArray(0);
  }

  // this gets done every time we change the size of the positions array
  void updateGL() {
    //std::cout << "inside Points.updateGL" << std::endl;

    // has this been init'd yet?
    if (glIsVertexArray(vao) == GL_FALSE) return;

    const size_t vlen = this->x[0].size()*sizeof(S);
    if (vlen > 0) {
      glBindVertexArray(vao);

      // Indicate and upload the data from CPU to GPU
      for (size_t i=0; i<Dimensions; ++i) {
        // the positions
        glBindBuffer(GL_ARRAY_BUFFER, vbo[i]);
        glBufferData(GL_ARRAY_BUFFER, vlen, this->x[i].data(), GL_DYNAMIC_DRAW);
      }

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {

        // just don't upload strengths or radii

      } else { // this->E is active or reactive

        // the strengths
        if (this->s) {
          glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
          glBufferData(GL_ARRAY_BUFFER, vlen, (*this->s).data(), GL_DYNAMIC_DRAW);
        }
        // the radii
        glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
        glBufferData(GL_ARRAY_BUFFER, vlen, r.data(), GL_DYNAMIC_DRAW);
      }

      glBindVertexArray(0);

      // must tell draw call how many elements are there
      num_uploaded = this->x[0].size();
    }
  }

  // OpenGL3 stuff to display points, called once per frame
  void drawGL(std::vector<float>& _projmat,
              RenderParams&       _rparams) {

    //std::cout << "inside Points.drawGL" << std::endl;

    // has this been init'd yet?
    if (glIsVertexArray(vao) == GL_FALSE) {
      initGL(_projmat, _rparams.pos_circ_color,
                       _rparams.neg_circ_color,
                       _rparams.default_color);
      updateGL();
    }

    if (num_uploaded > 0) {
      glBindVertexArray(vao);

      // get blending ready
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {

        // draw as small dots
        glUseProgram(draw_point_program);

        glEnableVertexAttribArray(quad_attribute_pt);

        // upload the current uniforms
        glUniformMatrix4fv(projmat_attribute_pt, 1, GL_FALSE, _projmat.data());
        glUniform4fv(def_color_attribute, 1, (const GLfloat *)_rparams.default_color);
        glUniform1f (unif_rad_attribute, (const GLfloat)(2.5f*_rparams.tracer_size));

        // the one draw call here
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, num_uploaded);

      } else { // this->E is active or reactive

        // draw as colored clouds
        glUseProgram(draw_blob_program);

        glEnableVertexAttribArray(quad_attribute_bl);

        // upload the current projection matrix
        glUniformMatrix4fv(projmat_attribute_bl, 1, GL_FALSE, _projmat.data());

        // upload the current color values
        glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_rparams.pos_circ_color);
        glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_rparams.neg_circ_color);
        glUniform1f (str_scale_attribute, (const GLfloat)(0.4f/max_strength));

        // the one draw call here
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, num_uploaded);
      }

      // return state
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glBindVertexArray(0);
    }
  }
#endif

  std::string to_string() const {
    std::string retstr = ElementBase<S>::to_string() + " Points";
    return retstr;
  }

protected:
  // additional state vector
  Vector<S> r;					// thickness/radius
  // in 3D, this is where elong and ug would be

private:
#ifdef USE_GL
  // OpenGL stuff
  GLuint vao, vbo[4];
  GLuint draw_blob_program, draw_point_program;
  GLsizei num_uploaded;
  GLint projmat_attribute_bl, projmat_attribute_pt, quad_attribute_bl, quad_attribute_pt;
  GLint def_color_attribute, pos_color_attribute, neg_color_attribute, str_scale_attribute;
  GLint unif_rad_attribute;
#endif
  float max_strength;
};

