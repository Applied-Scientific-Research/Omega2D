/*
 * Surfaces.h - Specialized class for surfaces in 2D
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega2D.h"
#include "VectorHelper.h"
#include "ElementBase.h"

#ifdef USE_GL
#include "OglHelper.h"
#include "ShaderHelper.h"
#include "glad.h"
#endif

#include <iostream>
#include <vector>
#include <array>
//#include <memory>
//#include <optional>
//#include <random>
#include <cassert>
#define _USE_MATH_DEFINES
//#include <cmath>
//#include <algorithm>


// 1-D elements
template <class S>
class Surfaces: public ElementBase<S> {
public:
  // constructor - accepts vector of vectors of (x,y,s) pairs
  //               makes one closed body for each outer vector
  //               each inner vector must have even number of floats
  //               and first pair must equal last pair to close the shape
  //               last parameter (s) is either fixed strength or boundary
  //               condition for next segment
  Surfaces(const std::vector<S>&   _x,
           const std::vector<Int>& _idx,
           const std::vector<S>&   _val,
           const elem_t _e, const move_t _m)
    : ElementBase<S>(0, _e, _m),
      max_strength(-1.0) {

    // make sure input arrays are correctly-sized
    assert(_x.size() % 2 == 0);
    assert(_idx.size() % 2 == 0);
    assert(_idx.size()/2 == _val.size());
    const size_t nnodes = _x.size() / 2;
    const size_t nsurfs = _idx.size() / 2;

    std::cout << "  new collection with " << nsurfs << " surface panels..." << std::endl;

    // pull out the node locations
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(nnodes);
      for (size_t i=0; i<nnodes; ++i) {
        this->x[d][i] = _x[2*i+d];
      }
    }

    // copy over the node indices (with a possible type change)
    bool idx_are_all_good = true;
    idx.resize(_idx.size());
    for (size_t i=0; i<2*nsurfs; ++i) {
      // make sure it exists in the nodes array
      if (_idx[i] >= nnodes) idx_are_all_good = false;
      idx[i] = _idx[i];
    }
    assert(idx_are_all_good);

    // now, depending on the element type, put the value somewhere
    if (this->E == active) {
      // value is a fixed strength for the segment
      Vector<S> new_s(_val.size());
      std::copy(_val.begin(), _val.end(), new_s.begin());
      //new_s.resize(nsurfs);
      this->s = std::move(new_s);

    } else if (this->E == reactive) {
      // value is a boundary condition
      bc.resize(_val.size());
      std::copy(_val.begin(), _val.end(), bc.begin());
      //bc.resize(nsurfs);
      //for (size_t i=0; i<nsurfs; ++i) {
      //  bc[i] = _val[i];
      //}
      // we still need strengths
      Vector<S> new_s(_val.size());
      this->s = std::move(new_s);

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnodes);
    }

    // debug print
    if (false) {
      std::cout << "Nodes" << std::endl;
      for (size_t i=0; i<nnodes; ++i) {
        std::cout << "  " << i << " " << this->x[0][i] << " " << this->x[1][i] << std::endl;
      }
      std::cout << "Segments" << std::endl;
      for (size_t i=0; i<nsurfs; ++i) {
        std::cout << "  " << i << " " << idx[2*i] << " " << idx[2*i+1] << std::endl;
      }
    }

    // need to reset the base class n
    this->n = nnodes;
  }

  size_t get_npanels() const { return idx.size()/2; }

  // callers should never have to change this array
  const std::vector<Int>& get_idx() const { return idx; }
  const std::vector<S>&   get_bcs() const { return bc; }


  // add more nodes and panels to this collection
  void add_new(const std::vector<S>&   _x,
               const std::vector<Int>& _idx,
               const std::vector<S>&   _val) {

    // remember old sizes of nodes and element arrays
    const size_t nnold = this->n;
    const size_t neold = get_npanels();

    // make sure input arrays are correctly-sized
    assert(_x.size() % 2 == 0);
    assert(_idx.size() % 2 == 0);
    assert(_idx.size()/2 == _val.size());
    const size_t nnodes = _x.size() / 2;
    const size_t nsurfs = _idx.size() / 2;

    std::cout << "  adding " << nsurfs << " new surface panels to collection..." << std::endl;

    // DON'T call the method in the base class, because we do things differently here
    //ElementBase<S>::add_new(_in);

    // pull out the node locations
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(nnold+nnodes);
      for (size_t i=0; i<nnodes; ++i) {
        this->x[d][nnold+i] = _x[2*i+d];
      }
    }

    // copy over the node indices, taking care to offset into the new array
    bool idx_are_all_good = true;
    idx.resize(2*neold + _idx.size());
    for (size_t i=0; i<2*nsurfs; ++i) {
      // make sure it exists in the nodes array
      if (_idx[i] >= nnold+nnodes) idx_are_all_good = false;
      idx[2*neold+i] = nnold + _idx[i];
    }
    assert(idx_are_all_good);

    // now, depending on the element type, put the value somewhere
    if (this->E == active) {
      // value is a fixed strength for the element
      //*(this->s).resize(neold+nsurfs);
      //(*this->s)[i]
      this->s->reserve(neold+nsurfs); 

    } else if (this->E == reactive) {
      // value is a boundary condition
      bc.reserve(neold+nsurfs); 
      bc.insert(bc.end(), _val.begin(), _val.end());
      // and we still need strengths
      this->s->reserve(neold+nsurfs); 

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnold+nnodes);
    }

    // debug print
    if (false) {
      std::cout << "Nodes" << std::endl;
      for (size_t i=0; i<nnold+nnodes; ++i) {
        std::cout << "  " << i << " " << this->x[0][i] << " " << this->x[1][i] << std::endl;
      }
      std::cout << "Segments" << std::endl;
      for (size_t i=0; i<neold+nsurfs; ++i) {
        std::cout << "  " << i << " " << idx[2*i] << " " << idx[2*i+1] << std::endl;
      }
    }

    // need to reset the base class n
    this->n += nnodes;
  }

/*
  // up-size all arrays to the new size, filling with sane values
  void resize(const size_t _nnew) {
    const size_t currn = this->n;
    //std::cout << "  inside Surfaces::resize with " << currn << " " << _nnew << std::endl;

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

  //
  // 1st order Euler advection and stretch
  //
  void move(const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_dt);

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
        max_strength = 0.1*thismax + 0.9*max_strength;
      }
      //std::cout << "  New max_strength is " << max_strength << std::endl;
    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }
*/

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
  }

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor,
              float*              _defcolor) {

    //std::cout << "inside Surfaces.initGL" << std::endl;
    std::cout << "inside Surfaces.initGL with E=" << this->E << " and M=" << this->M << std::endl;

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

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[2]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, idx.data(), GL_STATIC_DRAW);

    if (this->s) {
      glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
      glBufferData(GL_ARRAY_BUFFER, 0, (*this->s).data(), GL_STATIC_DRAW);
    }

    // Load and create the blob-drawing shader program
    draw_surface_line_prog = create_draw_surface_line_prog();

    // Now do the four arrays
    prepare_opengl_buffer(draw_surface_line_prog, 0, "px");
    prepare_opengl_buffer(draw_surface_line_prog, 1, "py");
    prepare_opengl_buffer(draw_surface_line_prog, 2, "rawstr");

    // and for the compute shaders! (later)

    // Get the location of the attributes that enters in the vertex shader
    projmat_attribute = glGetUniformLocation(draw_surface_line_prog, "Projection");

    // upload the projection matrix
    glUniformMatrix4fv(projmat_attribute, 1, GL_FALSE, _projmat.data());

    // locate where the colors and color scales go
    pos_color_attribute = glGetUniformLocation(draw_surface_line_prog, "pos_color");
    neg_color_attribute = glGetUniformLocation(draw_surface_line_prog, "neg_color");
    def_color_attribute = glGetUniformLocation(draw_surface_line_prog, "def_color");
    str_scale_attribute = glGetUniformLocation(draw_surface_line_prog, "str_scale");

    // send the current values
    glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_poscolor);
    glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_negcolor);
    glUniform4fv(def_color_attribute, 1, (const GLfloat *)_defcolor);
    glUniform1f (str_scale_attribute, (const GLfloat)1.0);
    //std::cout << "init pos color as " << _poscolor[0] << " " << _poscolor[1] << " " << _poscolor[2] << " " << _poscolor[3] << std::endl;

    // and indicate the fragment color output
    glBindFragDataLocation(draw_surface_line_prog, 0, "frag_color");

    glBindVertexArray(0);
  }

  // this gets done every time we change the size of the index array
  void updateGL() {
    //std::cout << "inside Surfaces.updateGL" << std::endl;

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

      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[2]);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Int)*idx.size(), idx.data(), GL_DYNAMIC_DRAW);

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {
        // just don't upload strengths
      } else { // this->E is active or reactive
        // the strengths
        if (this->s) {
          glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
          glBufferData(GL_ARRAY_BUFFER, vlen, (*this->s).data(), GL_DYNAMIC_DRAW);
        }
      }

      glBindVertexArray(0);

      // must tell draw call how many elements are there
      num_uploaded = idx.size();
    }
  }

  // OpenGL3 stuff to display points, called once per frame
  void drawGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor,
              float*              _defcolor,
              float               _tracersize) {

    //std::cout << "inside Surfaces.drawGL" << std::endl;

    // has this been init'd yet?
    if (glIsVertexArray(vao) == GL_FALSE) {
      initGL(_projmat, _poscolor, _negcolor, _defcolor);
      updateGL();
    }

    if (num_uploaded > 0) {
      glBindVertexArray(vao);

      // get blending ready
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      glLineWidth(2.0);

      // here is where we split on element type: active/reactive vs. inert
      //if (this->E == inert) {
      //} else { // this->E is active or reactive

      // draw as lines
      glUseProgram(draw_surface_line_prog);

      // upload the current projection matrix
      glUniformMatrix4fv(projmat_attribute, 1, GL_FALSE, _projmat.data());

      // upload the current color values
      glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_poscolor);
      glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_negcolor);
      glUniform4fv(def_color_attribute, 1, (const GLfloat *)_defcolor);
      glUniform1f (str_scale_attribute, (const GLfloat)max_strength);

      // the one draw call here
      glDrawElements(GL_LINES, num_uploaded, get_gl_type<Int>, 0);

      // return state
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glBindVertexArray(0);
    }
  }
#endif

  std::string to_string() const {
    std::string retstr = ElementBase<S>::to_string() + " Panels";
    return retstr;
  }

protected:
  // ElementBase.h has x, s, u
  std::vector<Int>	idx;	// indexes into the x array
  Vector<S>		bc;	// boundary condition for the elements

  //std::vector<std::pair<Int,Int>> body_idx;	// n, offset of rows in the BEM?

private:
#ifdef USE_GL
  // OpenGL stuff
  GLuint vao, vbo[4];
  GLuint draw_surface_line_prog;
  GLsizei num_uploaded;
  GLint projmat_attribute;
  GLint def_color_attribute, pos_color_attribute, neg_color_attribute, str_scale_attribute;
#endif
  float max_strength;
};

