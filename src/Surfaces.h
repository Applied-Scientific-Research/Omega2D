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
#include "RenderParams.h"
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
           const elem_t _e,
           const move_t _m,
           std::shared_ptr<Body> _bp)
    : ElementBase<S>(0, _e, _m, _bp),
      area(-1.0),
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

    // save them as untransformed if we are given a Body pointer
    if (_bp) {
      this->ux = this->x;
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

    // find geometric center
    if (this->M == bodybound) {
      set_geom_center();
    }
  }

  size_t get_npanels() const { return idx.size()/2; }
  const S get_area() const { return area; }
  const std::array<S,2> get_geom_center() const { return tc; }

  // callers should never have to change this array
  const std::vector<Int>& get_idx() const { return idx; }
  const Vector<S>&        get_bcs() const { return bc; }

  // find out the next row index in the BEM after this collection
  void set_first_row(const Int _i) { istart = _i; }
  const Int get_first_row() const { return istart; }
  const Int get_num_rows()  const { return bc.size() + (this->B ? 1 : 0); }
  const Int get_next_row()  const { return istart+get_num_rows(); }

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

    // save them as untransformed if we have a Body pointer
    if (this->B) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*this->ux)[d].resize(nnold+nnodes);
        for (size_t i=nnold; i<nnold+nnodes; ++i) {
          (*this->ux)[d][i] = this->x[d][i];
        }
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
      this->s->resize(neold+nsurfs); 

    } else if (this->E == reactive) {
      // value is a boundary condition
      bc.reserve(neold+nsurfs); 
      bc.insert(bc.end(), _val.begin(), _val.end());
      // and we still need strengths
      this->s->resize(neold+nsurfs); 

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

    // re-find geometric center
    if (this->M == bodybound) {
      set_geom_center();
    }
  }

  void add_body_motion(const S _factor, const double _time) {
    // no need to call base class now
    //ElementBase<S>::add_body_motion(_factor);

    // make sure we've calculated transformed center (we do this when we do area)
    assert(area > 0.0);
    // and we trust that we've transformed utc to tc

    // apply a factor times the body motion
    for (size_t i=0; i<this->get_n(); ++i) {

      // apply the translational velocity
      std::array<double,Dimensions> thisvel = this->B->get_vel(_time);
      for (size_t d=0; d<Dimensions; ++d) {
        this->u[d][i] += _factor * (float)thisvel[d];
      }

      // now compute the rotational velocity with respect to the geometric center
      double thisrotvel = this->B->get_rotvel(_time);
      // center of this panel
      Int id0 = idx[2*i];
      Int id1 = idx[2*i+1];
      // panel center
      const S xc = 0.5 * (this->x[0][id1] + this->x[0][id0]);
      const S yc = 0.5 * (this->x[1][id1] + this->x[1][id0]);
      // add rotational velocity
      this->u[0][i] -= _factor * (float)thisrotvel * (yc - tc[1]);
      this->u[1][i] += _factor * (float)thisrotvel * (xc - tc[0]);
    }
  }

  // augment the strengths with a value equal to that which accounts for
  //   the solid-body rotation of the object
  // NOTE: this needs to provide both the vortex AND source strengths!
  // AND we don't have the time - assume bodies have been transformed
  void add_rot_strengths(const S _constfac, const S _rotfactor) {

    // make sure we've calculated transformed center (we do this when we do area)
    assert(area > 0.0);
    // and we trust that we've transformed utc to tc

    // if no rotation, strengths, or no parent Body, then no problem!
    if (not this->B) return;
    if (not this->s) return;
    const S rotvel = (S)this->B->get_rotvel();
    //if (std::abs(rotvel) < std::numeric_limits<float>::epsilon()) return;

    // have we made ss yet?
    if (not ss) {
      // value is a fixed strength for the segment
      Vector<S> new_ss(this->s->size());
      ss = std::move(new_ss);
    }

    // what is the actual factor that we will add?
    const S factor = _constfac + rotvel*_rotfactor;

    // still here? let's do it. use the untransformed coordinates
    for (size_t i=0; i<get_npanels(); i++) {
      const size_t j   = idx[2*i];
      const size_t jp1 = idx[2*i+1];
      // vector from object geometric center to panel center
      const S dx = 0.5 * ((*this->ux)[0][j] + (*this->ux)[0][jp1]) - utc[0];
      const S dy = 0.5 * ((*this->ux)[1][j] + (*this->ux)[1][jp1]) - utc[1];
      // velocity of the panel center
      const S ui = -factor * dy;
      const S vi =  factor * dx;

      // panel tangential vector, fluid to the left, body to the right
      S panelx = (*this->ux)[0][jp1] - (*this->ux)[0][j];
      S panely = (*this->ux)[1][jp1] - (*this->ux)[1][j];
      const S panell = 1.0 / std::sqrt(panelx*panelx + panely*panely);
      panelx *= panell;
      panely *= panell;

      // the vortex strength - we ADD to the existing
      (*this->s)[i] += -1.0 * (ui*panelx + vi*panely);

      // the source strength
      (*ss)[i] += -1.0 * (ui*panely - vi*panelx);

      // debug print
      if (_rotfactor > 0.0 and false) {
        std::cout << "  panel " << i << " at " << dx << " " << dy
                  << " adds to vortex str " << (-1.0 * (ui*panelx + vi*panely))
                  << " and source str " << (-1.0 * (ui*panely - vi*panelx)) << std::endl;
      }
    }
  }

  // calculate the geometric center of all geometry in this object
  void set_geom_center() {

    // we must have an attached body and a set of untransformed coordinates
    assert(this->B);
    assert(this->ux);

    std::cout << "  inside Surfaces::set_geom_center with " << get_npanels() << " panels" << std::endl;
    // (*this->ux)[d][i]

    // iterate over panels, accumulating area and CM
    S asum = 0.0;
    S xsum = 0.0;
    S ysum = 0.0;
    for (size_t i=0; i<get_npanels(); i++) {
      const size_t j   = idx[2*i];
      const size_t jp1 = idx[2*i+1];
      // assume a triangle from 0,0 to two ends of each panel
      const S xc = (0.0 + (*this->ux)[0][j] + (*this->ux)[0][jp1]) / 3.0;
      const S yc = (0.0 + (*this->ux)[1][j] + (*this->ux)[1][jp1]) / 3.0;
      const S panelx = (*this->ux)[0][jp1] - (*this->ux)[0][j];
      const S panely = (*this->ux)[1][jp1] - (*this->ux)[1][j];
      // and the side lengths
      const S a = std::sqrt(std::pow((*this->ux)[0][j],2)+std::pow((*this->ux)[1][j],2));
      const S b = std::sqrt(std::pow(panelx,2)+std::pow(panely,2));
      const S c = std::sqrt(std::pow((*this->ux)[0][jp1],2)+std::pow((*this->ux)[1][jp1],2));
      //std::cout << "  panel " << i << " has side lens " << a << " " << b << " " << c << std::endl;
      // Heron's formula for the area
      const S hs = 0.5*(a+b+c);
      S thisarea = std::sqrt(hs*(hs-a)*(hs-b)*(hs-c));
      // negate area if the winding is backwards
      if ((*this->ux)[1][j]*panelx - (*this->ux)[0][j]*panely < 0.0) thisarea = -thisarea;
      // add this to the running sums
      //std::cout << "    and area " << thisarea << " and center " << xc << " " << yc << std::endl;
      asum += thisarea;
      xsum += xc*thisarea;
      ysum += yc*thisarea;
    }
    area = asum;
    utc[0] = xsum/area;
    utc[1] = ysum/area;

    std::cout << "    geom center is " << utc[0] << " " << utc[1] << " and area is " << area << std::endl;
  }

  // when transforming a body-bound object to a new time, we must also transform the geometric center
  void transform(const double _time) {
    // must explicitly call the method in the base class
    ElementBase<S>::transform(_time);

    // prepare for the transform
    std::array<double,Dimensions> thispos = this->B->get_pos();
    const double theta = this->B->get_orient();
    const S st = std::sin(theta);
    const S ct = std::cos(theta);

    // transform the utc to tc here
    tc[0] = (S)thispos[0] + utc[0]*ct - utc[1]*st;
    tc[1] = (S)thispos[1] + utc[0]*st + utc[1]*ct;
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
        max_strength = 0.1*thismax + 0.9*max_strength;
      }
      //std::cout << "  New max_strength is " << max_strength << std::endl;
    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }
*/

  // return a particle version of the panels (useful during Diffusion)
  std::vector<S> represent_as_particles(const S _offset, const S _vdelta) {

    // how many panels?
    const size_t num_pts = get_npanels();

    // init the output vector
    std::vector<S> px(num_pts*4);

    // outside is to the left walking from one point to the next
    // so go CW around the circle starting at theta=0 (+x axis)
    S oopanlen, along[2];

    for (size_t i=0; i<num_pts; i++) {
      Int id0 = idx[2*i];
      Int id1 = idx[2*i+1];
      // start at center of panel
      px[4*i+0] = 0.5 * (this->x[0][id1] + this->x[0][id0]);
      px[4*i+1] = 0.5 * (this->x[1][id1] + this->x[1][id0]);
      // push out a fixed distance
      along[0] = this->x[0][id1] - this->x[0][id0];
      along[1] = this->x[1][id1] - this->x[1][id0];
      // one over the panel length is useful
      oopanlen = 1.0 / std::sqrt(along[0]*along[0] + along[1]*along[1]);
      // this assumes properly resolved, vdelta and dt
      px[4*i+0] += _offset * -along[1] * oopanlen;
      px[4*i+1] += _offset *  along[0] * oopanlen;
      // complete the element with a strength and radius
      px[4*i+2] = (*this->s)[i] / oopanlen;
      px[4*i+3] = _vdelta;
      //std::cout << "  new part is " << px[4*i+0] << " " << px[4*i+1] << " " << px[4*i+2] << " " << px[4*i+3] << std::endl;
    }

    return px;
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

  // add and return the total circulation of all elements
  S get_total_circ(const double _time) {
    // do not call the parent
    if (this->B) {
      // we're attached to a body - great! what's the rotation rate?
      return 2.0 * area * (S)this->B->get_rotvel(_time);
    } else {
      // we are fixed, thus not rotating
      return 0.0;
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
              RenderParams&       _rparams) {

    //std::cout << "inside Surfaces.drawGL" << std::endl;

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

      glLineWidth(2.0);

      // here is where we split on element type: active/reactive vs. inert
      //if (this->E == inert) {
      //} else { // this->E is active or reactive

      // draw as lines
      glUseProgram(draw_surface_line_prog);

      // upload the current projection matrix
      glUniformMatrix4fv(projmat_attribute, 1, GL_FALSE, _projmat.data());

      // upload the current color values
      glUniform4fv(pos_color_attribute, 1, (const GLfloat *)_rparams.pos_circ_color);
      glUniform4fv(neg_color_attribute, 1, (const GLfloat *)_rparams.neg_circ_color);
      glUniform4fv(def_color_attribute, 1, (const GLfloat *)_rparams.default_color);
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
  Int istart;	// index of first entry in RHS vector and A matrix

  // the source strengths per unit length which represent the velocity
  //   influence of the volume vorticity of the parent body
  std::optional<Vector<S>> ss;

  S area;			// area of the body - for augmented BEM solution
  std::array<S,2> utc;		// untransformed geometric center
  std::array<S,2>  tc;		// transformed geometric center

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

