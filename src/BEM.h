/*
 * BEM.h - Library code for a 2D vortex boundary element solver
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Kernels.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>		// for BiCGSTAB and GMRES
#include <unsupported/Eigen/src/IterativeSolvers/GMRES.h>	// for GMRES

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <chrono>
#define _USE_MATH_DEFINES // Required by MSVC to define M_PI,etc. in <cmath>
#include <cmath>
#include <vector>

//
// Class to hold BEM parameters and temporaries
//
// templatized on 'S'torage type and 'I'ndex type
//
template <class S, class I>
class BEM {
public:
  BEM() : A_is_current(false), solver_initialized(false) {};

  bool is_A_current() { return A_is_current; }
  void panels_changed() { A_is_current = false; solver_initialized = false; }
  void assemble_influence_matrix(const std::vector<S>&, const std::vector<I>&);
  void set_rhs(std::vector<S>&);
  void solve();

  std::vector<S> getRhs();
  std::vector<S> getStrengths();

protected:

private:
  // the actual matrix equation
  Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<S, Eigen::Dynamic, 1> b;
  Eigen::Matrix<S, Eigen::Dynamic, 1> strengths;

  // is the A matrix current?
  bool A_is_current;
  bool solver_initialized;
};

//
// convert Eigen matrix to c++ vector
//
template <class S, class I>
std::vector<S> BEM<S,I>::getRhs() {
  std::vector<S> retval(b.size());
  retval.assign(b.data(), b.data()+b.size());
  return retval;
}

template <class S, class I>
std::vector<S> BEM<S,I>::getStrengths() {
  std::vector<S> retval(strengths.size());
  retval.assign(strengths.data(), strengths.data()+strengths.size());
  return retval;
}

//
// assemble the influence matrix
//
// galerkin or colocation approach?
//
template <class S, class I>
void BEM<S,I>::assemble_influence_matrix( const std::vector<S>& x, const std::vector<I>& idx) {

  // check for existence and size of A
  const I thisn = idx.size() / 2;

  // allocate space
  A.resize(thisn, thisn);
  //A.setRandom(num_rows,num_rows);
  //A.setIdentity(num_rows,num_rows);

  // calculate influences
  #pragma omp parallel for
  for (size_t j=0; j<thisn; j++) {

    for (size_t i=0; i<thisn; i++) {
      const I first  = idx[2*i];
      const I second = idx[2*i+1];
      const S x0 = x[2*first];
      const S y0 = x[2*first+1];
      const S x1 = x[2*second];
      const S y1 = x[2*second+1];

      // collocation point for panel i
      const S xi = 0.5 * (x1 + x0);
      const S yi = 0.5 * (y1 + y0);

      // influence of vortex panel j with unit circulation on center of panel i
      auto vel = vortex_panel_affects_point<S,S>(x[2*idx[2*j]],   x[2*idx[2*j]+1],
                                                 x[2*idx[2*j+1]], x[2*idx[2*j+1]+1],
                                                 1.0, xi, yi);

      // target panel vector
      const S panelx = x1 - x0;
      const S panely = y1 - y0;
      const S panell = std::sqrt(panelx*panelx + panely*panely);

      // dot product with tangent vector, applying normalization here
      A(i,j) = (vel[0]*panelx + vel[1]*panely) / panell;

      // (un-normalized) surface normal pointing into fluid
      //S normx = -panely;
      //S normy = panelx;

      // dot product with normal vector, applying normalization here
      //A(i,j) = (vel[0]*normx + vel[1]*normy) / panell;
    }

    // special case: self-influence
    A(j,j) = M_PI;
    //A(j,j) = 0.0;
  }

  //std::cout << "Here is the matrix A^T:\n" << A.topLeftCorner(6,6) << std::endl;
  //std::cout << "Here is the right hand side b:\n" << b.transpose() << std::endl;

  // scale all influences by constant
  A *= 1.0 / (2.0 * M_PI);

  A_is_current = true;
}

//
// set the rhs vector from a set of input velocities
// trying to make the input "const" is asking for trouble!
//
template <class S, class I>
void BEM<S,I>::set_rhs(std::vector<S>& _b) {
  b.resize(_b.size());
  b = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1>>(_b.data(), _b.size());
}

//
// Find the change in strength that would occur over one dt
//
template <class S, class I>
void BEM<S,I>::solve() {

  bool verbose = false;

  // initialize timers
  std::chrono::system_clock::time_point start, end;
  std::chrono::duration<double> elapsed_seconds;

  // ensure that the solution vector is the right size
  strengths.resizeLike(b);

  //std::cout << "A is " << A.size() << std::endl;
  //std::cout << "b is " << b.size() << std::endl;
  //std::cout << "x is " << strengths.size() << std::endl;

  // construction of A matrix takes:
  // bem2d -circle -d=0.01  in 13m cycles
  // bem2d -circle -d=0.001 in 1.1b cycles

  // trying various dense solvers
  start = std::chrono::system_clock::now();

  // n=100 in 2.7 mil cycles
  // n=1000 in 1060 mil cycles
  // bem2d -circle -d=0.01  in 21m cycles with 4.83514e-07 L2
  //strengths = A.fullPivLu().solve(b);

  // bem2d -circle -d=0.01  in 14m cycles with 5.79192e-07 L2
  //strengths = A.colPivHouseholderQr().solve(b);

  // bem2d -circle -d=0.01  in 3.6m cycles with 4.97579e-07 L2 - WOW!
  // bem2d -circle -d=0.001  in 1.1b cycles with 1.44188e-06 L2 - WOW! (N=3141)
  //strengths = A.partialPivLu().solve(b);

  // try an iterative solver!
  // bem2d -circle -d=0.01  in 0.4m cycles with 4.07142e-07
  // bem2d -circle -d=0.001  in 33m cycles with 1.20702e-06 L2 - WOW! (N=3141)
  // note that square takes 6 iterations (circle takes 1) and 490m cycles
  //Eigen::BiCGSTAB<Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> > solver;

  static Eigen::GMRES<Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> > solver(A);

  if (not solver_initialized) {

    // if A changes, we need to re-run this
    start = std::chrono::system_clock::now();
    solver.compute(A);
    end = std::chrono::system_clock::now();

    elapsed_seconds = end-start;
    printf("    solver.compute:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());

    solver_initialized = true;
  }

  // note that BiCGSTAB accepts a preconditioner as a template arg
  // note that solveWithGuess() can seed the solution with last step's solution!
  start = std::chrono::system_clock::now();
  strengths = solver.solve(b);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  printf("    solver.solve:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());

  if (verbose) printf("    num iterations:     %d\n", (uint32_t)solver.iterations());
  if (verbose) printf("    estimated error: %g\n", solver.error());

  start = std::chrono::system_clock::now();
  double relative_error = (A*strengths - b).norm() / b.norm(); // norm() is L2 norm
  if (verbose) printf("    L2 norm of error is %g\n", relative_error);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  if (verbose) printf("    solver.error:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());

/*
  // make sure all vector sizes are identical
  assert(x.size()==y.size());
  assert(x.size()==r.size());
  assert(x.size()==newr.size());
  assert(x.size()==s.size());
  assert(x.size()==ds.size());
  int32_t n = x.size();
*/
}

