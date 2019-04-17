/*
 * BEM.h - Library code for a 2D vortex boundary element solver
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "VectorHelper.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>		// for BiCGSTAB and GMRES
#include <unsupported/Eigen/src/IterativeSolvers/GMRES.h>	// for GMRES

#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <chrono>
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
  void just_made_A() { A_is_current = true; }
  void panels_changed() { A_is_current = false; solver_initialized = false; }
  void reset();
  void set_block(const size_t, const size_t, const size_t, const size_t, const Vector<S>&);
  void set_rhs(std::vector<S>&);
  void set_rhs(const size_t, const size_t, std::vector<S>&);
  void solve();

  std::vector<S> getRhs();
  std::vector<S> getStrengths();
  Vector<S> get_str(const size_t, const size_t);

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

// remove any memory and reset flags
template <class S, class I>
void BEM<S,I>::reset() {
  A_is_current = false;
  solver_initialized = false;
  A.resize(1,1);
  b.resize(1);
  strengths.resize(1);
}

//
// Convert Eigen matrix to c++ vector
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

template <class S, class I>
Vector<S> BEM<S,I>::get_str(const size_t cstart, const size_t ncols) {
  assert(cstart >= 0 && "Start index is negative");
  assert(ncols >= 0 && "Column count is negative");
  assert(cstart+ncols <= (size_t)strengths.size() && "Asking for more values than are present in array");

  Vector<S> retval(ncols);
  retval.assign(strengths.data()+cstart, strengths.data()+cstart+ncols);
  return retval;
}

//
// Set a block in the A matrix from the given vector of coefficients
//
template <class S, class I>
void BEM<S,I>::set_block(const size_t rstart, const size_t nrows,
                         const size_t cstart, const size_t ncols,
                         const Vector<S>& _in) {

  //std::cout << "    putting data into A matrix at " << rstart << ":" << (rstart+nrows) << " "
  //                                                  << cstart << ":" << (cstart+ncols) << std::endl;

  // allocate space
  const size_t new_rows = std::max((size_t)(A.rows()), (size_t)(rstart+nrows));
  const size_t new_cols = std::max((size_t)(A.cols()), (size_t)(cstart+ncols));
  //std::cout << "    resizing A to " << new_rows << " rows and " << new_cols << " cols" << std::endl;
  A.conservativeResize(new_rows, new_cols);

  size_t iptr = 0;
  for (size_t j=0; j<ncols; ++j) {
    for (size_t i=0; i<nrows; ++i) {
      A(i+rstart,j+cstart) = _in[iptr++];
    }
  }
}

//
// Set the rhs vector from a set of input velocities
// trying to make the input "const" is asking for trouble!
//
template <class S, class I>
void BEM<S,I>::set_rhs(std::vector<S>& _b) {
  b.resize(_b.size());
  b = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1>>(_b.data(), _b.size());
}

template <class S, class I>
void BEM<S,I>::set_rhs(const size_t cstart, const size_t ncols, std::vector<S>& _b) {
  // use cstart and ncols to put the _b vector in a specific place
  std::cout << "    putting data into rhs vector from " << cstart << " to " << (cstart+ncols) << std::endl;
  assert(_b.size() == ncols && "Input array size does not match");
  const size_t new_n = std::max((size_t)(b.size()), (size_t)(cstart+ncols));
  b.conservativeResize(new_n);
  // but for now, assume _b is the whole thing
  //b = Eigen::Map<Eigen::Matrix<S, Eigen::Dynamic, 1>>(_b.data(), _b.size());
  // put the pieces where they belong
  for (size_t i=0; i<ncols; ++i) b[cstart+i] = _b[i];

  //std::cout << "RHS size is " << b.size() << std::endl;
  //std::cout << "RHS vector is" << std::endl;
  //std::cout << b << std::endl;
}


//
// Find the change in strength that would occur over one dt
//
template <class S, class I>
void BEM<S,I>::solve() {

  bool verbose = false;

  // ensure that the solution vector is the right size
  strengths.resizeLike(b);

  //std::cout << "A is " << A.size() << std::endl;
  //std::cout << "b is " << b.size() << std::endl;
  //std::cout << "x is " << strengths.size() << std::endl;

  // the Eigen solver object - persistent from call to call
  static Eigen::GMRES<Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> > solver(A);

  if (not solver_initialized) {

    // if A changes, we need to re-run this
    auto istart = std::chrono::system_clock::now();
    solver.compute(A);
    auto iend = std::chrono::system_clock::now();

    std::chrono::duration<double> ielapsed_seconds = iend-istart;
    printf("    solver.init:\t[%.6f] cpu seconds\n", (float)ielapsed_seconds.count());

    solver_initialized = true;
  }

  // note that BiCGSTAB accepts a preconditioner as a template arg
  // note that solveWithGuess() can seed the solution with last step's solution!

  // here is the matrix solution
  auto start = std::chrono::system_clock::now();
  strengths = solver.solve(b);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    solver.solve:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());

  if (false) {
    const size_t nr = 20;
    //const size_t nr = b.size();
    std::cout << "Matrix equation is" << std::endl;
    std::cout << A.block(0,0,nr,8) << std::endl;
    std::cout << "RHS vector is" << std::endl;
    std::cout << b.head(nr) << std::endl;
    std::cout << "Solution vector is" << std::endl;
    std::cout << strengths.head(nr) << std::endl;
  }

  if (verbose) printf("    num iterations:     %d\n", (uint32_t)solver.iterations());
  if (verbose) printf("    estimated error: %g\n", solver.error());

  // find L2 norm of error
  start = std::chrono::system_clock::now();
  double relative_error = (A*strengths - b).norm() / b.norm(); // norm() is L2 norm
  if (verbose) printf("    L2 norm of error is %g\n", relative_error);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  if (verbose) printf("    solver.error:\t[%.6f] cpu seconds\n", (float)elapsed_seconds.count());
}

