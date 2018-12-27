// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================
// Copyright shivanandvp (shivanandvp.oss@gmail.com)

// ==========
// References
// ==========
// [1] Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."
//***REMOVED***  Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.
// [2] Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
//***REMOVED***  Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

// ===============================
// Includes and global definitions
// ===============================
#include "nlmpbvp_impl.hpp" // Implementation of the solver function
using namespace Eigen;***REMOVED******REMOVED***// For matrix and vector data types and operations
// ===============================

// ==================
// Function "nlmpBVP"
// ==================
template <typename T>
BVPSolution<T> nlmpBVP(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state is evaluated
***REMOVED*** RowVectorXm<T> tBC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // tBC***REMOVED******REMOVED******REMOVED******REMOVED***= row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm<T> oxt1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXm<T> dxBydt(T t, VectorXm<T> x),***REMOVED***// dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXm<T> BCResidues(MatrixXm<T> xBC),  // BCResidues***REMOVED***  = a function that defines the boundary condition residues at state vectors xBC***REMOVED******REMOVED******REMOVED***  -- (nx1) 
***REMOVED*** IVAMParameters<T> ivamParameters***REMOVED******REMOVED******REMOVED*** // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** );
// ==================

// ===================
// Function "nlmpBVP2"
// ===================
template <typename T> BVPSolution<T> nlmpBVP2(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state can be evaluated
***REMOVED*** RowVectorXm<T> tBC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tBC***REMOVED******REMOVED******REMOVED******REMOVED***= row vector of values at which boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** MatrixXm<T> oxt1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // oxt1***REMOVED******REMOVED******REMOVED***  = matrix of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx(m-1))***REMOVED*** 
***REMOVED*** VectorXm<T> dxBydt(T t, VectorXm<T> x),***REMOVED*** // dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED*** VectorXm<T> BCResidues(MatrixXm<T> xBCL,***REMOVED***// BCResidues***REMOVED***  = a function that defines the boundary condition residues at the nodal state vectors... -- (n(m-1)x1)
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***MatrixXm<T> xBCR),  //***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  ... on the left and right side of integration intervals, xBCL and xBCR***REMOVED***  
***REMOVED*** IVAMParameters<T> ivamParameters***REMOVED******REMOVED******REMOVED***  // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** );
// ===================