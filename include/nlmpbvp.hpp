// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================

// ===============================
// Includes and global definitions
// ===============================
#include "nlmpbvp_impl.hpp"***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Implementation of the solver function
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For matrix and vector data types and operations
// ===============================

// ==================
// Function "nlmpBVP"
// ==================
template <typename T>
BVPSolution<T> nlmpBVP(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state is evaluated
***REMOVED*** RowVectorXm<T> tBC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // tBC***REMOVED******REMOVED******REMOVED******REMOVED***= row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm<T> oxt1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXm<T> dxBydt(T t, VectorXm<T> x),  // dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXm<T> BCResidues(MatrixXm<T> xBC), // BCResidues***REMOVED***  = a function that defines the boundary condition residues at state vectors xBC  -- (nx1) 
***REMOVED*** IVAMParameters<T> ivamParameters***REMOVED******REMOVED******REMOVED***// ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** );
// ==================
