// Author: shivanandvp (shivanandvp.oss@gmail.com)

// ===============================
// Includes and global definitions
// ===============================
#include <Eigen/Eigen>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For the matrix and vector classes
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  //
#include <Eigen/MPRealSupport>
using namespace mpfr;
using MatrixXm = Matrix<mpreal, Dynamic, Dynamic>;
using VectorXm = Matrix<mpreal, Dynamic, 1>;
using RowVectorXm = Matrix<mpreal, 1, Dynamic>; // For the convenience of declaring row vectors
// ===============================

// ========================
// Structure IVAMParameters
// ========================
// Parameters for the Initial Value Adjusting Method (IVAM)
struct IVAMParameters{
***REMOVED*** mpreal EPSILON; // EPSILON = the state perturbation parameter to probe the differential equation system with
***REMOVED*** mpreal ALPHA;***REMOVED***// ALPHA***REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** mpreal SIGMA;***REMOVED***// SIGMA***REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** mpreal BETA;***REMOVED*** // BETA***REMOVED*** = the deflation factor
};
// ========================

// =====================
// Structure BVPSolution
// =====================
struct BVPSolution{
***REMOVED*** RowVectorXm t;***REMOVED*** // t***REMOVED*** = row vector of values at which state vectors x are evaluated***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xnGrid)
***REMOVED*** MatrixXm x;***REMOVED******REMOVED*** // x***REMOVED*** = state vectors at all grid points t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nxnGrid)  
***REMOVED*** RowVectorXm tBC; // tBC = row vector of boundary condition values at which state vectors xBC are evaluated -- (1xm)***REMOVED***
***REMOVED*** MatrixXm xBC;***REMOVED*** // xBC = state vectors at boundary conditions tBC***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nxm)
};
// =====================

// ===================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state is evaluated
***REMOVED*** RowVectorXm tBC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // tBC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm oxt1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // oxt1***REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXm dxBydt(mpreal t, VectorXm x), // dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXm BCResidues(MatrixXm xBC),***REMOVED*** // BCResidues***REMOVED***  = a function that defines the boundary condition residues at state vectors xBC -- (nx1) 
***REMOVED*** IVAMParameters ivamParameters***REMOVED******REMOVED******REMOVED*** // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** );
// ===================
