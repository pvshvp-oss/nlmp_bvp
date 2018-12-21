// ===============================
// Includes and global definitions
// ===============================
#include <Eigen/Eigen>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For the matrix and vector classes
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  //
using RowVectorXd = Matrix<double, 1, Dynamic>; // For the convenience of declaring row vectors
===============================

// ========================
// Structure IVAMParameters
// ========================
// Parameters for the Initial Value Adjusting Method (IVAM)
struct IVAMParameters{
***REMOVED*** double EPSILON; // EPSILON = the state perturbation parameter to probe the differential equation system with
***REMOVED*** double ALPHA;***REMOVED***// ALPHA***REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** double SIGMA;***REMOVED***// SIGMA***REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** double BETA;***REMOVED*** // BETA***REMOVED*** = the deflation factor
}
// ========================

// =====================
// Structure BVPSolution
// =====================
struct BVPSolution{
***REMOVED*** RowVectorXd t;***REMOVED*** // t***REMOVED*** = row vector of values at which state vectors x are evaluated***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xnGrid)
***REMOVED*** MatrixXd x;***REMOVED******REMOVED*** // x***REMOVED*** = state vectors at all grid points t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nxnGrid)  
***REMOVED*** RowVectorXd t_BC; // t_BC = row vector of boundary condition values at which state vectors x_BC are evaluated -- (1xm)***REMOVED***
***REMOVED*** MatrixXd x_BC;***REMOVED*** // x_BC = state vectors at boundary conditions t_BC***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nxm)
}
// =====================

// ===================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state is evaluated
***REMOVED*** RowVectorXd t_BC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXd _0x_t1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // _0x_t1***REMOVED******REMOVED******REMOVED***= column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXd dxBydt(double t, VectorXd x), // dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXd BCResidue(MatrixXd x_BC)***REMOVED******REMOVED***// BCResidue***REMOVED******REMOVED***= a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
***REMOVED*** const IVAMParameters ivamParameters***REMOVED*** // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** );
// ===================