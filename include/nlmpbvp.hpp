// ========================================
// Author: Shivanand Pattanshetti 
// Email : shivanand.pattanshetti@gmail.com
// ========================================

// ===============================
// Includes and global definitions
// ===============================
#include "nlmpbvp_impl.hpp"                                           // Implementation of the solver function
using namespace Eigen;                                                 // For matrix and vector data types and operations
// ===============================

// ==================
// Function "nlmpBVP"
// ==================
template <typename T>
BVPSolution<T> nlmpBVP(
    int n,                                   // n              = the number of differential equations = the number of boundary conditions
    int m,                                   // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                               // nGrid          = the number of points at which the state is evaluated
    RowVectorXm<T> tBC,                      // tBC            = row vector of values at which the boundary conditions are specified           -- (1xm)
    VectorXm<T> oxt1,                        // oxt1           = column vector of the guessed initial state                                    -- (nx1)    
    VectorXm<T> dxBydt(T t, VectorXm<T> x),  // dxBydt         = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXm<T> BCResidues(MatrixXm<T> xBC), // BCResidues     = a function that defines the boundary condition residues at state vectors xBC  -- (nx1) 
    IVAMParameters<T> ivamParameters         // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    );
// ==================
