// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================
// Copyright shivanandvp (shivanandvp.oss@gmail.com)

// ==========
// References
// ==========
// [1] Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."
//     Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.
// [2] Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
//     Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

// ===============================
// Includes and global definitions
// ===============================
#include "nlmpbvp_impl.hpp" // Implementation of the solver functions
using namespace Eigen;      // For matrix and vector data types and operations
// ===============================

// ==================
// Function "nlmpBVP"
// ==================
template <typename T>
BVPSolution<T> nlmpBVP(
    int n,                                    // n              = the number of differential equations = the number of boundary conditions
    int m,                                    // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                                // nGrid          = the number of points at which the state is evaluated
    RowVectorXm<T> tBC,                       // tBC            = row vector of values at which the boundary conditions are specified                    -- (1xm)
    VectorXm<T> oxt1,                         // oxt1           = column vector of the guessed initial state                                             -- (nx1)    
    VectorXm<T> dxBydt(T t, VectorXm<T> x),   // dxBydt         = a function that defines the derivative of a state vector x at t                        -- (nx1)
    VectorXm<T> BCResidues(MatrixXm<T> xBC),  // BCResidues     = a function that defines the boundary condition residues at state vectors xBC           -- (nx1) 
    IVAMParameters<T> ivamParameters          // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    );
// ==================

// ===================
// Function "nlmpBVP2"
// ===================
template <typename T> BVPSolution<T> nlmpBVP2(
    int n,                                     // n              = the number of differential equations = the number of boundary conditions
    int m,                                     // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                                 // nGrid          = the number of points at which the state can be evaluated
    RowVectorXm<T> tBC,                        // tBC            = row vector of values at which boundary conditions are specified                       -- (1xm)
    MatrixXm<T> oxt1,                          // oxt1           = matrix of the guessed initial state                                                   -- (nx(m-1))    
    VectorXm<T> dxBydt(T t, VectorXm<T> x),    // dxBydt         = a function that defines the derivative of a state vector x at t                       -- (nx1)
    VectorXm<T> BCResidues(MatrixXm<T> xBCL,   // BCResidues     = a function that defines the boundary condition residues at the nodal state vectors... -- (n(m-1)x1)
                           MatrixXm<T> xBCR),  //                    ... on the left and right side of integration intervals, xBCL and xBCR     
    IVAMParameters<T> ivamParameters           // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    );
// ===================