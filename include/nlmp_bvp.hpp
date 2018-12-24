// Author: Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)

// ===============================
// Includes and global definitions
// ===============================
#include <Eigen/Eigen>                          // For the matrix and vector classes
using namespace Eigen;                          //
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
    mpreal EPSILON; // EPSILON = the state perturbation parameter to probe the differential equation system with
    mpreal ALPHA;   // ALPHA   = the relaxation factor to scale the adjustment to the initial condition
    mpreal SIGMA;   // SIGMA   = the tolerance for error outside which the solver needs to  iterate further. 
    mpreal BETA;    // BETA    = the deflation factor
};
// ========================

// =====================
// Structure BVPSolution
// =====================
struct BVPSolution{
    RowVectorXm t;    // t    = row vector of values at which state vectors x are evaluated                       -- (1xnGrid)
    MatrixXm x;       // x    = state vectors at all grid points t                                                -- (nxnGrid)  
    RowVectorXm tBC; // tBC = row vector of boundary condition values at which state vectors xBC are evaluated -- (1xm)   
    MatrixXm xBC;    // xBC = state vectors at boundary conditions tBC                                         -- (nxm)
};
// =====================

// ===================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
    int n,                                 // n              = the number of differential equations = the number of boundary conditions
    int m,                                 // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                             // nGrid          = the number of points at which the state is evaluated
    RowVectorXm tBC,                      // tBC           = row vector of values at which the boundary conditions are specified           -- (1xm)
    VectorXm oxt1,                      // oxt1        = column vector of the guessed initial state                                    -- (nx1)    
    VectorXm dxBydt(mpreal t, VectorXm x), // dxBydt         = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXm BCResidues(MatrixXm xBC),    // BCResidues     = a function that defines the boundary condition residues at state vectors xBC -- (nx1) 
    IVAMParameters ivamParameters          // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    );
// ===================
