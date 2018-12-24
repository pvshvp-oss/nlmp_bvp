// Author: Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)

// ===============================
// Includes and global definitions
// ===============================
#include <Eigen/Eigen>                          // For the matrix and vector classes
using namespace Eigen;                          //
using RowVectorXd = Matrix<double, 1, Dynamic>; // For the convenience of declaring row vectors
// ===============================

// ========================
// Structure IVAMParameters
// ========================
// Parameters for the Initial Value Adjusting Method (IVAM)
struct IVAMParameters{
    double EPSILON; // EPSILON = the state perturbation parameter to probe the differential equation system with
    double ALPHA;   // ALPHA   = the relaxation factor to scale the adjustment to the initial condition
    double SIGMA;   // SIGMA   = the tolerance for error outside which the solver needs to  iterate further. 
    double BETA;    // BETA    = the deflation factor
};
// ========================

// =====================
// Structure BVPSolution
// =====================
struct BVPSolution{
    RowVectorXd t;    // t    = row vector of values at which state vectors x are evaluated                       -- (1xnGrid)
    MatrixXd x;       // x    = state vectors at all grid points t                                                -- (nxnGrid)  
    RowVectorXd tBC; // tBC = row vector of boundary condition values at which state vectors xBC are evaluated -- (1xm)   
    MatrixXd xBC;    // xBC = state vectors at boundary conditions tBC                                         -- (nxm)
};
// =====================

// ===================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
    int n,                                 // n              = the number of differential equations = the number of boundary conditions
    int m,                                 // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                             // nGrid          = the number of points at which the state is evaluated
    RowVectorXd tBC,                      // tBC           = row vector of values at which the boundary conditions are specified           -- (1xm)
    VectorXd oxt1,                      // oxt1        = column vector of the guessed initial state                                    -- (nx1)    
    VectorXd dxBydt(double t, VectorXd x), // dxBydt         = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXd BCResidues(MatrixXd xBC),    // BCResidues     = a function that defines the boundary condition residues at state vectors xBC -- (nx1) 
    IVAMParameters ivamParameters          // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    );
// ===================