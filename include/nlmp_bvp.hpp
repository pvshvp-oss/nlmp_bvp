// ===============================
// Includes and global definitions
// ===============================
#include <Eigen/Eigen>                          // For the matrix and vector classes
using namespace Eigen;                          //
using RowVectorXd = Matrix<double, 1, Dynamic>; // For the convenience of declaring row vectors
===============================

// ===================
// Function "nlmp_bvp"
// ===================
int nlmp_bvp(
    int n,                                 // n         = the number of differential equations = the number of boundary conditions
    int m,                                 // m         = the number of nodes at which boundary conditions are specified
    int nGrid,                             // nGrid     = the number of points at which the state is evaluated
    RowVectorXd t_BC,                      // t_BC      = row vector of values at which the boundary conditions are specified           -- (1xm)
    VectorXd _0x_t1,                       // _0x_t1    = column vector of the guessed initial state                                    -- (nx1)    
    VectorXd dxBydt(double t, VectorXd x), // dxBydt    = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXd BCResidue(MatrixXd x_BC)      // BCResidue = a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
    );
// ===================

// =====================
// Structure BVPSolution
// =====================
struct BVPSolution{
    RowVectorXd t;    // t    = row vector of values at which state vectors x are evaluated                       -- (1xnGrid)
    MatrixXd x;       // x    = state vectors at all grid points t                                                -- (nxnGrid)  
    RowVectorXd t_BC; // t_BC = row vector of boundary condition values at which state vectors x_BC are evaluated -- (1xm)   
    MatrixXd x_BC;    // x_BC = state vectors at boundary conditions t_BC                                         -- (nxm)
}
// =====================