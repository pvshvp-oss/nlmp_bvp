// ========================
// Includes and definitions
// ========================
#include<iostream>
#include<nlmp_bvp.hpp>
using namespace std;
using namespace Eigen;
using RowVectorXd = Matrix<double, 1, Dynamic>;
// ========================

// ===========================================================================
// Functions that define the differential equation and its boundary conditions
// ===========================================================================
VectorXd dFunction(double t, VectorXd x){
***REMOVED*** VectorXd dxdt(2);
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = 0.01*(1.0 - pow(x(0),2.0))*x(1) - x(0);
***REMOVED*** return dxdt;
}
VectorXd BCFunction(MatrixXd BC){
***REMOVED*** VectorXd residues(2);
***REMOVED*** residues(0) = 4.0*BC(0,0) + pow(BC(1,1),2.0) - 8.0;
***REMOVED*** residues(1) = pow(BC(0,0),2.0) + 2.0*BC(1,1) - 5.0;
***REMOVED*** return residues;
}
// ============================================================================

// =================
// The main function
// =================
int main(int argc, char **argv)
{
***REMOVED*** VectorXd startingState(2);
***REMOVED*** RowVectorXd tNodes(2);
***REMOVED*** tNodes << 0.0, 1.0;
***REMOVED*** startingState << 0.5, 0.5;
***REMOVED*** nlmp_bvp(2, startingState, tNodes, dFunction, BCFunction);
***REMOVED*** return 0;
}
// =================