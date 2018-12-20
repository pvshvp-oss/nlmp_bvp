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
***REMOVED*** VectorXd dxdt(6);
***REMOVED*** double r = sqrt(pow(x(0),2) + pow(x(2),2) + pow(x(4),2));
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = -x(0) / pow(r,3);
***REMOVED*** dxdt(2) = x(3);
***REMOVED*** dxdt(3) = -x(2) / pow(r,3);
***REMOVED*** dxdt(4) = x(5);
***REMOVED*** dxdt(5) = -x(4) / pow(r,3);
***REMOVED*** return dxdt;
}
VectorXd BCFunction(MatrixXd BC){
***REMOVED*** VectorXd residues(6);
***REMOVED*** residues(0) = BC(0,0) - 1.076;
***REMOVED*** residues(1) = BC(1,0) + pow(BC(1,0),2) + BC(1,1) + BC(1,2) + BC(1,3) + 2.053292953504164;
***REMOVED*** residues(2) = BC(2,0) + BC(3,0) - 0.472283099142472;
***REMOVED*** residues(3) = pow(BC(3,1),2) + BC(4,2) - 1.040204804411078;
***REMOVED*** residues(4) = BC(2,3) + BC(4,3) - 1.57661; 
***REMOVED*** residues(5) = BC(5,3) + 0.03407297218269353;
***REMOVED*** return residues;
}
// ============================================================================

// =================
// The main function
// =================
int main(int argc, char **argv)
{
***REMOVED*** VectorXd startingState(6);
***REMOVED*** RowVectorXd tNodes(4);
***REMOVED*** tNodes << 0.0, 0.8, 1.4, 2;
***REMOVED*** startingState << 1.076, 0.538, 0.0, 0.288, 0.0, 0.49883;
***REMOVED*** nlmp_bvp(6, startingState, tNodes, dFunction, BCFunction);
***REMOVED*** return 0;
}
// =================