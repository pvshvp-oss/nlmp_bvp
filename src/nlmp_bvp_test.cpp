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
    VectorXd dxdt(6);
    double r = sqrt(pow(x(0),2) + pow(x(2),2) + pow(x(4),2));
    dxdt(0) = x(1);
    dxdt(1) = -x(0) / pow(r,3);
    dxdt(2) = x(3);
    dxdt(3) = -x(2) / pow(r,3);
    dxdt(4) = x(5);
    dxdt(5) = -x(4) / pow(r,3);
    return dxdt;
}
VectorXd BCFunction(MatrixXd BC){
    VectorXd residues(6);
    residues(0) = BC(0,0) - 1.076;
    residues(1) = BC(1,0) + pow(BC(1,0),2) + BC(1,1) + BC(1,2) + BC(1,3) + 2.053292953504164;
    residues(2) = BC(2,0) + BC(3,0) - 0.472283099142472;
    residues(3) = pow(BC(3,1),2) + BC(4,2) - 1.040204804411078;
    residues(4) = BC(2,3) + BC(4,3) - 1.57661; 
    residues(5) = BC(5,3) + 0.03407297218269353;
    return residues;
}
// ============================================================================

// =================
// The main function
// =================
int main(int argc, char **argv)
{
    VectorXd startingState(6);
    RowVectorXd tNodes(4);
    tNodes << 0.0, 0.8, 1.4, 2;
    startingState << 1.076, 0.538, 0.0, 0.288, 0.0, 0.49883;
    nlmp_bvp(6, startingState, tNodes, dFunction, BCFunction);
    return 0;
}
// =================