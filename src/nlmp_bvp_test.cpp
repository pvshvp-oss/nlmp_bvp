// ========================
// Includes and definitions
// ========================
#include<iostream>
#include<nlmp_bvp.hpp>
using namespace std;
// ========================

// ===========================================================================
// Functions that define the differential equation and its boundahy conditions
// ===========================================================================
VectorXd dFunction(double t, VectorXd x){
    VectorXd dxdt(2);
    dxdt(0) = x(1);
    dxdt(1) = 0.01*(1.0 - pow(x(0),2.0))*x(1) - x(0);
    return dxdt;
}
VectorXd BCFunction(MatrixXd BC){
    VectorXd residues(2);
    residues(0) = 4.0*BC(0,0) + pow(BC(1,1),2.0) - 8.0;
    residues(1) = pow(BC(0,0),2.0) + 2.0*BC(1,1) - 5.0;
    return residues;
}
// ============================================================================

// =================
// The main function
// =================
int main(int argc, char **argv)
{
    VectorXd startingState;
    VectorXd tNodes;
    nlmp_bvp(2, startingState, tNodes, dFunction, BCFunction);
    cout << "Hello" << endl;
    return 0;
}
// =================