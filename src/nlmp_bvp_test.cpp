// ===============================
// Includes and global definitions
// ===============================
#include<iostream>                              // For the cout statements
#include<nlmp_bvp.hpp>                          // For the function declarations
using namespace std;                            //
using namespace Eigen;                          //
using RowVectorXd = Matrix<double, 1, Dynamic>; // For the convenience of declaring row vectors
// ===============================

// ================
// Functions dxBydt
// ================
// dxBydt = a function that defines the derivative of a state vector x at t -- (nx1)
VectorXd dxBydt(double t, VectorXd x){ 
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
// ================

// ===================
// Functions BCResidue
// ===================
// BCResidue = a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
VectorXd BCResidue(MatrixXd BC){
    VectorXd residues(6);
    residues(0) = BC(0,0) - 1.076;
    residues(1) = BC(1,0) + pow(BC(1,0),2) + BC(1,1) + BC(1,2) + BC(1,3) + 2.053292953504164;
    residues(2) = BC(2,0) + BC(3,0) - 0.472283099142472;
    residues(3) = pow(BC(3,1),2) + BC(4,2) - 1.040204804411078;
    residues(4) = BC(2,3) + BC(4,3) - 1.57661; 
    residues(5) = BC(5,3) + 0.03407297218269353;
    return residues;
}
// ===================

// =================
// The main function
// =================
// This is where the program execution begins
int main(
    int argc,   // argc = the number of arguments passed to the program
    char **argv // argv = an array of strings that are passed to  the program 
    ){

    cout<<endl;
    cout<<"Program started..."<<endl;

    // Variable declarations   
    RowVectorXd t_BC(4);          // t_BC        = row vector of values at which the boundary conditions are specified -- (1xm)
    VectorXd _0x_t1(6);           // _0x_t1      = column vector of the guessed initial state                          -- (nx1)
    BVPSolution bvpSolution;      // bvpSolution = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters ivamParameters // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    // Variable assignments
    t_BC   << 0.0, 0.8, 1.4, 2;
    _0x_t1 << 1.07600,
              0.53800,
              0.00000,
              0.28800,
              0.00000,
              0.49883;    
    ivamParameters.EPSILON = 1e-8;  // EPSILON = the state perturbation parameter to probe the differential equation system with
    ivamParameters.ALPHA   = 1;     // ALPHA   = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters.SIGMA   = 1e-14; // SIGMA   = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters.BETA    = 1e-3;  // BETA    = the deflation factor

    cout<<"Initiating the BVP solver..."<<endl;
    
    // Solve the boundary value problem
    bvpSolution = nlmp_bvp(6, 4, 500, t_BC, _0x_t1, dxBydt, BCResidue, ivamParameters);

    cout<<"Done solving the BVP..."<<endl;

    // Print the output
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"t = ";
    cout<<t<<endl;
    cout<<endl;
    cout<<"x = "<<endl;
    cout<<x<<endl;
    cout<<endl;    
    cout<<"t_BC = ";
    cout<<t_BC<<endl;
    cout<<endl;
    cout<<"x_BC = "<<endl;
    cout<<x_BC<<endl;
    cout<<endl;    

    return 0;
}
// =================