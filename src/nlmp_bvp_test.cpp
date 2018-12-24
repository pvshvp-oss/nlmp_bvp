// Author: Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)

// ===============================
// Includes and global definitions
// ===============================
#include<iostream>                              // For the cout statements
#include<cmath>
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
    // VectorXd dxdt(6);
    // double r = sqrt(pow(x(0),2) + pow(x(2),2) + pow(x(4),2));
    // dxdt(0) = x(1);
    // dxdt(1) = -x(0) / pow(r,3);
    // dxdt(2) = x(3);
    // dxdt(3) = -x(2) / pow(r,3);
    // dxdt(4) = x(5);
    // dxdt(5) = -x(4) / pow(r,3);

    VectorXd dxdt(2);
    dxdt(0) = x(1);
    dxdt(1) = -fabs(x(0));

    return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at state vectors xBC -- (nx1) 
VectorXd BCResidues(MatrixXd xBC){
    // VectorXd residues(6);
    // residues(0) = xBC(0,0) - 1.076;
    // residues(1) = xBC(1,0) + pow(xBC(1,0),2) + xBC(1,1) + xBC(1,2) + xBC(1,3) + 2.053292953504164;
    // residues(2) = xBC(2,0) + xBC(3,0) - 0.472283099142472;
    // residues(3) = pow(xBC(3,1),2) + xBC(4,2) - 1.040204804411078;
    // residues(4) = xBC(2,3) + xBC(4,3) - 1.57661; 
    // residues(5) = xBC(5,3) + 0.03407297218269353;

    VectorXd residues(2);
    residues(0) = xBC(0,0) - 0;
    residues(1) = xBC(0,1) - 4;

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
    // RowVectorXd t_BC(4);           // t_BC           = row vector of values at which the boundary conditions are specified -- (1xm)
    // VectorXd oxt1(6);              // oxt1        = column vector of the guessed initial state                          -- (nx1)
    VectorXd oxt1(2);
    RowVectorXd tBC(2);
    BVPSolution bvpSolution;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters ivamParameters; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    // Variable definitions
    // t_BC    << 0.0, 0.8, 1.4, 2;
    // oxt1 << 1.07600,
    //            0.53800,
    //            0.00000,
    //            0.28800,
    //            0.00000,
    //            0.49883;  
    // oxt1<< -1,
    //         0;  
    tBC << 0.0, 4.0;
    oxt1<< -1,
           0;
   
    ivamParameters.EPSILON = 1e-10; // EPSILON = the state perturbation parameter to probe the differential equation system with
    ivamParameters.ALPHA   = 1.0;   // ALPHA   = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters.SIGMA   = 1e-14; // SIGMA   = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters.BETA    = 1e-3;  // BETA    = the deflation factor

    cout<<"tBC = "<<tBC<<endl;
    cout<<"oxt1 = "<<endl<<oxt1<<endl;

    cout<<"Initiating the BVP solver..."<<endl;
    
    // Solve the boundary value problem
    // bvpSolution = nlmp_bvp(6, 4, 501, tBC, oxt1, dxBydt, BCResidues, ivamParameters);
    bvpSolution = nlmp_bvp(2, 2, 11, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

    cout<<"Done solving the BVP..."<<endl;

    // Print the output
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"t = ";
    cout<<bvpSolution.t<<endl;
    cout<<endl;
    cout<<"x = "<<endl;
    cout<<bvpSolution.x<<endl;
    cout<<endl;    
    cout<<"tBC = ";
    cout<<bvpSolution.tBC<<endl;
    cout<<endl;
    cout<<"xBC = "<<endl;
    cout<<bvpSolution.xBC<<endl;
    cout<<endl;    

    return 0;
}
// =================