// Author: shivanandvp (shivanandvp.oss@gmail.com)

// ===============================
// Includes and global definitions
// ===============================
#include<iostream>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// For the cout statements
#include<cmath>
#include<nlmp_bvp.hpp>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For the function declarations
using namespace std;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** //
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  //
using RowVectorXd = Matrix<double, 1, Dynamic>; // For the convenience of declaring row vectors
// ===============================

// ================
// Functions dxBydt
// ================
// dxBydt = a function that defines the derivative of a state vector x at t -- (nx1)
VectorXd dxBydt(double t, VectorXd x){ 
***REMOVED*** // VectorXd dxdt(6);
***REMOVED*** // double r = sqrt(pow(x(0),2) + pow(x(2),2) + pow(x(4),2));
***REMOVED*** // dxdt(0) = x(1);
***REMOVED*** // dxdt(1) = -x(0) / pow(r,3);
***REMOVED*** // dxdt(2) = x(3);
***REMOVED*** // dxdt(3) = -x(2) / pow(r,3);
***REMOVED*** // dxdt(4) = x(5);
***REMOVED*** // dxdt(5) = -x(4) / pow(r,3);

***REMOVED*** VectorXd dxdt(2);
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = -fabs(x(0));

***REMOVED*** return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at state vectors xBC -- (nx1) 
VectorXd BCResidues(MatrixXd xBC){
***REMOVED*** // VectorXd residues(6);
***REMOVED*** // residues(0) = xBC(0,0) - 1.076;
***REMOVED*** // residues(1) = xBC(1,0) + pow(xBC(1,0),2) + xBC(1,1) + xBC(1,2) + xBC(1,3) + 2.053292953504164;
***REMOVED*** // residues(2) = xBC(2,0) + xBC(3,0) - 0.472283099142472;
***REMOVED*** // residues(3) = pow(xBC(3,1),2) + xBC(4,2) - 1.040204804411078;
***REMOVED*** // residues(4) = xBC(2,3) + xBC(4,3) - 1.57661; 
***REMOVED*** // residues(5) = xBC(5,3) + 0.03407297218269353;

***REMOVED*** VectorXd residues(2);
***REMOVED*** residues(0) = xBC(0,0) - 0;
***REMOVED*** residues(1) = xBC(0,1) - 4;

***REMOVED*** return residues;
}
// ===================

// =================
// The main function
// =================
// This is where the program execution begins
int main(
***REMOVED*** int argc,***REMOVED***// argc = the number of arguments passed to the program
***REMOVED*** char **argv // argv = an array of strings that are passed to  the program 
***REMOVED*** ){

***REMOVED*** cout<<endl;
***REMOVED*** cout<<"Program started..."<<endl;

***REMOVED*** // Variable declarations***REMOVED***
***REMOVED*** // RowVectorXd t_BC(4);***REMOVED******REMOVED******REMOVED***  // t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified -- (1xm)
***REMOVED*** // VectorXd oxt1(6);***REMOVED******REMOVED******REMOVED******REMOVED***  // oxt1***REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED*** VectorXd oxt1(2);
***REMOVED*** RowVectorXd tBC(2);
***REMOVED*** BVPSolution bvpSolution;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters ivamParameters; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** // Variable definitions
***REMOVED*** // t_BC***REMOVED*** << 0.0, 0.8, 1.4, 2;
***REMOVED*** // oxt1 << 1.07600,
***REMOVED*** //***REMOVED******REMOVED******REMOVED******REMOVED***0.53800,
***REMOVED*** //***REMOVED******REMOVED******REMOVED******REMOVED***0.00000,
***REMOVED*** //***REMOVED******REMOVED******REMOVED******REMOVED***0.28800,
***REMOVED*** //***REMOVED******REMOVED******REMOVED******REMOVED***0.00000,
***REMOVED*** //***REMOVED******REMOVED******REMOVED******REMOVED***0.49883;  
***REMOVED*** // oxt1<< -1,
***REMOVED*** //***REMOVED******REMOVED******REMOVED***0;  
***REMOVED*** tBC << 0.0, 4.0;
***REMOVED*** oxt1<< -1,
***REMOVED******REMOVED******REMOVED***  0;
***REMOVED***
***REMOVED*** ivamParameters.EPSILON = 1e-10; // EPSILON = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters.ALPHA***REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters.SIGMA***REMOVED***= 1e-14; // SIGMA***REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters.BETA***REMOVED*** = 1e-3;  // BETA***REMOVED*** = the deflation factor

***REMOVED*** cout<<"tBC = "<<tBC<<endl;
***REMOVED*** cout<<"oxt1 = "<<endl<<oxt1<<endl;

***REMOVED*** cout<<"Initiating the BVP solver..."<<endl;
***REMOVED*** 
***REMOVED*** // Solve the boundary value problem
***REMOVED*** // bvpSolution = nlmp_bvp(6, 4, 501, tBC, oxt1, dxBydt, BCResidues, ivamParameters);
***REMOVED*** bvpSolution = nlmp_bvp(2, 2, 11, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

***REMOVED*** cout<<"Done solving the BVP..."<<endl;

***REMOVED*** // Print the output
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"t = ";
***REMOVED*** cout<<bvpSolution.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"x = "<<endl;
***REMOVED*** cout<<bvpSolution.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"tBC = ";
***REMOVED*** cout<<bvpSolution.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"xBC = "<<endl;
***REMOVED*** cout<<bvpSolution.xBC<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 

***REMOVED*** return 0;
}
// =================