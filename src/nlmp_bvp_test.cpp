// ===============================
// Includes and global definitions
// ===============================
#include<iostream>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// For the cout statements
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
// ================

// ===================
// Functions BCResidue
// ===================
// BCResidue = a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
VectorXd BCResidue(MatrixXd BC){
***REMOVED*** VectorXd residues(6);
***REMOVED*** residues(0) = BC(0,0) - 1.076;
***REMOVED*** residues(1) = BC(1,0) + pow(BC(1,0),2) + BC(1,1) + BC(1,2) + BC(1,3) + 2.053292953504164;
***REMOVED*** residues(2) = BC(2,0) + BC(3,0) - 0.472283099142472;
***REMOVED*** residues(3) = pow(BC(3,1),2) + BC(4,2) - 1.040204804411078;
***REMOVED*** residues(4) = BC(2,3) + BC(4,3) - 1.57661; 
***REMOVED*** residues(5) = BC(5,3) + 0.03407297218269353;
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
***REMOVED*** RowVectorXd t_BC(4);***REMOVED***  // t_BC***REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified -- (1xm)
***REMOVED*** VectorXd _0x_t1(6);***REMOVED******REMOVED***// _0x_t1***REMOVED******REMOVED***= column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED*** BVPSolution bvpSolution; // bvpSolution = the structure in which the solutions of the boundary value problem will be saved

***REMOVED*** // Variable assignments
***REMOVED*** t_BC***REMOVED***<< 0.0, 0.8, 1.4, 2;
***REMOVED*** _0x_t1 << 1.07600,
***REMOVED******REMOVED******REMOVED******REMOVED***  0.53800,
***REMOVED******REMOVED******REMOVED******REMOVED***  0.00000,
***REMOVED******REMOVED******REMOVED******REMOVED***  0.28800,
***REMOVED******REMOVED******REMOVED******REMOVED***  0.00000,
***REMOVED******REMOVED******REMOVED******REMOVED***  0.49883;

***REMOVED*** cout<<"Initiating the BVP solver..."<<endl;
***REMOVED*** 
***REMOVED*** // Solve the boundary value problem
***REMOVED*** bvpSolution = nlmp_bvp(6, 4, 500, t_BC, _0x_t1, dxBydt, BCResidue);

***REMOVED*** cout<<"Done solving the BVP..."<<endl;

***REMOVED*** // Print the output
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"t = ";
***REMOVED*** cout<<t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"x = "<<endl;
***REMOVED*** cout<<x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"t_BC = ";
***REMOVED*** cout<<t_BC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"x_BC = "<<endl;
***REMOVED*** cout<<x_BC<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 

***REMOVED*** return 0;
}
// =================