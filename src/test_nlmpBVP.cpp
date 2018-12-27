// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================

// ============================
// Credits for example problems
// ============================
/* Boundary Value Problem 1 */
// "Using Initial Guess to Indicate Desired Solution"
// “bvp4c.” Mathworks Documentation, Mathworks, www.mathworks.com/help/matlab/ref/bvp4c.html

/* Boundary value problem 2 */
// Example 2
// Tr, Ramesh. (2017). A novel method for solving multipoint boundary value problems.
// Global Journal of Pure and Applied Mathematics. 13. 850-857. 
// ============================

// ===============================
// Includes and global definitions
// ===============================
#include <iostream>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// For the cout statements
#include "nlmpbvp.hpp"***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// For the boundary value problem solver function declarations
#include <Eigen/MPRealSupport>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For arbitrary precision computation
using namespace std;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For cout
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// For matrix and vector data types and operations
using namespace mpfr;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For arbitrary precision computation
// ===============================

// ================
// Functions dxBydt
// ================
// dxBydt = a function that defines the derivative of a state vector x at t -- (nx1)
VectorXm<mpreal> dxBydt(mpreal t, VectorXm<mpreal> x){ 

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** // VectorXm<mpreal> dxdt(2);
***REMOVED*** // dxdt(0) = x(1);
***REMOVED*** // dxdt(1) = -fabs(x(0));

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** // VectorXm<mpreal> dxdt(3);
***REMOVED*** // dxdt(0) = x(1);
***REMOVED*** // dxdt(1) = x(2);
***REMOVED*** // dxdt(2) = 25*x(1) - 1;

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** VectorXm<mpreal> dxdt(2);
***REMOVED*** dxdt(0) =  x(1);
***REMOVED*** dxdt(1) = -x(0);

***REMOVED*** return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at state vectors xBC -- (nx1) 
VectorXm<mpreal> BCResidues(MatrixXm<mpreal> xBC){
***REMOVED*** 
***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** // VectorXm<mpreal> residues(2);
***REMOVED*** // residues(0) = xBC(0,0) - 0;
***REMOVED*** // residues(1) = xBC(0,1) + 2;

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** VectorXm<mpreal> residues(3);
***REMOVED*** residues(0) = xBC(1,0) - 0;
***REMOVED*** residues(1) = xBC(1,2) - 0;
***REMOVED*** residues(2) = xBC(0,1) - 0;

***REMOVED*** return residues;
}
// ===================

// ====================
// Functions BCResidues
// ====================
// BCResidues***REMOVED***  = a function that defines the boundary condition residues...  -- (n(m-1)x1)
//***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***...at the left and right state vectors xBCL and xBCR
VectorXm<mpreal> BCResidues(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** VectorXm<mpreal> residues(10);
***REMOVED*** residues(0) = xBCL(0,0)-0;
***REMOVED*** residues(1) = xBCL(1,0)-1;
***REMOVED*** residues(2) = xBCR(0,0) - xBCL(0,1) - 1;
***REMOVED*** residues(3) = xBCR(1,0) - xBCL(1,1) - 0;
***REMOVED*** residues(4) = xBCR(1,1) - xBCL(1,2) + 1;
***REMOVED*** residues(5) = xBCR(0,1) - xBCL(0,2) + 0;
***REMOVED*** residues(6) = xBCR(0,2) - xBCL(0,3) - 1;
***REMOVED*** residues(7) = xBCR(1,2) - xBCL(1,3) - (sqrt(3)-1);
***REMOVED*** residues(8) = xBCR(0,3) - xBCL(0,4) - 0;
***REMOVED*** residues(9) = xBCR(1,3) - xBCL(1,4) - 0;

***REMOVED*** return residues;
}


// =================
// The main function
// =================
// This is where the program execution begins
int main(
***REMOVED*** int argc,***REMOVED***// argc = the number of arguments passed to the program
***REMOVED*** char **argv // argv = an array of strings that are passed to  the program 
***REMOVED*** ){

***REMOVED*** mpreal::set_default_prec(64); // Set the number of digits of precision you want for computations

***REMOVED*** cout<<endl;
***REMOVED*** cout<<"===================================================================="<<endl;
***REMOVED*** cout<<"Test: Non-linear multipoint boundary value problem solver (nlmpBVP)"<<endl;
***REMOVED*** cout<<"===================================================================="<<endl;
***REMOVED*** cout<<"Copyright shivanandvp (shivanandvp.oss@gmail.com)"<<endl;

***REMOVED*** // Variable declarations***REMOVED***

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** // RowVectorXm<mpreal> tBC(2);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** // VectorXm<mpreal>***REMOVED***oxt1(2);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** // RowVectorXm<mpreal> tBC(3);***REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** // VectorXm<mpreal>***REMOVED***oxt1(3);***REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** RowVectorXm<mpreal> tBC(6);
***REMOVED*** MatrixXm<mpreal>***REMOVED***oxt1(2,5); 

***REMOVED*** BVPSolution<mpreal> bvpSolution;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** // Variable definitions

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** // tBC  << 0.0, 4.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** // oxt1 <<  1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1) 
***REMOVED******REMOVED******REMOVED******REMOVED***//  0;  
***REMOVED*** // tBC  << 0.0, 4.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** // oxt1 << -1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1) 
***REMOVED*** //***REMOVED******REMOVED******REMOVED*** 0;

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** // tBC  << 0.0, 0.5, 1.0;***REMOVED******REMOVED******REMOVED*** // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** // oxt1 << 1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED*** //***REMOVED******REMOVED******REMOVED***1,
***REMOVED*** //***REMOVED******REMOVED******REMOVED***1;

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** tBC  << 0.0, mpfr::const_pi()/6, mpfr::const_pi()/3, mpfr::const_pi()/2, 2*mpfr::const_pi()/3, mpfr::const_pi();
***REMOVED*** oxt1 <<  0.1, 0.1, 0.4, 0.8, 0.9,
***REMOVED******REMOVED******REMOVED******REMOVED***-0.6, 0.1, 0.9, 2.1, 0.8;

***REMOVED*** // Assign the parameters for IVAM
***REMOVED*** ivamParameters.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters.printDebug = true; // printDebug = specify whether debug messages should be output to the console

***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1<<endl;

***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** 
***REMOVED*** // Solve the boundary value problem

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** // bvpSolution = nlmpBVP<mpreal>(2, 2, 101, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** // bvpSolution = nlmpBVP<mpreal>(3, 3, 101, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** bvpSolution = nlmpBVP2<mpreal>(2, 6, 12*10+1, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

***REMOVED*** cout<<endl<<"Done solving the BVP..."<<endl;

***REMOVED*** // Print the output
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution.xBC<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 

***REMOVED*** return 0;
}
// =================