# nlmp_bvp

**Nonlinear Multipoint Boundary Value Problem Solver** in `C++` based on:  
>>>
[1] Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."  
Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.  
[2] Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."  
Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.  
This was in response to a need for a suitable C++ alternative to bvp4c in Mathworks MATLAB.  
>>>

```cpp
// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================
// Copyright shivanandvp (shivanandvp.oss@gmail.com)

// ==========
// References
// ==========
// [1] Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."
//***REMOVED***  Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.
// [2] Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
//***REMOVED***  Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

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

/* Boundary value problem 3 */
// Example 1
// Welsh, Wayne, and Takeo Ojika.
// "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
// Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

/* Boundary value problem 4 */
// Example 3
// Welsh, Wayne, and Takeo Ojika.
// "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
// Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.
// ============================


// ===============================
// Includes and global definitions
// ===============================
#include <iostream>***REMOVED******REMOVED******REMOVED******REMOVED*** // For the cout statements
#include "nlmpbvp.hpp"***REMOVED******REMOVED******REMOVED*** // For the boundary value problem solver function declarations
#include <Eigen/MPRealSupport>  // For arbitrary precision computation
using namespace std;***REMOVED******REMOVED******REMOVED******REMOVED***// For cout
using namespace Eigen;***REMOVED******REMOVED******REMOVED*** // For matrix and vector data types and operations
using namespace mpfr;***REMOVED******REMOVED******REMOVED***  // For arbitrary precision computation
// ===============================

// ================
// Functions dxBydt
// ================
// dxBydt = a function that defines the derivative of a state vector x at t -- (nx1)

/* Boundary Value Problem 1 */
VectorXm<mpreal> dxBydt_BVP1(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(2);
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = -fabs(x(0));
***REMOVED*** return dxdt;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> dxBydt_BVP2(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(3);
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = x(2);
***REMOVED*** dxdt(2) = 25*x(1) - 1;
***REMOVED*** return dxdt;
}

/* Boundary Value Problem 3 */
VectorXm<mpreal> dxBydt_BVP3(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(2);
***REMOVED*** dxdt(0) =  x(1);
***REMOVED*** dxdt(1) = -x(0);
***REMOVED*** return dxdt;
}

/* Boundary Value Problem 4 */
VectorXm<mpreal> dxBydt_BVP4(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(2);
***REMOVED*** dxdt(0) =  x(1);
***REMOVED*** dxdt(1) = -x(0);
***REMOVED*** return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at nodal state vectors xBC -- (nx1) 

/* Boundary Value Problem 1 */
VectorXm<mpreal> BCResidues_BVP1(MatrixXm<mpreal> xBC){
***REMOVED*** VectorXm<mpreal> residues(2);
***REMOVED*** residues(0) = xBC(0,0) - 0;
***REMOVED*** residues(1) = xBC(0,1) + 2;
***REMOVED*** return residues;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> BCResidues_BVP2(MatrixXm<mpreal> xBC){
***REMOVED*** VectorXm<mpreal> residues(3);
***REMOVED*** residues(0) = xBC(1,0) - 0;
***REMOVED*** residues(1) = xBC(1,2) - 0;
***REMOVED*** residues(2) = xBC(0,1) - 0;
***REMOVED*** return residues;
}


// BCResidues***REMOVED***  = a function that defines the boundary condition residues at the nodal state vectors... -- (n(m-1)x1)
//***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***... on the left and right side of integration intervals, xBCL and xBCR

/* Boundary Value Problem 3 */
VectorXm<mpreal> BCResidues_BVP3(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
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

/* Boundary Value Problem 4 */
VectorXm<mpreal> BCResidues_BVP4(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
***REMOVED*** VectorXm<mpreal> residues(8);
***REMOVED*** residues(0) = xBCL(0,0)-0;
***REMOVED*** residues(1) = xBCL(1,0)-1;
***REMOVED*** residues(2) = pow(xBCL(0,1),2)*pow(xBCL(1,1),3) - exp(xBCR(1,1))*pow(xBCR(1,0),2)*pow(xBCL(1,3),2) + pow(xBCR(1,3),2)*pow(xBCR(0,2),2) - 0.1859767072;
***REMOVED*** residues(3) = pow(xBCR(0,1),2)*pow(xBCR(0,0),3)*pow(xBCL(0,3),2) + pow(xBCL(1,2),2)*exp(xBCR(1,2)) + pow(xBCL(0,2),2)*pow(xBCR(0,3),2) - 0.1261677772;
***REMOVED*** residues(4) = xBCR(0,0) - xBCL(0,1) + 0;
***REMOVED*** residues(5) = xBCR(1,0) - xBCL(1,1) + 0;
***REMOVED*** residues(6) = xBCR(0,2) - xBCL(0,3) + 0;
***REMOVED*** residues(7) = xBCR(1,2) - xBCL(1,3) + 0;
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
***REMOVED*** cout<<"============================================================================="<<endl;
***REMOVED*** cout<<"Test: Non-linear multipoint boundary value problem solver (nlmpBVP, nlmpBVP2)"<<endl;
***REMOVED*** cout<<"============================================================================="<<endl;
***REMOVED*** cout<<"Copyright shivanandvp (shivanandvp.oss@gmail.com)"<<endl;

***REMOVED*** // Variable declarations***REMOVED***

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP1(2);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm<mpreal>***REMOVED***oxt1_BVP1(2);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP2(3);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm<mpreal>***REMOVED***oxt1_BVP2(3);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP3(6);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)  
***REMOVED*** MatrixXm<mpreal> oxt1_BVP3(2,5);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = matrix of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))

***REMOVED*** /* Boundary Value Problem 4 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP4(5);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** MatrixXm<mpreal> oxt1_BVP4(2,4);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = matrix of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP1;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP1; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP2;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP2; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP3;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP3; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP4;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP4; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** // Variable definitions

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** tBC_BVP1  << 0.0, 4.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** oxt1_BVP1 <<  1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1) 
***REMOVED******REMOVED******REMOVED******REMOVED*** 0;  
***REMOVED*** // tBC_BVP1  << 0.0, 4.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** // oxt1_BVP1 << -1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1) 
***REMOVED*** //***REMOVED******REMOVED******REMOVED*** 0;

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** tBC_BVP2  << 0.0, 0.5, 1.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** oxt1_BVP2 << 1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED******REMOVED******REMOVED***1,
***REMOVED******REMOVED******REMOVED******REMOVED***1;

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** tBC_BVP3  << 0.0, mpfr::const_pi()/6, mpfr::const_pi()/3, mpfr::const_pi()/2, 2*mpfr::const_pi()/3, mpfr::const_pi();
***REMOVED*** // oxt1 = a matrix of the guessed initial state on the left side of each integration interval***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))
***REMOVED*** oxt1_BVP3 <<  0.1, 0.1, 0.4, 0.8, 0.9,
***REMOVED******REMOVED******REMOVED******REMOVED***-0.6, 0.1, 0.9, 2.1, 0.8;

***REMOVED*** /* Boundary Value Problem 4 */
***REMOVED*** // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** tBC_BVP4  << 0.0, mpfr::const_pi()/4, mpfr::const_pi()/2, 3*mpfr::const_pi()/4, mpfr::const_pi();
***REMOVED*** // oxt1 = a matrix of the guessed initial state on the left side of each integration interval***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))
***REMOVED*** oxt1_BVP4 << -0.1, 0.7, 0.1,-0.7,
***REMOVED******REMOVED******REMOVED******REMOVED*** 1.1, 0.7,-1.1,-0.7;

***REMOVED*** // Assign the parameters for IVAM

***REMOVED*** ivamParameters_BVP1.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP1.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP1.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP1.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP1.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** ivamParameters_BVP2.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP2.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP2.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP2.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP2.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** ivamParameters_BVP3.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP3.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP3.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP3.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP3.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** ivamParameters_BVP4.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP4.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP4.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP4.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP4.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 1"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP1<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP1<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP1 = nlmpBVP<mpreal>(2, 2, 101, tBC_BVP1, oxt1_BVP1, dxBydt_BVP1, BCResidues_BVP1, ivamParameters_BVP1);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP1.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP1.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP1.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP1.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED***  /* Boundary Value Problem 2 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 2"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP2<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP2<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP2 = nlmpBVP<mpreal>(3, 3, 101, tBC_BVP2, oxt1_BVP2, dxBydt_BVP2, BCResidues_BVP2, ivamParameters_BVP2);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP2.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP2.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP2.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP2.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 3"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP3<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP3<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP3 = nlmpBVP2<mpreal>(2, 6, 12*10+1, tBC_BVP3, oxt1_BVP3, dxBydt_BVP3, BCResidues_BVP3, ivamParameters_BVP3);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP3.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP3.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP3.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP3.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED*** /* Boundary Value Problem 4 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 4"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP4<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP4<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP4 = nlmpBVP2<mpreal>(2, 5, 10*10+1, tBC_BVP4, oxt1_BVP4, dxBydt_BVP4, BCResidues_BVP4, ivamParameters_BVP4);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP4.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP4.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP4.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP4.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED*** cout<<endl<<"Program ended..."<<endl<<endl;

***REMOVED*** return 0;
}
// =================
```

