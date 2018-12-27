// ========================================
// Author: Shivanand Pattanshetti 
// Email : shivanand.pattanshetti@gmail.com
// ========================================
// Copyright Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)

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
// ============================

// ===============================
// Includes and global definitions
// ===============================
#include <iostream>             // For the cout statements
#include "nlmpbvp.hpp"          // For the boundary value problem solver function declarations
#include <Eigen/MPRealSupport>  // For arbitrary precision computation
using namespace std;            // For cout
using namespace Eigen;          // For matrix and vector data types and operations
using namespace mpfr;           // For arbitrary precision computation
// ===============================

// ================
// Functions dxBydt
// ================
// dxBydt = a function that defines the derivative of a state vector x at t                                -- (nx1)
VectorXm<mpreal> dxBydt(mpreal t, VectorXm<mpreal> x){ 

    /* Boundary Value Problem 1 */
    // VectorXm<mpreal> dxdt(2);
    // dxdt(0) = x(1);
    // dxdt(1) = -fabs(x(0));

    /* Boundary Value Problem 2 */
    // VectorXm<mpreal> dxdt(3);
    // dxdt(0) = x(1);
    // dxdt(1) = x(2);
    // dxdt(2) = 25*x(1) - 1;

    /* Boundary Value Problem 3 */
    VectorXm<mpreal> dxdt(2);
    dxdt(0) =  x(1);
    dxdt(1) = -x(0);

    return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at nodal state vectors xBC         -- (nx1) 
VectorXm<mpreal> BCResidues(MatrixXm<mpreal> xBC){
    
    /* Boundary Value Problem 1 */
    // VectorXm<mpreal> residues(2);
    // residues(0) = xBC(0,0) - 0;
    // residues(1) = xBC(0,1) + 2;

    /* Boundary Value Problem 2 */
    VectorXm<mpreal> residues(3);
    residues(0) = xBC(1,0) - 0;
    residues(1) = xBC(1,2) - 0;
    residues(2) = xBC(0,1) - 0;

    return residues;
}
// ===================

// ====================
// Functions BCResidues
// ====================
// BCResidues     = a function that defines the boundary condition residues at the nodal state vectors...  -- (n(m-1)x1)
//                  ... on the left and right side of integration intervals, xBCL and xBCR
VectorXm<mpreal> BCResidues(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){

    /* Boundary Value Problem 3 */
    VectorXm<mpreal> residues(10);
    residues(0) = xBCL(0,0)-0;
    residues(1) = xBCL(1,0)-1;
    residues(2) = xBCR(0,0) - xBCL(0,1) - 1;
    residues(3) = xBCR(1,0) - xBCL(1,1) - 0;
    residues(4) = xBCR(1,1) - xBCL(1,2) + 1;
    residues(5) = xBCR(0,1) - xBCL(0,2) + 0;
    residues(6) = xBCR(0,2) - xBCL(0,3) - 1;
    residues(7) = xBCR(1,2) - xBCL(1,3) - (sqrt(3)-1);
    residues(8) = xBCR(0,3) - xBCL(0,4) - 0;
    residues(9) = xBCR(1,3) - xBCL(1,4) - 0;

    return residues;
}


// =================
// The main function
// =================
// This is where the program execution begins
int main(
    int argc,   // argc = the number of arguments passed to the program
    char **argv // argv = an array of strings that are passed to  the program 
    ){

    mpreal::set_default_prec(64); // Set the number of digits of precision you want for computations

    cout<<endl;
    cout<<"===================================================================="<<endl;
    cout<<"Test: Non-linear multipoint boundary value problem solver (nlmpBVP)"<<endl;
    cout<<"===================================================================="<<endl;
    cout<<"Copyright Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)"<<endl;

    // Variable declarations   

    /* Boundary Value Problem 1 */
    // RowVectorXm<mpreal> tBC(2);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    // VectorXm<mpreal>   oxt1(2);            // oxt1           = column vector of the guessed initial state                                       -- (nx1)

    /* Boundary Value Problem 2 */
    // RowVectorXm<mpreal> tBC(3);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    // VectorXm<mpreal>   oxt1(3);            // oxt1           = column vector of the guessed initial state                                       -- (nx1)

    /* Boundary Value Problem 3 */
    RowVectorXm<mpreal> tBC(6);
    MatrixXm<mpreal>   oxt1(2,5); 

    BVPSolution<mpreal> bvpSolution;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    // Variable definitions

    /* Boundary Value Problem 1 */
    // tBC  << 0.0, 4.0;                     // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    // oxt1 <<  1,                           // oxt1 = column vector of the guessed initial state                                        -- (nx1) 
            //  0;  
    // tBC  << 0.0, 4.0;                     // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    // oxt1 << -1,                           // oxt1 = column vector of the guessed initial state                                        -- (nx1) 
    //          0;

    /* Boundary Value Problem 2 */
    // tBC  << 0.0, 0.5, 1.0;                // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    // oxt1 << 1,                            // oxt1 = column vector of the guessed initial state                                        -- (nx1)
    //         1,
    //         1;

    /* Boundary Value Problem 3 */
    // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    tBC  << 0.0, mpfr::const_pi()/6, mpfr::const_pi()/3, mpfr::const_pi()/2, 2*mpfr::const_pi()/3, mpfr::const_pi();
    // oxt1 = a matrix of the guessed initial state on the left side of each integration interval                                        -- (nx(m-1))
    oxt1 <<  0.1, 0.1, 0.4, 0.8, 0.9,
            -0.6, 0.1, 0.9, 2.1, 0.8;

    // Assign the parameters for IVAM
    ivamParameters.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1<<endl;

    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    
    // Solve the boundary value problem

    /* Boundary Value Problem 1 */
    // bvpSolution = nlmpBVP<mpreal>(2, 2, 101, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

    /* Boundary Value Problem 2 */
    // bvpSolution = nlmpBVP<mpreal>(3, 3, 101, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

    /* Boundary Value Problem 3 */
    bvpSolution = nlmpBVP2<mpreal>(2, 6, 12*10+1, tBC, oxt1, dxBydt, BCResidues, ivamParameters);

    cout<<endl<<"Done solving the BVP..."<<endl;

    // Print the output
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution.xBC<<endl;
    cout<<endl;    

    return 0;
}
// =================