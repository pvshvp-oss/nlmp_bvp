// ========================================
// Author: Shivanand Pattanshetti 
// Email : shivanand.pattanshetti@gmail.com
// ========================================
// Copyright Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)

// ==========
// References
// ==========
// [1] Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."
//     Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.
// [2] Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
//     Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

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
// dxBydt = a function that defines the derivative of a state vector x at t -- (nx1)

/* Boundary Value Problem 1 */
VectorXm<mpreal> dxBydt_BVP1(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(2);
    dxdt(0) = x(1);
    dxdt(1) = -fabs(x(0));
    return dxdt;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> dxBydt_BVP2(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(3);
    dxdt(0) = x(1);
    dxdt(1) = x(2);
    dxdt(2) = 25*x(1) - 1;
    return dxdt;
}

/* Boundary Value Problem 3 */
VectorXm<mpreal> dxBydt_BVP3(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(2);
    dxdt(0) =  x(1);
    dxdt(1) = -x(0);
    return dxdt;
}

/* Boundary Value Problem 4 */
VectorXm<mpreal> dxBydt_BVP4(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(2);
    dxdt(0) =  x(1);
    dxdt(1) = -x(0);
    return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at nodal state vectors xBC -- (nx1) 

/* Boundary Value Problem 1 */
VectorXm<mpreal> BCResidues_BVP1(MatrixXm<mpreal> xBC){
    VectorXm<mpreal> residues(2);
    residues(0) = xBC(0,0) - 0;
    residues(1) = xBC(0,1) + 2;
    return residues;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> BCResidues_BVP2(MatrixXm<mpreal> xBC){
    VectorXm<mpreal> residues(3);
    residues(0) = xBC(1,0) - 0;
    residues(1) = xBC(1,2) - 0;
    residues(2) = xBC(0,1) - 0;
    return residues;
}


// BCResidues     = a function that defines the boundary condition residues at the nodal state vectors... -- (n(m-1)x1)
//                  ... on the left and right side of integration intervals, xBCL and xBCR

/* Boundary Value Problem 3 */
VectorXm<mpreal> BCResidues_BVP3(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
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

/* Boundary Value Problem 4 */
VectorXm<mpreal> BCResidues_BVP4(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
    VectorXm<mpreal> residues(8);
    residues(0) = xBCL(0,0)-0;
    residues(1) = xBCL(1,0)-1;
    residues(2) = pow(xBCL(0,1),2)*pow(xBCL(1,1),3) - exp(xBCR(1,1))*pow(xBCR(1,0),2)*pow(xBCL(1,3),2) + pow(xBCR(1,3),2)*pow(xBCR(0,2),2) - 0.1859767072;
    residues(3) = pow(xBCR(0,1),2)*pow(xBCR(0,0),3)*pow(xBCL(0,3),2) + pow(xBCL(1,2),2)*exp(xBCR(1,2)) + pow(xBCL(0,2),2)*pow(xBCR(0,3),2) - 0.1261677772;
    residues(4) = xBCR(0,0) - xBCL(0,1) + 0;
    residues(5) = xBCR(1,0) - xBCL(1,1) + 0;
    residues(6) = xBCR(0,2) - xBCL(0,3) + 0;
    residues(7) = xBCR(1,2) - xBCL(1,3) + 0;
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
    cout<<"============================================================================="<<endl;
    cout<<"Test: Non-linear multipoint boundary value problem solver (nlmpBVP, nlmpBVP2)"<<endl;
    cout<<"============================================================================="<<endl;
    cout<<"Copyright Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)"<<endl;

    // Variable declarations   

    /* Boundary Value Problem 1 */
    RowVectorXm<mpreal> tBC_BVP1(2);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    VectorXm<mpreal>   oxt1_BVP1(2);            // oxt1           = column vector of the guessed initial state                                       -- (nx1)

    /* Boundary Value Problem 2 */
    RowVectorXm<mpreal> tBC_BVP2(3);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    VectorXm<mpreal>   oxt1_BVP2(3);            // oxt1           = column vector of the guessed initial state                                       -- (nx1)

    /* Boundary Value Problem 3 */
    RowVectorXm<mpreal> tBC_BVP3(6);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)  
    MatrixXm<mpreal> oxt1_BVP3(2,5);            // oxt1           = matrix of the guessed initial state                                              -- (nx(m-1))

    /* Boundary Value Problem 4 */
    RowVectorXm<mpreal> tBC_BVP4(5);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    MatrixXm<mpreal> oxt1_BVP4(2,4);            // oxt1           = matrix of the guessed initial state                                              -- (nx(m-1))

    BVPSolution<mpreal> bvpSolution_BVP1;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP1; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    BVPSolution<mpreal> bvpSolution_BVP2;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP2; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    BVPSolution<mpreal> bvpSolution_BVP3;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP3; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    BVPSolution<mpreal> bvpSolution_BVP4;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP4; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    // Variable definitions

    /* Boundary Value Problem 1 */
    tBC_BVP1  << 0.0, 4.0;                     // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    oxt1_BVP1 <<  1,                           // oxt1 = column vector of the guessed initial state                                        -- (nx1) 
             0;  
    // tBC_BVP1  << 0.0, 4.0;                     // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    // oxt1_BVP1 << -1,                           // oxt1 = column vector of the guessed initial state                                        -- (nx1) 
    //          0;

    /* Boundary Value Problem 2 */
    tBC_BVP2  << 0.0, 0.5, 1.0;                // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    oxt1_BVP2 << 1,                            // oxt1 = column vector of the guessed initial state                                        -- (nx1)
            1,
            1;

    /* Boundary Value Problem 3 */
    // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    tBC_BVP3  << 0.0, mpfr::const_pi()/6, mpfr::const_pi()/3, mpfr::const_pi()/2, 2*mpfr::const_pi()/3, mpfr::const_pi();
    // oxt1 = a matrix of the guessed initial state on the left side of each integration interval                                        -- (nx(m-1))
    oxt1_BVP3 <<  0.1, 0.1, 0.4, 0.8, 0.9,
            -0.6, 0.1, 0.9, 2.1, 0.8;

    /* Boundary Value Problem 4 */
    // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    tBC_BVP4  << 0.0, mpfr::const_pi()/4, mpfr::const_pi()/2, 3*mpfr::const_pi()/4, mpfr::const_pi();
    // oxt1 = a matrix of the guessed initial state on the left side of each integration interval                                        -- (nx(m-1))
    oxt1_BVP4 << -0.1, 0.7, 0.1,-0.7,
             1.1, 0.7,-1.1,-0.7;

    // Assign the parameters for IVAM

    ivamParameters_BVP1.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP1.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP1.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP1.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP1.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    ivamParameters_BVP2.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP2.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP2.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP2.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP2.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    ivamParameters_BVP3.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP3.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP3.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP3.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP3.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    ivamParameters_BVP4.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP4.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP4.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP4.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP4.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    /* Boundary Value Problem 1 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 1"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP1<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP1<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP1 = nlmpBVP<mpreal>(2, 2, 101, tBC_BVP1, oxt1_BVP1, dxBydt_BVP1, BCResidues_BVP1, ivamParameters_BVP1);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP1.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP1.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP1.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP1.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

     /* Boundary Value Problem 2 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 2"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP2<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP2<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP2 = nlmpBVP<mpreal>(3, 3, 101, tBC_BVP2, oxt1_BVP2, dxBydt_BVP2, BCResidues_BVP2, ivamParameters_BVP2);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP2.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP2.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP2.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP2.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

    /* Boundary Value Problem 3 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 3"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP3<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP3<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP3 = nlmpBVP2<mpreal>(2, 6, 12*10+1, tBC_BVP3, oxt1_BVP3, dxBydt_BVP3, BCResidues_BVP3, ivamParameters_BVP3);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP3.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP3.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP3.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP3.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

    /* Boundary Value Problem 4 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 4"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP4<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP4<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP4 = nlmpBVP2<mpreal>(2, 5, 10*10+1, tBC_BVP4, oxt1_BVP4, dxBydt_BVP4, BCResidues_BVP4, ivamParameters_BVP4);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP4.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP4.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP4.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP4.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

    cout<<endl<<"Program ended..."<<endl<<endl;

    return 0;
}
// =================