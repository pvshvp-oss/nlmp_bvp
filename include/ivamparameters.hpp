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

// ========================
// Structure IVAMParameters
// ========================
// Parameters for the Initial Value Adjusting Method (IVAM)
template <typename T> struct IVAMParameters{
    T EPSILON;        // EPSILON    = the state perturbation parameter to probe the differential equation system with
    T ALPHA;          // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    T SIGMA;          // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    T BETA;           // BETA       = the deflation factor
    bool printDebug;  // printDebug = specify whether debug messages should be output to the console
};
// ========================