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

// ========================
// Structure IVAMParameters
// ========================
// Parameters for the Initial Value Adjusting Method (IVAM)
template <typename T> struct IVAMParameters{
***REMOVED*** T EPSILON;***REMOVED******REMOVED***  // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** T ALPHA;***REMOVED******REMOVED******REMOVED*** // ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** T SIGMA;***REMOVED******REMOVED******REMOVED*** // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** T BETA;***REMOVED******REMOVED******REMOVED***  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** bool printDebug;  // printDebug = specify whether debug messages should be output to the console
};
// ========================