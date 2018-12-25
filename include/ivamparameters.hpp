// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================

// ========================
// Structure IVAMParameters
// ========================
// Parameters for the Initial Value Adjusting Method (IVAM)
template <typename T> struct IVAMParameters{
***REMOVED*** T EPSILON; // EPSILON = the state perturbation parameter to probe the differential equation system with
***REMOVED*** T ALPHA;***REMOVED***// ALPHA***REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** T SIGMA;***REMOVED***// SIGMA***REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** T BETA;***REMOVED*** // BETA***REMOVED*** = the deflation factor
};
// ========================