// ========================================
// Author: Shivanand Pattanshetti 
// Email : shivanand.pattanshetti@gmail.com
// ========================================

// ========================
// Structure IVAMParameters
// ========================
// Parameters for the Initial Value Adjusting Method (IVAM)
template <typename T> struct IVAMParameters{
    T EPSILON; // EPSILON = the state perturbation parameter to probe the differential equation system with
    T ALPHA;   // ALPHA   = the relaxation factor to scale the adjustment to the initial condition
    T SIGMA;   // SIGMA   = the tolerance for error outside which the solver needs to  iterate further. 
    T BETA;    // BETA    = the deflation factor
};
// ========================