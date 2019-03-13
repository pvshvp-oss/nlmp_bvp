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

// ===============================
// Includes and global definitions
// ===============================
#include <Eigen/Eigen>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For matrix and vector math
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For matrix and vector data types and operations
template <typename T> using MatrixXm = Matrix<T, Dynamic, Dynamic>; // Dynamic size matrix with the 'T' data-type
template <typename T> using VectorXm = Matrix<T, Dynamic, 1>;***REMOVED******REMOVED*** // Dynamic size vector with the 'T' data-type
template <typename T> using RowVectorXm = Matrix<T, 1, Dynamic>;***REMOVED*** // Dynamic size row vector with the 'T' data-type
using RowVectorXi = Matrix<int, 1, Dynamic>;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// Dynamic size row vector with the 'integer' data-type
// ===============================

// =====================
// Structure BVPSolution
// =====================
template <typename T> struct BVPSolution{
***REMOVED*** RowVectorXm<T>***REMOVED***t; // t***REMOVED*** = row vector of values at which state vectors x are evaluated***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (1xnGrid)
***REMOVED*** MatrixXm<T>***REMOVED******REMOVED***x; // x***REMOVED*** = state vectors at all grid points t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nxnGrid)  
***REMOVED*** RowVectorXm<T> tBC; // tBC  = row vector of boundary condition values at which state vectors xBC are evaluated -- (1xm)***REMOVED***
***REMOVED*** MatrixXm<T>***REMOVED*** xBC; // xBC  = state vectors at boundary conditions tBC***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nxm)
};
// =====================