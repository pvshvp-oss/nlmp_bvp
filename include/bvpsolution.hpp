// ========================================
// Author: Shivanand Pattanshetti 
// Email : shivanand.pattanshetti@gmail.com
// ========================================
// Copyright Shivanand Pattanshetti (shivanand.pattanshetti@gmail.com)

// ===============================
// Includes and global definitions
// ===============================
#include <Eigen/Eigen>                                              // For matrix and vector math
using namespace Eigen;                                              // For matrix and vector data types and operations
template <typename T> using MatrixXm = Matrix<T, Dynamic, Dynamic>; // Dynamic size matrix with the 'T' data-type
template <typename T> using VectorXm = Matrix<T, Dynamic, 1>;       // Dynamic size vector with the 'T' data-type
template <typename T> using RowVectorXm = Matrix<T, 1, Dynamic>;    // Dynamic size row vector with the 'T' data-type
using RowVectorXi = Matrix<int, 1, Dynamic>;                        // Dynamic size row vector with the 'integer' data-type
// ===============================

// =====================
// Structure BVPSolution
// =====================
template <typename T> struct BVPSolution{
    RowVectorXm<T>   t; // t    = row vector of values at which state vectors x are evaluated                      -- (1xnGrid)
    MatrixXm<T>      x; // x    = state vectors at all grid points t                                               -- (nxnGrid)  
    RowVectorXm<T> tBC; // tBC  = row vector of boundary condition values at which state vectors xBC are evaluated -- (1xm)   
    MatrixXm<T>    xBC; // xBC  = state vectors at boundary conditions tBC                                         -- (nxm)
};
// =====================