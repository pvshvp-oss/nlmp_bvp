// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================

// ===============================
// Includes and global definitions
// ===============================
#include <cmath>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // C++ analog to math.h
#include <boost/numeric/odeint.hpp>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For the initial value problem (IVP) solver
#include "ivamparameters.hpp"***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For the struct "IVAMParameters"
#include "bvpsolution.hpp"***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For the struct "BVPSolution"
using namespace std;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For cout
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For matrix and vector data types and operations
using namespace boost::numeric::odeint;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// For the initial value problem (IVP) solver***REMOVED***  
template <typename T> using StepperType = runge_kutta_dopri5<***REMOVED******REMOVED******REMOVED***  // For the convenience of specifying the stepper type for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** VectorXm<T>,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // the state vector type for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** T,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// the state variable value type for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** VectorXm<T>,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // the type for the derivative of the state vector x for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** T,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// the type for independent variable t for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** vector_space_algebra***REMOVED******REMOVED******REMOVED******REMOVED***// the type of algebra to be done for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** >;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// 
template <typename T> const T INF = std::numeric_limits<T>::infinity(); // INF = infinity***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
// ===============================

// ==================
// Function "nlmpBVP"
// ==================
template <typename T> 
BVPSolution<T> nlmpBVP(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state is evaluated
***REMOVED*** RowVectorXm<T> tBC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // tBC***REMOVED******REMOVED******REMOVED******REMOVED***= row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm<T> oxt1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXm<T> dxBydt(T t, VectorXm<T> x),  // dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXm<T> BCResidues(MatrixXm<T> xBC), // BCResidues***REMOVED***  = a function that defines the boundary condition residues at state vectors xBC  -- (nx1) 
***REMOVED*** IVAMParameters<T> ivamParameters***REMOVED******REMOVED******REMOVED***// ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** ){  

***REMOVED******REMOVED***  // Variable declarations  
***REMOVED******REMOVED***  int omega = 0.5; 
***REMOVED******REMOVED***  int j;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // j***REMOVED******REMOVED******REMOVED***= the inner iterating variable for IVAM***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- [0,n-1]
***REMOVED******REMOVED***  int k;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // k***REMOVED******REMOVED******REMOVED***= the outer iterating variable for IVAM***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- [0,Inf)
***REMOVED******REMOVED***  int iCol;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // iCol***REMOVED******REMOVED***= the column index of the x solution for the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- [0,nGrid-1]
***REMOVED******REMOVED***  int iColP;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // iColP***REMOVED***  = the column index of the perturbed x solution for the IVP solver***REMOVED***
***REMOVED******REMOVED***  T h;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // h***REMOVED******REMOVED******REMOVED***= the stepper step size for the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- [0,nGrid-1] 
***REMOVED******REMOVED***  T t0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // t0***REMOVED******REMOVED***  = the first boundary value of the independent variable
***REMOVED******REMOVED***  T tm;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // tm***REMOVED******REMOVED***  = the last boundary value of the independent variable
***REMOVED******REMOVED***  T kepsilonj;***REMOVED******REMOVED******REMOVED******REMOVED***// kepsilonj = the perturbation parameter for a state variable at a particular iteration k
***REMOVED******REMOVED***  T kalpha;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// kalpha***REMOVED*** = the relaxation factor at a particular iteration k to scale the adjustment to the initial condition 
***REMOVED******REMOVED***  T kG;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // kG***REMOVED******REMOVED***  = the Root Mean Square (RMS) error of boundary residues at a particular iteration k
***REMOVED******REMOVED***  T kGPrev;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// kGPrev***REMOVED*** = the Root Mean Square (RMS) error of boundary residues at the previous iteration k-1
***REMOVED******REMOVED***  RowVectorXm<T> tSol(nGrid);***REMOVED***  // tSol***REMOVED******REMOVED***= the independent variable t over the whole grid in the solution of the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (1xnGrid)
***REMOVED******REMOVED***  MatrixXm<T> xSol(n,nGrid);***REMOVED******REMOVED***// xSol***REMOVED******REMOVED***= the state vector x integrated over the whole grid in the solution of the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nxnGrid)
***REMOVED******REMOVED***  RowVectorXm<T> tSolPert(nGrid); // tSolPert  = the independent variable t over the whole grid in the perturbed solution of the IVP solver***REMOVED******REMOVED******REMOVED***-- (1xnGrid)***REMOVED*** 
***REMOVED******REMOVED***  MatrixXm<T> xSolPert(n,nGrid);  // xSolPert  = the state vector x integrated over the whole grid in the perturbed solution of the IVP solver***REMOVED******REMOVED***-- (nxnGrid)
***REMOVED******REMOVED***  RowVectorXi BCCols(m);***REMOVED******REMOVED*** // BCCols***REMOVED*** = the columns in the grid that correspond to boundary values***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED******REMOVED***  VectorXm<T> kxt1(n);***REMOVED******REMOVED******REMOVED******REMOVED***// kxt1***REMOVED******REMOVED***= the computed initial state vector in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED******REMOVED***  VectorXm<T> kxt1Prev(n);***REMOVED******REMOVED***  // kxt1Prev  = the computed initial state vector in the previous (k-1)-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED***  VectorXm<T> kxt1P(n);***REMOVED******REMOVED******REMOVED***  // kxt1***REMOVED******REMOVED***= the computed perturbed initial state vector in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED***  VectorXm<T> xt1(n);***REMOVED******REMOVED******REMOVED******REMOVED*** // xt1***REMOVED******REMOVED*** = the computed initial state vector to be input to the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED******REMOVED***  VectorXm<T> xt1P(n);***REMOVED******REMOVED******REMOVED******REMOVED***// xt1P***REMOVED******REMOVED***= the computed perturbed initial state vector to be input to the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED******REMOVED***  VectorXm<T> kg(n);***REMOVED******REMOVED******REMOVED******REMOVED***  // kg***REMOVED******REMOVED***  = the boundary condition residues in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED***  VectorXm<T> kgj(n);***REMOVED******REMOVED******REMOVED******REMOVED*** // kgj***REMOVED******REMOVED*** = the j-th boundary condition perturbed system residues in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED******REMOVED***  MatrixXm<T> kS(n,n);***REMOVED******REMOVED******REMOVED******REMOVED***// kS***REMOVED******REMOVED***  = the adjusting matrix for correcting the initial condition k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nxn) 
***REMOVED******REMOVED***  BVPSolution<T> bvpSolution;

***REMOVED******REMOVED***  // Variable definitions
***REMOVED******REMOVED***  t0***REMOVED***  = tBC(0);***REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED***  tm***REMOVED***  = tBC(m-1);
***REMOVED******REMOVED***  h***REMOVED******REMOVED***= (tm - t0)/(nGrid-1);
***REMOVED******REMOVED***  BCCols = ((tBC-t0*RowVectorXm<T>::Ones(m))/h).template array().template round().template cast<int>();

***REMOVED******REMOVED***  cout<<"Boundary nodes correspond to the below columns: "<<endl<<BCCols<<endl<<endl;
***REMOVED******REMOVED***  
***REMOVED******REMOVED***  // Wrapper function to be called by the IVP solver to retrieve the definitions for the differential equations
***REMOVED******REMOVED***  auto dxBydtWrapper = [dxBydt] // Captured variables
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  (const VectorXm<T> &x, VectorXm<T> &dxdt, T t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dxBydt(t, x);
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  // Observer to store the solutions of the IVP Solver
***REMOVED******REMOVED***  auto storeSol = [&iCol, &tSol, &xSol] // Captured variables
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***(const VectorXm<T> &x , const T t){
***REMOVED******REMOVED******REMOVED******REMOVED***tSol(iCol)***REMOVED***  = t;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***xSol.col(iCol) = x;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***++iCol;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  // Observer to store the perturbed solutions of the IVP Solver
***REMOVED******REMOVED***  auto storeSolP = [&iColP, &tSolPert, &xSolPert] // Captured variables
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** (const VectorXm<T> &x , const T t){
***REMOVED******REMOVED******REMOVED******REMOVED***tSolPert(iColP)***REMOVED***  = t;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***xSolPert.col(iColP) = x;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***++iColP;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  };***REMOVED******REMOVED*** 
 
***REMOVED******REMOVED***  // Initialize the required variables
***REMOVED******REMOVED***  kGPrev***REMOVED***= INF<T>;
***REMOVED******REMOVED***  kalpha***REMOVED***= ivamParameters.ALPHA;
***REMOVED******REMOVED***  k***REMOVED******REMOVED***  = 0;***REMOVED******REMOVED*** 
***REMOVED******REMOVED***  kxt1***REMOVED***  = oxt1; 
***REMOVED******REMOVED***  kxt1Prev = kxt1;

***REMOVED******REMOVED***  // Solve the initial value problem for the first time 
***REMOVED******REMOVED***  xt1***REMOVED******REMOVED***= kxt1; // Assign the current initial condition state vector to a dummy variable
***REMOVED******REMOVED***  iCol***REMOVED***  = 0;***REMOVED*** // Set the solution column index to 0 before the IVP solver starts integrating
***REMOVED******REMOVED***  integrate_const(StepperType<T>(), dxBydtWrapper, xt1, t0, tm, h, storeSol); 
***REMOVED******REMOVED***  kg = BCResidues(xSol(Eigen::all, BCCols));
***REMOVED******REMOVED***  kG = kg.norm()/sqrt(n);

***REMOVED******REMOVED***  while(kG > ivamParameters.SIGMA /* When the error is more than the max. threshold */){***REMOVED******REMOVED***

***REMOVED******REMOVED******REMOVED******REMOVED***// Adjust the relaxation parameter to control the rate of convergence
***REMOVED******REMOVED******REMOVED******REMOVED***if(kG < 0.01*kGPrev) {
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Increasing alpha..."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kalpha = fmin(1.2*kalpha, 1.0); 
***REMOVED******REMOVED******REMOVED******REMOVED***} else if(kG >= kGPrev){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Decreasing alpha..."<<endl;***REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kalpha = 0.8*kalpha;***REMOVED******REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***} else{
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***}

***REMOVED******REMOVED******REMOVED******REMOVED***// Inner loop to perturb each state variable separately and find the normalized change in the boundary condition residues
***REMOVED******REMOVED******REMOVED******REMOVED***for(j = 0; j < n; j++){***REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Determine the perturbation parameter
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kepsilonj = fmax(ivamParameters.EPSILON, fabs(ivamParameters.EPSILON * kxt1(j)));
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // kepsilonj = ivamParameters.EPSILON;
***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Perturb the initial conditions***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kxt1P = kxt1 + kepsilonj*MatrixXm<T>::Identity(n,n).col(j);
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Solve the perturbed initial value problem***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** xt1P = kxt1P; // Assign the current initial condition state vector to a dummy variable***REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** iColP  = 0;***REMOVED***// Set the solution column index to 0 before the IVP solver starts integrating***REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** integrate_const(StepperType<T>(), dxBydtWrapper, xt1P, t0, tm, h, storeSolP); 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kgj =  BCResidues(xSolPert(Eigen::all, BCCols));

***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Compute one column of the adjusting matrix
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Each column corresponds to the change in boundary condition residues due to a corresponding perturbation in *one* state variable, normalized with respect to the perturbation
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kS.col(j) = (kgj - kg)*pow(kepsilonj,-1);  
***REMOVED******REMOVED******REMOVED******REMOVED***}

***REMOVED******REMOVED******REMOVED******REMOVED***kalpha = 1;
***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the linarized adjusting equation***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***kxt1Prev = kxt1; 
***REMOVED******REMOVED******REMOVED******REMOVED***kxt1 = kxt1 - kS.colPivHouseholderQr().solve(kalpha*kg);

***REMOVED******REMOVED******REMOVED******REMOVED***// Start the next iteration
***REMOVED******REMOVED******REMOVED******REMOVED***kGPrev = kG;
***REMOVED******REMOVED******REMOVED******REMOVED***++k;

***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the initial value problem***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***xt1***REMOVED***  = kxt1; // Assign the current initial condition state vector to a dummy variable
***REMOVED******REMOVED******REMOVED******REMOVED***iCol***REMOVED***  = 0;***REMOVED***// Set the solution column index to 0 before the IVP solver starts integrating
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(StepperType<T>(), dxBydtWrapper, xt1, t0, tm, h, storeSol); 
***REMOVED******REMOVED******REMOVED******REMOVED***kg = BCResidues(xSol(Eigen::all, BCCols));
***REMOVED******REMOVED******REMOVED******REMOVED***kG = kg.norm()/sqrt(n);

***REMOVED******REMOVED******REMOVED******REMOVED***if(k >= 1000){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"[WARNING]: The solution did not converge after 1000 iterations. Terminating the process."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** break;
***REMOVED******REMOVED******REMOVED******REMOVED***}
***REMOVED******REMOVED***  }  

***REMOVED******REMOVED***  cout<<"Ran "<<k<<" iteration(s)."<<endl;

***REMOVED******REMOVED***  bvpSolution.t***REMOVED***= tSol;
***REMOVED******REMOVED***  bvpSolution.x***REMOVED***= xSol;
***REMOVED******REMOVED***  bvpSolution.tBC = tBC;
***REMOVED******REMOVED***  bvpSolution.xBC = xSol(Eigen::all, BCCols);
***REMOVED******REMOVED***  return bvpSolution;
***REMOVED*** }
// ==================