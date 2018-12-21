// ===============================
// Includes and global definitions
// ===============================
#include <cmath>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// C++ analog to math.h
#include <Eigen/Eigen>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// For matrix and vector math
#include <boost/numeric/odeint.hpp>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For the initial value problem (IVP) solver
#include <nlmp_bvp.hpp>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For the function declarations***REMOVED***
using namespace std;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  //
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***//
using namespace boost::numeric::odeint;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** //
using RowVectorXi = Matrix<int, 1, Dynamic>;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For the convenience of declaring row vectors
using RowVectorXd = Matrix<double, 1, Dynamic>;***REMOVED******REMOVED******REMOVED******REMOVED***  // For the convenience of declaring row vectors
using StepperType = runge_kutta_dopri5<***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For the convenience of specifying the stepper type for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** VectorXd,***REMOVED******REMOVED******REMOVED******REMOVED***// the state vector type for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** double,***REMOVED******REMOVED******REMOVED******REMOVED***  // the state variable value type for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** VectorXd,***REMOVED******REMOVED******REMOVED******REMOVED***// the type for the derivative of the state vector x for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** double,***REMOVED******REMOVED******REMOVED******REMOVED***  // the type for independent variable t for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** vector_space_algebra // the type of algebra to be done for the IVP solver
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** >;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // 
const double INF = std::numeric_limits<double>::infinity();  // INF = infinity***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
// ===============================

// ==================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state is evaluated
***REMOVED*** RowVectorXd t_BC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXd _0_x_t1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // _0_x_t1***REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXd dxBydt(double t, VectorXd x), // dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXd BCResidues(MatrixXd x_BC),***REMOVED*** // BCResidues***REMOVED***  = a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
***REMOVED*** IVAMParameters ivamParameters***REMOVED******REMOVED******REMOVED*** // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** ){  

***REMOVED******REMOVED***  // Variable declarations***REMOVED******REMOVED***  
***REMOVED******REMOVED***  int j;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // j***REMOVED******REMOVED******REMOVED******REMOVED***= the inner iterating variable for IVAM***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- [0,n-1]
***REMOVED******REMOVED***  int k;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // k***REMOVED******REMOVED******REMOVED******REMOVED***= the outer iterating variable for IVAM***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- [0,Inf)
***REMOVED******REMOVED***  int iCol;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // iCol***REMOVED******REMOVED******REMOVED***= the column index of the x solution for the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- [0,nGrid-1]
***REMOVED******REMOVED***  int iColP;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// iColP***REMOVED******REMOVED***  = the column index of the perturbed x solution for the IVP solver***REMOVED***
***REMOVED******REMOVED***  double h;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // h***REMOVED******REMOVED******REMOVED******REMOVED***= the stepper step size for the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- [0,nGrid-1] 
***REMOVED******REMOVED***  double t0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// t0***REMOVED******REMOVED******REMOVED***  = the first boundary value of the independent variable
***REMOVED******REMOVED***  double tm;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tm***REMOVED******REMOVED******REMOVED***  = the last boundary value of the independent variable
***REMOVED******REMOVED***  double _k_epsilon_j;***REMOVED******REMOVED***  // _k_epsilon_j = the perturbation parameter for a state variable at a particular iteration k
***REMOVED******REMOVED***  double _k_alpha;***REMOVED******REMOVED******REMOVED******REMOVED***// _k_alpha***REMOVED***  = the relaxation factor at a particular iteration k to scale the adjustment to the initial condition 
***REMOVED******REMOVED***  double _k_G;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // _k_G***REMOVED******REMOVED******REMOVED***= the Root Mean Square (RMS) error of boundary residues at a particular iteration k
***REMOVED******REMOVED***  double _k_GPrev;***REMOVED******REMOVED******REMOVED******REMOVED***// _k_GPrev***REMOVED***  = the Root Mean Square (RMS) error of boundary residues at the previous iteration k-1
***REMOVED******REMOVED***  RowVectorXd tSol(nGrid);***REMOVED*** // tSol***REMOVED******REMOVED******REMOVED***= the independent variable t over the whole grid in the solution of the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (1xnGrid)
***REMOVED******REMOVED***  MatrixXd xSol(n,nGrid);***REMOVED***  // xSol***REMOVED******REMOVED******REMOVED***= the state vector x integrated over the whole grid in the solution of the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nxnGrid)
***REMOVED******REMOVED***  RowVectorXd tsolP(nGrid);***REMOVED***// tsolP***REMOVED******REMOVED***  = the independent variable t over the whole grid in the perturbed solution of the IVP solver***REMOVED******REMOVED******REMOVED***-- (1xnGrid)***REMOVED*** 
***REMOVED******REMOVED***  MatrixXd xSolP(n,nGrid);***REMOVED*** // xSolP***REMOVED******REMOVED***  = the state vector x integrated over the whole grid in the perturbed solution of the IVP solver***REMOVED******REMOVED***-- (nxnGrid)
***REMOVED******REMOVED***  RowVectorXi BCCols(m);***REMOVED******REMOVED***// BCCols***REMOVED******REMOVED*** = the columns in the grid that correspond to boundary values***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED******REMOVED***  VectorXd _k_x_t1(n);***REMOVED******REMOVED***  // _k_x_t1***REMOVED******REMOVED***= the computed initial state vector in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED******REMOVED***  VectorXd _k_x_t1P(n);***REMOVED******REMOVED*** // _k_x_t1***REMOVED******REMOVED***= the computed perturbed initial state vector in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED***  VectorXd x_t1(n);***REMOVED******REMOVED******REMOVED***  // x_t1***REMOVED******REMOVED******REMOVED***= the computed initial state vector to be input to the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED******REMOVED***  VectorXd x_t1P(n);***REMOVED******REMOVED******REMOVED*** // x_t1P***REMOVED******REMOVED***  = the computed perturbed initial state vector to be input to the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED******REMOVED***  VectorXd _k_g(n);***REMOVED******REMOVED******REMOVED***  // _k_g***REMOVED******REMOVED******REMOVED***= the boundary condition residues in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED***  VectorXd _k_g_j(n);***REMOVED******REMOVED******REMOVED***// _k_g_j***REMOVED******REMOVED*** = the j-th boundary condition perturbed system residues in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED******REMOVED***  MatrixXd _k_S(n,n);***REMOVED******REMOVED******REMOVED***// _k_S***REMOVED******REMOVED******REMOVED***= the adjusting matrix for correcting the initial condition k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nxn) 
***REMOVED******REMOVED***  BVPSolution bvpSolution;

***REMOVED******REMOVED***  // Variable definitions***REMOVED******REMOVED***  
***REMOVED******REMOVED***  h***REMOVED******REMOVED***= (tm - t0)/(nGrid-1);
***REMOVED******REMOVED***  t0***REMOVED***  = t_BC(0);
***REMOVED******REMOVED***  tm***REMOVED***  = t_BC(m-1);
***REMOVED******REMOVED***  BCCols = ((t_BC-t0*RowVectorXd::Ones(m))/h).array().round().cast<int>();

***REMOVED******REMOVED***  // Wrapper function to be called by the IVP solver to retrieve the definitions for the differential equations
***REMOVED******REMOVED***  auto dxBydtWrapper = [dxBydt] // Captured variables
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  (const VectorXd &x, VectorXd &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dxBydt(t, x);
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  // Observer to store the solutions of the IVP Solver
***REMOVED******REMOVED***  auto storeSol = [&iCol, &tSol, &xSol] // Captured variables
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***(const VectorXd &x , const double t){
***REMOVED******REMOVED******REMOVED******REMOVED***tSol(iCol)***REMOVED***  = t;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***xSol.col(iCol) = x;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***++iCol;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  // Observer to store the perturbed solutions of the IVP Solver
***REMOVED******REMOVED***  auto storeSolP = [&iColP, &tsolP, &xSolP] // Captured variables
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** (const VectorXd &x , const double t){
***REMOVED******REMOVED******REMOVED******REMOVED***tsolP(iColP)***REMOVED***  = t;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***xSolP.col(iColP) = x;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***++iColP;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  };***REMOVED******REMOVED*** 
 
***REMOVED******REMOVED***  _k_GPrev = INF;
***REMOVED******REMOVED***  _k_alpha = ivamParameters.ALPHA;
***REMOVED******REMOVED***  k***REMOVED******REMOVED***  = 0;***REMOVED******REMOVED*** // Set k for the first time
***REMOVED******REMOVED***  _k_x_t1  = _0_x_t1; // Assign the initial condition state vector for the 0-th iteration

***REMOVED******REMOVED***  // Solve the initial value problem for the first time 
***REMOVED******REMOVED***  x_t1***REMOVED***  = _k_x_t1; // Assign the current initial condition state vector to a dummy variable
***REMOVED******REMOVED***  iCol***REMOVED***  = 0;***REMOVED******REMOVED*** // Set the solution column index to 0 before the IVP solver starts integrating
***REMOVED******REMOVED***  integrate_const(StepperType(), dxBydtWrapper, x_t1, t0, tm, h, storeSol); 
***REMOVED******REMOVED***  _k_g = BCResidues(xSol(Eigen::all, BCCols));
***REMOVED******REMOVED***  _k_G = _k_g.norm()/sqrt(n);
***REMOVED******REMOVED***  while(_k_G > ivamParameters.SIGMA){  
***REMOVED******REMOVED******REMOVED******REMOVED***if(_k_G < 0.1*_k_GPrev) {
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** _k_alpha = fmin(1.2*_k_alpha, 1.0);***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED***} else if(_k_G >= _k_GPrev){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** _k_alpha = 0.8*_k_alpha;
***REMOVED******REMOVED******REMOVED******REMOVED***}***REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***for(j = 0; j < n; j++){***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Determine the perturbation parameter
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** _k_epsilon_j = fmax(ivamParameters.EPSILON, fabs(ivamParameters.EPSILON * _k_x_t1(j)));
***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Perturb the initial conditions***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** _k_x_t1P = _k_x_t1 + _k_epsilon_j*MatrixXd::Identity(n,n).col(j);
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Solve the perturbed initial value problem***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** x_t1P = _k_x_t1P; // Assign the current initial condition state vector to a dummy variable***REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** iColP  = 0;***REMOVED******REMOVED***  // Set the solution column index to 0 before the IVP solver starts integrating***REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** integrate_const(StepperType(), dxBydtWrapper, x_t1P, t0, tm, h, storeSolP); 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** _k_g_j =  BCResidues(xSolP(Eigen::all, BCCols));

***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Compute a column of the adjusting matrix***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** _k_S.col(j) = (_k_g_j- _k_g)/_k_epsilon_j;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED***}

***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the linarized adjusting equation
***REMOVED******REMOVED******REMOVED******REMOVED***_k_x_t1 = _k_S.colPivHouseholderQr().solve(-_k_alpha*_k_g) + _k_x_t1;

***REMOVED******REMOVED******REMOVED******REMOVED***_k_GPrev = _k_G;
***REMOVED******REMOVED******REMOVED******REMOVED***++k;

***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the initial value problem***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***x_t1***REMOVED***  = _k_x_t1; // Assign the current initial condition state vector to a dummy variable
***REMOVED******REMOVED******REMOVED******REMOVED***iCol***REMOVED***  = 0;***REMOVED******REMOVED*** // Set the solution column index to 0 before the IVP solver starts integrating
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(StepperType(), dxBydtWrapper, x_t1, t0, tm, h, storeSol); 
***REMOVED******REMOVED******REMOVED******REMOVED***_k_g = BCResidues(xSol(Eigen::all, BCCols));
***REMOVED******REMOVED******REMOVED******REMOVED***_k_G = _k_g.norm()/sqrt(n);***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  }  
***REMOVED******REMOVED***  bvpSolution.t***REMOVED*** = tSol;
***REMOVED******REMOVED***  bvpSolution.x***REMOVED*** = xSol;
***REMOVED******REMOVED***  bvpSolution.t_BC = t_BC;
***REMOVED******REMOVED***  bvpSolution.x_BC = xSol(Eigen::all, BCCols);
***REMOVED******REMOVED***  return bvpSolution;
***REMOVED*** }
// ===================


