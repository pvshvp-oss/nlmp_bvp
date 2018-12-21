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
// ===============================

// ===================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// n***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// m***REMOVED******REMOVED******REMOVED******REMOVED***  = the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // nGrid***REMOVED******REMOVED******REMOVED*** = the number of points at which the state is evaluated
***REMOVED*** RowVectorXd t_BC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXd _0_x_t1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // _0_x_t1***REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXd dxBydt(double t, VectorXd x), // dxBydt***REMOVED******REMOVED******REMOVED***= a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXd BCResidue(MatrixXd x_BC)***REMOVED******REMOVED***// BCResidue***REMOVED******REMOVED***= a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
***REMOVED*** const IVAMParameters ivamParameters***REMOVED*** // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
***REMOVED*** ){  

***REMOVED******REMOVED***  // Variable declarations***REMOVED******REMOVED***  
***REMOVED******REMOVED***  int j;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // j***REMOVED******REMOVED******REMOVED******REMOVED***= the inner iterating variable for IVAM***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- [0,n-1]
***REMOVED******REMOVED***  int k;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // k***REMOVED******REMOVED******REMOVED******REMOVED***= the outer iterating variable for IVAM***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- [0,Inf)
***REMOVED******REMOVED***  int iCol;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // iCol***REMOVED******REMOVED******REMOVED***= the column index of the x solution for the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- [0,nGrid-1]
***REMOVED******REMOVED***  int iColP;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // iColP***REMOVED******REMOVED***  = the column index of the perturbed x solution for the IVP solver***REMOVED***
***REMOVED******REMOVED***  double h;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // h***REMOVED******REMOVED******REMOVED******REMOVED***= the stepper step size for the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- [0,nGrid-1] 
***REMOVED******REMOVED***  double t0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // t0***REMOVED******REMOVED******REMOVED***  = the first boundary value of the independent variable
***REMOVED******REMOVED***  double tm;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // tm***REMOVED******REMOVED******REMOVED***  = the last boundary value of the independent variable
***REMOVED******REMOVED***  double _k_epsilon_J;***REMOVED******REMOVED******REMOVED***// _k_epsilon_J = the perturbation parameter for a state variable at a particular iteration k
***REMOVED******REMOVED***  double _k_alpha;***REMOVED******REMOVED******REMOVED******REMOVED*** // _k_alpha***REMOVED***  = the relaxation factor at a particular iteration k to scale the adjustment to the initial condition 
***REMOVED******REMOVED***  double _k_G;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // _k_G***REMOVED******REMOVED******REMOVED***= the Root Mean Square (RMS) error of boundary residues at a particular iteration k
***REMOVED******REMOVED***  double _k_GPrev;***REMOVED******REMOVED******REMOVED******REMOVED*** // _k_GPrev***REMOVED***  = the Root Mean Square (RMS) error of boundary residues at the previous iteration k-1
***REMOVED******REMOVED***  RowVectorXd _k_tSol(nGrid);  // _k_tSol***REMOVED******REMOVED***= the independent variable t over the whole grid in the solution of the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (1xnGrid)
***REMOVED******REMOVED***  MatrixXd _k_xSol(n,nGrid);***REMOVED***// _k_xSol***REMOVED******REMOVED***= the state vector x integrated over the whole grid in the solution of the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nxnGrid)
***REMOVED******REMOVED***  RowVectorXd _t_solPJ(nGrid); // _t_solPJ***REMOVED***  = the independent variable t over the whole grid in the perturbed solution of the IVP solver***REMOVED******REMOVED******REMOVED***-- (1xnGrid)***REMOVED*** 
***REMOVED******REMOVED***  MatrixXd _x_solPJ(n,nGrid);  // _x_solPJ***REMOVED***  = the state vector x integrated over the whole grid in the perturbed solution of the IVP solver***REMOVED******REMOVED***-- (nxnGrid)
***REMOVED******REMOVED***  VectorXd _k_x_t1(n);***REMOVED******REMOVED******REMOVED***// _k_x_t1***REMOVED******REMOVED***= the computed initial state vector in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED******REMOVED***  VectorXd _k_x_t1P(n);***REMOVED******REMOVED***  // _k_x_t1***REMOVED******REMOVED***= the computed perturbed initial state vector in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED***  VectorXd x_t1(n);***REMOVED******REMOVED******REMOVED******REMOVED***// x_t1***REMOVED******REMOVED******REMOVED***= the computed initial state vector to be input to the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED******REMOVED***  VectorXd x_t1P(n);***REMOVED******REMOVED******REMOVED***  // x_t1P***REMOVED******REMOVED***  = the computed perturbed initial state vector to be input to the IVP solver***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  -- (nx1)
***REMOVED******REMOVED***  VectorXd _k_g(n);***REMOVED******REMOVED******REMOVED******REMOVED***// _k_g***REMOVED******REMOVED******REMOVED***= the boundary condition residues in the k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED***  MatrixXd _k_S(n,n);***REMOVED******REMOVED******REMOVED*** // _k_S***REMOVED******REMOVED******REMOVED***= the adjusting matrix for correcting the initial condition k-th iteration***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nxn) 

***REMOVED******REMOVED***  // Variable definitions
***REMOVED******REMOVED***  h = (tm - t0)/nGrid;
***REMOVED******REMOVED***  t0 = t_BC(0);
***REMOVED******REMOVED***  tm = t_BC(m-1);

***REMOVED******REMOVED***  // Capture function calls by the ODEInt library for differentials and convert it to a custom form 
***REMOVED******REMOVED***  auto dFunctionWrapper = [dxBydt] (const VectorXd &x, VectorXd &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dxBydt(t, x);
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  // Observer to handle the solutions of the IVP Solver
***REMOVED******REMOVED***  auto odeIObserver = [k, n, &IVPIColumnIndex, &IVPIXSolutions, &IVPITSolutions] (const VectorXd &x , const double t){
***REMOVED******REMOVED******REMOVED******REMOVED***IVPIXSolutions.col(IVPIColumnIndex) = x;
***REMOVED******REMOVED******REMOVED******REMOVED***IVPITSolutions(IVPIColumnIndex) = t;
***REMOVED******REMOVED******REMOVED******REMOVED***++IVPIColumnIndex;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  };***REMOVED*** 
***REMOVED******REMOVED***  auto odePObserver = [k, n, &IVPPColumnIndex, &IVPPXSolutions, &IVPPTSolutions] (const VectorXd &x , const double t){
***REMOVED******REMOVED******REMOVED******REMOVED***IVPPXSolutions.col(IVPPColumnIndex) = x;
***REMOVED******REMOVED******REMOVED******REMOVED***IVPPTSolutions(IVPPColumnIndex) = t;
***REMOVED******REMOVED******REMOVED******REMOVED***++IVPPColumnIndex;
***REMOVED******REMOVED***  };***REMOVED******REMOVED***

***REMOVED******REMOVED***  auto getBCs = [n, m, t_BC] (RowVectorXd IVPTSolutions, MatrixXd xSolutions) -> MatrixXd{
***REMOVED******REMOVED******REMOVED******REMOVED***return xSolutions(Eigen::all, ((t_BC-t_BC(0)*RowVectorXd::Ones(m))/h).array().round().cast<int>());
***REMOVED******REMOVED***  };***REMOVED******REMOVED***  

***REMOVED******REMOVED***  // Solve the initial value problem for the first time  
***REMOVED******REMOVED***  kX0Temp = kX0;
***REMOVED******REMOVED***  cout<<"Starting I solver..."<<endl;
***REMOVED******REMOVED***  integrate_const(StepperType(), dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver); 
***REMOVED******REMOVED***  if(IVPIColumnIndex < nGrid){
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"OMG! I Solver is still running..."<<" k = "<<-1<<"... IVPIColumnIndex = "<<IVPIColumnIndex<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***return 0;
***REMOVED******REMOVED***  }***REMOVED******REMOVED***  
***REMOVED******REMOVED***  IVPIColumnIndex = 0;  
***REMOVED******REMOVED***  gkX0 = BCResidue(getBCs(IVPITSolutions, IVPIXSolutions));
***REMOVED******REMOVED***  // cout<<"gkX0 = "<<endl<<gkX0<<endl;
***REMOVED******REMOVED***  kP1G = gkX0.norm();
***REMOVED******REMOVED***  kG = kP1G;
***REMOVED******REMOVED***  while(kP1G > SIGMA){  
***REMOVED******REMOVED******REMOVED******REMOVED***if(kP1G < 0.1*kG) {
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"Going too fast. Changing ALPHA from "<<kAlpha<<" to "<<fmin(1.2*kAlpha, 1.0)<<"..."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kAlpha = fmin(1.2*kAlpha, 1.0);***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED***} else if(kP1G >= kG){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"Oops. Error increased. To make it go faster, changing ALPHA from "<<kAlpha<<" to "<<0.8*kAlpha<<"..."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kAlpha = 0.8*kAlpha;
***REMOVED******REMOVED******REMOVED******REMOVED***}***REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***for(int j = 0; j < n; j++){***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Determine the perturbation parameter
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** //kEpsilon = max(EPSION, abs(EPSION * kX0(j)));
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kEpsilon = EPSION;
***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Perturb the initial conditions***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(j);
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Solve the perturbed initial value problem***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** integrate_const(StepperType(), dFunctionWrapper, pX, t0, tm, h, odePObserver);
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** if(IVPPColumnIndex < nGrid){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  cout<<"OMG! P Solver is still running..."<<" k = "<<k<<"... IVPPColumnIndex = "<<IVPPColumnIndex<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  return 0;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** }***REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** IVPPColumnIndex = 0;

***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Compute a column of the adjusting matrix***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** S.col(j) = (BCResidue(getBCs(IVPPTSolutions, IVPPXSolutions))- gkX0)/kEpsilon;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // cout<<"gkxp = "<<endl<<BCResidue(getBCs(IVPPTSolutions, IVPPXSolutions))<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***}
***REMOVED******REMOVED******REMOVED******REMOVED***// cout<<"S = "<<endl<<S<<endl;

***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the linarized adjusting equation
***REMOVED******REMOVED******REMOVED******REMOVED***kX0 = S.colPivHouseholderQr().solve(-kAlpha*gkX0) + kX0;
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"Change in x0 = "<<endl<<S.completeOrthogonalDecomposition().solve(-kAlpha*gkX0)<<endl;***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***//kX0 = kX0 - kAlpha*S.inverse()*gkX0;

***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the initial value problem***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***kX0Temp = kX0;
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(StepperType(), dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver);  
***REMOVED******REMOVED******REMOVED******REMOVED***if(IVPIColumnIndex < nGrid){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"OMG! I Solver is still running..."<<" k = "<<-1<<"... IVPIColumnIndex = "<<IVPIColumnIndex<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** return 0;
***REMOVED******REMOVED******REMOVED******REMOVED***}***REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED***IVPIColumnIndex = 0; 
***REMOVED******REMOVED******REMOVED******REMOVED***gkX0 = BCResidue(getBCs(IVPITSolutions, IVPIXSolutions));

***REMOVED******REMOVED******REMOVED******REMOVED***kG = kP1G;
***REMOVED******REMOVED******REMOVED******REMOVED***kP1G = gkX0.norm()/sqrt(n);
***REMOVED******REMOVED******REMOVED******REMOVED***if(kP1G > kG){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"-> kP1G > kG... kP1G = "<<kP1G<<"... kG = "<<kG<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** //cout<<"S = "<<endl<<S<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** //cout<<"S^-1 = "<<endl<<S.inverse()<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***}
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"kP1G = "<<kP1G<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***++k;
***REMOVED******REMOVED***  }  
***REMOVED******REMOVED***  cout<<"t(0) = "<<endl<<IVPITSolutions.col(0)<<endl;
***REMOVED******REMOVED***  cout<<"t(nGrid-1) = "<<endl<<IVPITSolutions.col(nGrid-1)<<endl;
***REMOVED******REMOVED***  cout<<"x(0) = "<<endl<<IVPIXSolutions.col(0)<<endl;
***REMOVED******REMOVED***  cout<<"x(nGrid-1) = "<<endl<<IVPIXSolutions.col(nGrid-1)<<endl;
***REMOVED******REMOVED***  //x0 = x0 + EPSION*MatrixXd::Identity(n,n).col(0);
***REMOVED******REMOVED***  // pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(0);
***REMOVED******REMOVED***  // integrate_const(IVPPStepper, dFunctionWrapper, x0, t0, tm, h, odePObserver);
***REMOVED******REMOVED***  // IVPPColumnIndex = 0; 

***REMOVED******REMOVED***  // Perturb the initial conditions***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  // x0 = x0 + kEpsilon*MatrixXd::Identity(n,n).col(1);***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED***  // Solve the perturbed initial value problem***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  // integrate_const(IVPPStepper, dFunctionWrapper, x0, t0, tm, h, odePObserver);
***REMOVED******REMOVED***  // IVPPColumnIndex = 0;
***REMOVED******REMOVED***  // cout<<"Check BCs = "<<endl<<BCResidue(getBCs(IVPPTSolutions, IVPPXSolutions))<<endl;
***REMOVED******REMOVED***  return 0;
***REMOVED*** }
// ===================


