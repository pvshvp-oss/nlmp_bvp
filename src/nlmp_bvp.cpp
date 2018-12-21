// ===============================
// Includes and global definitions
// ===============================
#include <cmath>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // C++ analog to math.h
#include <Eigen/Eigen>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // For matrix and vector math
#include <boost/numeric/odeint.hpp>***REMOVED******REMOVED******REMOVED******REMOVED*** // For the initial value problem solver
#include <nlmp_bvp.hpp>***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // For the function declarations***REMOVED***
using namespace std;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** //
using namespace Eigen;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  //
using namespace boost::numeric::odeint;***REMOVED******REMOVED******REMOVED***//
using RowVectorXi = Matrix<int, 1, Dynamic>;***REMOVED*** // For the convenience of declaring row vectors
using RowVectorXd = Matrix<double, 1, Dynamic>; // For the convenience of declaring row vectors
const double EPSION = 1e-8;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// EPSION = the state perturbation parameter to probe the differential equation system with
const double ALPHA = 1;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // ALPHA  = the relaxation factor to scale the adjustment to the initial condition
const double SIGMA = 1e-14;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// SIGMA  = the tolerance for error outside which the solver needs to  iterate further. 
const double BETA = 1e-3;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // BETA***REMOVED***= the deflation factor
// ===============================

// ===================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
***REMOVED*** int n,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// n***REMOVED******REMOVED******REMOVED***= the number of differential equations = the number of boundary conditions
***REMOVED*** int m,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// m***REMOVED******REMOVED******REMOVED***= the number of nodes at which boundary conditions are specified
***REMOVED*** int nGrid,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // nGrid***REMOVED***  = the number of points at which the state is evaluated
***REMOVED*** RowVectorXd t_BC,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // t_BC***REMOVED******REMOVED***= row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXd _0x_t1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***  // _0x_t1***REMOVED*** = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)***REMOVED*** 
***REMOVED*** VectorXd dxBydt(double t, VectorXd x), // dxBydt***REMOVED*** = a function that defines the derivative of a state vector x at t***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)
***REMOVED*** VectorXd BCResidue(MatrixXd x_BC)***REMOVED******REMOVED***// BCResidue = a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
***REMOVED*** ){  
***REMOVED******REMOVED***  int m;
***REMOVED******REMOVED***  int k = 0;***REMOVED***  
***REMOVED******REMOVED***  int IVPIColumnIndex = 0;
***REMOVED******REMOVED***  int IVPPColumnIndex = 0;***REMOVED******REMOVED***  
***REMOVED******REMOVED***  int nGrid = 0;  
***REMOVED******REMOVED***  double kEpsilon = 0;
***REMOVED******REMOVED***  double kAlpha = ALPHA;
***REMOVED******REMOVED***  double t0;
***REMOVED******REMOVED***  double tm;***REMOVED*** 
***REMOVED******REMOVED***  double kG, kP1G;
***REMOVED******REMOVED***  m = t_BC.cols();  
***REMOVED******REMOVED***  RowVectorXd IVPITSolutions, IVPPTSolutions;
***REMOVED******REMOVED***  VectorXd kX0;
***REMOVED******REMOVED***  VectorXd kX0Temp;
***REMOVED******REMOVED***  VectorXd pX;
***REMOVED******REMOVED***  VectorXd gkX0;
***REMOVED******REMOVED***  MatrixXd S;***REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  MatrixXd IVPIXSolutions, IVPPXSolutions; 
***REMOVED******REMOVED***  runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> IVPIStepper;
***REMOVED******REMOVED***  runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> IVPPStepper;

***REMOVED******REMOVED***  t0 = t_BC(0);
***REMOVED******REMOVED***  tm = t_BC(m-1);
***REMOVED******REMOVED***  nGrid = floor((tm - t0)/h) + 1;
***REMOVED******REMOVED***  IVPIXSolutions.resize(n, nGrid);
***REMOVED******REMOVED***  IVPITSolutions.resize(1, nGrid);
***REMOVED******REMOVED***  IVPPXSolutions.resize(n, nGrid);
***REMOVED******REMOVED***  IVPPTSolutions.resize(1, nGrid);
***REMOVED******REMOVED***  kX0.resize(n,1);
***REMOVED******REMOVED***  kX0Temp.resize(n,1);
***REMOVED******REMOVED***  gkX0.resize(n,1);
***REMOVED******REMOVED***  pX.resize(n,1);
***REMOVED******REMOVED***  S.resize(n,n);***REMOVED******REMOVED*** 
***REMOVED******REMOVED***  kX0 = x0;

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
***REMOVED******REMOVED***  integrate_const(IVPIStepper, dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver); 
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
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** integrate_const(IVPPStepper, dFunctionWrapper, pX, t0, tm, h, odePObserver);
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
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(IVPIStepper, dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver);  
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


