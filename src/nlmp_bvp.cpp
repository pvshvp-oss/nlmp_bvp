// ========================
// Includes and definitions
// ========================
#include <cmath>
#include <algorithm>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
// #include <boost/multiprecision/eigen.hpp>
#include <nlmp_bvp.hpp>
#include <mlinterp.hpp>
using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;
using namespace mlinterp;
using RowVectorXd = Matrix<double, 1, Dynamic>;
using RowVectorXi = Matrix<int, 1, Dynamic>;

const double h = 4e-3;
const double epsilon = 1e-10;
const double alpha = 1;
const double sigma = 1e-14;
const double beta = 1e-3;
// ========================

// ================
// The BVP function
// ================
int nlmp_bvp(
***REMOVED*** int n,
***REMOVED*** VectorXd x0,
***REMOVED*** RowVectorXd tNodes,
***REMOVED*** VectorXd dFunction(double t, VectorXd x),
***REMOVED*** VectorXd BCFunction(MatrixXd BC)
***REMOVED*** ){  
***REMOVED******REMOVED***  int m;
***REMOVED******REMOVED***  int k = 0;***REMOVED***  
***REMOVED******REMOVED***  int IVPIColumnIndex = 0;
***REMOVED******REMOVED***  int IVPPColumnIndex = 0;***REMOVED******REMOVED***  
***REMOVED******REMOVED***  int nSamples = 0;  
***REMOVED******REMOVED***  double kEpsilon = 0;
***REMOVED******REMOVED***  double kAlpha = alpha;
***REMOVED******REMOVED***  double t0;
***REMOVED******REMOVED***  double tm;***REMOVED*** 
***REMOVED******REMOVED***  double kG, kP1G;
***REMOVED******REMOVED***  m = tNodes.cols();  
***REMOVED******REMOVED***  RowVectorXd IVPITSolutions, IVPPTSolutions;
***REMOVED******REMOVED***  RowVectorXi boundaryColumns;
***REMOVED******REMOVED***  VectorXd kX0;
***REMOVED******REMOVED***  VectorXd kX0Temp;
***REMOVED******REMOVED***  VectorXd pX;
***REMOVED******REMOVED***  VectorXd gkX0;
***REMOVED******REMOVED***  MatrixXd S;***REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  MatrixXd IVPIXSolutions, IVPPXSolutions; 
***REMOVED******REMOVED***  runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> IVPIStepper;
***REMOVED******REMOVED***  runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> IVPPStepper;

***REMOVED******REMOVED***  t0 = tNodes(0);
***REMOVED******REMOVED***  tm = tNodes(m-1);
***REMOVED******REMOVED***  nSamples = floor((tm - t0)/h) + 1;
***REMOVED******REMOVED***  IVPIXSolutions.resize(n, nSamples);
***REMOVED******REMOVED***  IVPITSolutions.resize(1, nSamples);
***REMOVED******REMOVED***  IVPPXSolutions.resize(n, nSamples);
***REMOVED******REMOVED***  IVPPTSolutions.resize(1, nSamples);
***REMOVED******REMOVED***  boundaryColumns.resize(1, m);
***REMOVED******REMOVED***  kX0.resize(n,1);
***REMOVED******REMOVED***  kX0Temp.resize(n,1);
***REMOVED******REMOVED***  gkX0.resize(n,1);
***REMOVED******REMOVED***  pX.resize(n,1);
***REMOVED******REMOVED***  S.resize(n,n);***REMOVED******REMOVED*** 
***REMOVED******REMOVED***  kX0 = x0;

***REMOVED******REMOVED***  // Capture function calls by the ODEInt library for differentials and convert it to a custom form 
***REMOVED******REMOVED***  auto dFunctionWrapper = [dFunction] (const VectorXd &x, VectorXd &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dFunction(t, x);
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  // Observer to handle the solutions of the IVP Solver
***REMOVED******REMOVED***  auto odeIObserver = [n, &IVPIColumnIndex, &IVPIXSolutions, &IVPITSolutions] (const VectorXd &x , const double t){
***REMOVED******REMOVED******REMOVED******REMOVED***IVPIXSolutions.col(IVPIColumnIndex) = x;
***REMOVED******REMOVED******REMOVED******REMOVED***IVPITSolutions(IVPIColumnIndex) = t;
***REMOVED******REMOVED******REMOVED******REMOVED***++IVPIColumnIndex;
***REMOVED******REMOVED***  };***REMOVED*** 
***REMOVED******REMOVED***  auto odePObserver = [n, &IVPPColumnIndex, &IVPPXSolutions, &IVPPTSolutions] (const VectorXd &x , const double t){
***REMOVED******REMOVED******REMOVED******REMOVED***IVPPXSolutions.col(IVPPColumnIndex) = x;
***REMOVED******REMOVED******REMOVED******REMOVED***IVPPTSolutions(IVPPColumnIndex) = t;
***REMOVED******REMOVED******REMOVED******REMOVED***++IVPPColumnIndex;
***REMOVED******REMOVED***  };***REMOVED******REMOVED***

***REMOVED******REMOVED***  auto getBCs = [n, m, tNodes] (RowVectorXd IVPTSolutions, MatrixXd xSolutions) -> MatrixXd{
***REMOVED******REMOVED******REMOVED******REMOVED***RowVectorXi BCIndices;
***REMOVED******REMOVED******REMOVED******REMOVED***BCIndices.resize(1, m);
***REMOVED******REMOVED******REMOVED******REMOVED***BCIndices = ((tNodes-tNodes(0)*RowVectorXd::Ones(m))/h).array().round().cast<int>();
***REMOVED******REMOVED******REMOVED******REMOVED***return xSolutions(Eigen::all, BCIndices);
***REMOVED******REMOVED***  };***REMOVED******REMOVED***  

***REMOVED******REMOVED***  // Solve the initial value problem for the first time  
***REMOVED******REMOVED***  kX0Temp = kX0;
***REMOVED******REMOVED***  integrate_const(IVPIStepper, dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver); 
***REMOVED******REMOVED***  IVPIColumnIndex = 0;  
***REMOVED******REMOVED***  gkX0 = BCFunction(getBCs(IVPITSolutions, IVPIXSolutions));
***REMOVED******REMOVED***  // cout<<"gkX0 = "<<endl<<gkX0<<endl;
***REMOVED******REMOVED***  kP1G = gkX0.norm();
***REMOVED******REMOVED***  kG = kP1G;
***REMOVED******REMOVED***  while(kP1G > sigma){  
***REMOVED******REMOVED******REMOVED******REMOVED***if(kP1G < 0.1*kG) {
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"Going too fast. Changing alpha from "<<kAlpha<<" to "<<min(1.2*kAlpha, 1.0)<<"..."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kAlpha = min(1.2*kAlpha, 1.0);***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED***} else if(kP1G >= kG){
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"Oops. Error increased. To make it go faster, changing alpha from "<<kAlpha<<" to "<<0.8*kAlpha<<"..."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kAlpha = 0.8*kAlpha;
***REMOVED******REMOVED******REMOVED******REMOVED***}***REMOVED***  
***REMOVED******REMOVED******REMOVED******REMOVED***for(int j = 0; j < n; j++){***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Determine the perturbation parameter
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kEpsilon = max(epsilon, abs(epsilon * kX0(j)));
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** //kEpsilon = epsilon;
***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Perturb the initial conditions***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(j);
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Solve the perturbed initial value problem***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** integrate_const(IVPPStepper, dFunctionWrapper, pX, t0, tm, h, odePObserver);
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** IVPPColumnIndex = 0;

***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Compute a column of the adjusting matrix***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** S.col(j) = (BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))- gkX0)/kEpsilon;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // cout<<"gkxp = "<<endl<<BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***}
***REMOVED******REMOVED******REMOVED******REMOVED***// cout<<"S = "<<endl<<S<<endl;

***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the linarized adjusting equation
***REMOVED******REMOVED******REMOVED******REMOVED***kX0 = S.colPivHouseholderQr().solve(-kAlpha*gkX0) + kX0;

***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the initial value problem***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***kX0Temp = kX0;
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(IVPIStepper, dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver);***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***IVPIColumnIndex = 0; 
***REMOVED******REMOVED******REMOVED******REMOVED***gkX0 = BCFunction(getBCs(IVPITSolutions, IVPIXSolutions));

***REMOVED******REMOVED******REMOVED******REMOVED***kG = kP1G;
***REMOVED******REMOVED******REMOVED******REMOVED***kP1G = gkX0.norm()/sqrt(n);
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"kP1G = "<<kP1G<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***++k;
***REMOVED******REMOVED***  }  
***REMOVED******REMOVED***  cout<<"x(0) = "<<endl<<IVPIXSolutions.col(0)<<endl;
***REMOVED******REMOVED***  cout<<"x(nSamples-1) = "<<endl<<IVPIXSolutions.col(nSamples-1)<<endl;
***REMOVED******REMOVED***  //x0 = x0 + epsilon*MatrixXd::Identity(n,n).col(0);
***REMOVED******REMOVED***  // pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(0);
***REMOVED******REMOVED***  // integrate_const(IVPPStepper, dFunctionWrapper, x0, t0, tm, h, odePObserver);
***REMOVED******REMOVED***  // IVPPColumnIndex = 0; 

***REMOVED******REMOVED***  // Perturb the initial conditions***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  // x0 = x0 + kEpsilon*MatrixXd::Identity(n,n).col(1);***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** 
***REMOVED******REMOVED***  // Solve the perturbed initial value problem***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  // integrate_const(IVPPStepper, dFunctionWrapper, x0, t0, tm, h, odePObserver);
***REMOVED******REMOVED***  // IVPPColumnIndex = 0;
***REMOVED******REMOVED***  // cout<<"Check BCs = "<<endl<<BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))<<endl;
***REMOVED******REMOVED***  return 0;
***REMOVED*** }
// ================


