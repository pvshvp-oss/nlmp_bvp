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

const double h = 2e-3;
const double epsilon = 10^(-8);
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
***REMOVED******REMOVED***  double t0;
***REMOVED******REMOVED***  double tm;***REMOVED*** 
***REMOVED******REMOVED***  m = tNodes.cols();  
***REMOVED******REMOVED***  RowVectorXd IVPITSolutions, IVPPTSolutions;
***REMOVED******REMOVED***  RowVectorXi boundaryColumns;
***REMOVED******REMOVED***  VectorXd kX0;
***REMOVED******REMOVED***  VectorXd pX;
***REMOVED******REMOVED***  MatrixXd S;***REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED***  MatrixXd IVPIXSolutions, IVPPXSolutions; 
***REMOVED******REMOVED***  runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> stepper;

***REMOVED******REMOVED***  t0 = tNodes(0);
***REMOVED******REMOVED***  tm = tNodes(m-1);
***REMOVED******REMOVED***  nSamples = (int)((tm - t0)/h) + 1;
***REMOVED******REMOVED***  IVPIXSolutions.resize(n, nSamples);
***REMOVED******REMOVED***  IVPITSolutions.resize(1, nSamples);
***REMOVED******REMOVED***  IVPPXSolutions.resize(n, nSamples);
***REMOVED******REMOVED***  IVPPTSolutions.resize(1, nSamples);
***REMOVED******REMOVED***  boundaryColumns.resize(1, m);
***REMOVED******REMOVED***  kX0.resize(n,1);
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
***REMOVED******REMOVED******REMOVED******REMOVED***RowVectorXi BCIndices, BCIndicesNext;
***REMOVED******REMOVED******REMOVED******REMOVED***MatrixXd BCs;
***REMOVED******REMOVED******REMOVED******REMOVED***BCs.resize(n,m);
***REMOVED******REMOVED******REMOVED******REMOVED***BCIndices.resize(1, m);
***REMOVED******REMOVED******REMOVED******REMOVED***BCIndices = ((tNodes-tNodes(0)*RowVectorXd::Ones(m))/h).cast<int>();
***REMOVED******REMOVED******REMOVED******REMOVED***BCIndicesNext = BCIndices+ RowVectorXi::Ones(m);
***REMOVED******REMOVED******REMOVED******REMOVED***BCIndices = ((tNodes.array() - IVPTSolutions(Eigen::all, BCIndices).array()) < (IVPTSolutions(Eigen::all, BCIndicesNext).array() - tNodes.array())).select(BCIndices, BCIndicesNext); 
***REMOVED******REMOVED******REMOVED******REMOVED***return xSolutions(Eigen::all, BCIndices);
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  do{***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***// Solve the initial value problem***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(stepper, dFunctionWrapper, x0, t0, tm, h, odeIObserver);***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED***IVPIColumnIndex = 0;  

***REMOVED******REMOVED******REMOVED******REMOVED***// cout<<"Regular IVP solutions"<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***// cout<<"x = "<<IVPIXSolutions<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***// cout<<"t = "<<IVPITSolutions<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***// cout<<"*************************"<<endl;

***REMOVED******REMOVED******REMOVED******REMOVED***for(int j = 0; j < n; j++){***REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Determine the perturbation parameter
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** kEpsilon = max(epsilon, abs(epsilon * kX0(j,0)));
***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Perturb the initial conditions***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(j);

***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // Solve the perturbed initial value problem***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** integrate_const(stepper, dFunctionWrapper, pX, t0, tm, h, odePObserver);
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** IVPPColumnIndex = 0;

***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // cout<<"Perturbed IVP solutions"<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // cout<<"x = "<<IVPPXSolutions<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // cout<<"t = "<<IVPPTSolutions<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // cout<<"*************************"<<endl;  
***REMOVED******REMOVED*** 
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // VectorXi boundaryColumns = (tNodes/h).cast<int>();
***REMOVED******REMOVED******REMOVED******REMOVED***
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** S.col(j) = (BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))- BCFunction(getBCs(IVPITSolutions, IVPIXSolutions)));
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** cout<<"S = "<<S<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** break;
***REMOVED******REMOVED******REMOVED******REMOVED***}
***REMOVED******REMOVED******REMOVED******REMOVED***break;
***REMOVED******REMOVED***  }while(true);***REMOVED*** 
***REMOVED******REMOVED***  return 0;
***REMOVED*** }
// ================


