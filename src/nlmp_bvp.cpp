// ========================
// Includes and definitions
// ========================
#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <cmath>
using namespace Eigen;
using namespace boost::numeric::odeint;
using StateType = VectorXd;

const double STEPPER_STEP = 1e-5;
// ========================

// ================
// The BVP function
// ================
int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** VectorXd startingState,
***REMOVED*** VectorXd tNodes,
***REMOVED*** VectorXd dFunction(double t, VectorXd x),
***REMOVED*** VectorXd BCFunction(MatrixXd BC)
***REMOVED*** ){  
***REMOVED******REMOVED***  
***REMOVED******REMOVED***  auto dFunctionWrapper = [dFunction] (const StateType &x, StateType &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dFunction(t, x);
***REMOVED******REMOVED***  };
***REMOVED******REMOVED***  
***REMOVED******REMOVED***  int k = 0, j = 1;  

***REMOVED******REMOVED***  StateType startingStateDEMO;
***REMOVED******REMOVED***  startingStateDEMO << 1.0, 2.0;

***REMOVED******REMOVED***  runge_kutta_dopri5<StateType,double,StateType,double,vector_space_algebra> stepper;

***REMOVED******REMOVED***  do{
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(stepper, dFunctionWrapper, /*startingState*/ startingStateDEMO, tNodes(1), tNodes(tNodes.rows() - 1), STEPPER_STEP);
***REMOVED******REMOVED***  }while(true);
***REMOVED******REMOVED***  return 0;

***REMOVED*** }
// ================


