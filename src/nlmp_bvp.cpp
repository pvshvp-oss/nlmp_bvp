#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <cmath>

using namespace Eigen;
using namespace boost::numeric::odeint;
typedef std::vector<double> state_type;

const double STEPPER_STEP = 1e-5;

int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** VectorXd startingState,
***REMOVED*** VectorXd tNodes,
***REMOVED*** VectorXd (*dFunction) (double t, VectorXd x),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** ){  
***REMOVED******REMOVED***  auto dFunctionWrapper = [/*&dFunction*/] (const state_type &x, state_type &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***// dxdt = dFunction(t, x);
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt.push_back(1.0);
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt.push_back(2.0);
***REMOVED******REMOVED***  };
***REMOVED******REMOVED***  int k = 0, j = 1;  

***REMOVED******REMOVED***  state_type startingStateDEMO;
***REMOVED******REMOVED***  startingStateDEMO.push_back(1.0);
***REMOVED******REMOVED***  startingStateDEMO.push_back(2.0);

***REMOVED******REMOVED***  do{
***REMOVED******REMOVED******REMOVED******REMOVED***integrate(dFunctionWrapper, /*startingState*/ startingStateDEMO, tNodes(1), tNodes(tNodes.rows() - 1), STEPPER_STEP);
***REMOVED******REMOVED***  }while(true);
***REMOVED******REMOVED***  return 0;
***REMOVED*** }



