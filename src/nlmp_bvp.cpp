#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>

using namespace Eigen;
using namespace boost::numeric::odeint;
typedef VectorXd state_type;

const double STEPPER_STEP = 1e-5;

int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** VectorXd startingState,
***REMOVED*** VectorXd tNodes,
***REMOVED*** VectorXd (*dFunction) (double t, VectorXd x),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** ){  
***REMOVED******REMOVED***  auto dFunctionWrapper = [] (const VectorXd &x, VectorXd &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dFunction(t, x);
***REMOVED******REMOVED***  }
***REMOVED******REMOVED***  int k = 0, j = 1;  
***REMOVED******REMOVED***  do{
***REMOVED******REMOVED******REMOVED******REMOVED***integrate(dFunctionWrapper, startingState, tNodes(1), tNodes(tNodes.rows() - 1), STEPPER_STEP);
***REMOVED******REMOVED***  }while(true);
***REMOVED******REMOVED***  return 0;
***REMOVED*** }



