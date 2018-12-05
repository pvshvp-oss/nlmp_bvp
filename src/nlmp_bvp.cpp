#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>

using namespace Eigen;
using namespace boost::numeric::odeint;
typedef VectorXd state_type;

int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** VectorXd (*dFunction) (double t, VectorXd x, int intervalID),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** ){  
***REMOVED******REMOVED***  auto dFunctionWrapper = [] (const VectorXd &x, VectorXd &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dFunction(t, x, 0 /*TODO: Interval ID*/);
***REMOVED******REMOVED***  }

***REMOVED******REMOVED***  int k = 0, j = 1;  
***REMOVED******REMOVED***  do{

***REMOVED******REMOVED***  }while(true);
***REMOVED******REMOVED***  return 0;
***REMOVED*** }



