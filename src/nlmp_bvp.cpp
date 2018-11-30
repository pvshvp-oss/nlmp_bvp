#include <nlmp_bvp.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace Eigen;
using namespace boost::numeric::odeint;

int nlmp_bvp(
***REMOVED*** int nEquations, int nIntervals,
***REMOVED*** VectorXd (*dFunction) (double t, int intervalID),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** ){***REMOVED******REMOVED***
***REMOVED******REMOVED***  int k = 0, j = 1;  
***REMOVED******REMOVED***  do{

***REMOVED******REMOVED***  }while(true);
***REMOVED******REMOVED***  return 0;
***REMOVED*** }



