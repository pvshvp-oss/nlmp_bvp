#include <nlmp_bvp.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>

using namespace Eigen;
using namespace boost::numeric::odeint;
typedef Matrix<double, 1, Dynamic> ArrayXd;
typedef Matrix<int, 1, Dynamic> ArrayXi;

int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** ArrayXd samples, ArrayXi nodeID,
***REMOVED*** VectorXd (*dFunction) (double t, int intervalID),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** ){  
***REMOVED******REMOVED***  int nSamples = samples.size();
***REMOVED******REMOVED***  int nIntervals = nodeID.size() - 1;
***REMOVED******REMOVED***  int k = 0, j = 1;  
***REMOVED******REMOVED***  do{

***REMOVED******REMOVED***  }while(true);
***REMOVED******REMOVED***  return 0;
***REMOVED*** }



