#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

int nlmp_bvp(
***REMOVED*** int nEquations, int nIntervals,
***REMOVED*** VectorXd (*dFunction) (double t, int intervalID),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** );