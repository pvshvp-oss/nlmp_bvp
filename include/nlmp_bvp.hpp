#include <Eigen/Eigen>

using namespace Eigen;

int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** VectorXd (*dFunction) (double t, VectorXd x, int intervalID),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** );