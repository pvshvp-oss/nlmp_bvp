#include <Eigen/Eigen>

using namespace Eigen;

int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** VectorXd startingState,
***REMOVED*** VectorXd tNodes,
***REMOVED*** VectorXd (*dFunction) (double t, VectorXd x),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** );