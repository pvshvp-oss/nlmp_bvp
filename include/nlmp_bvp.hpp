#include <Eigen/Eigen>

using namespace Eigen;
using RowVectorXd = Matrix<double, 1, Dynamic>;

int nlmp_bvp(
***REMOVED*** int nEquations,
***REMOVED*** VectorXd startingState,
***REMOVED*** RowVectorXd tNodes,
***REMOVED*** VectorXd (*dFunction) (double t, VectorXd x),
***REMOVED*** VectorXd (*BCFunction) (MatrixXd BC)
***REMOVED*** );