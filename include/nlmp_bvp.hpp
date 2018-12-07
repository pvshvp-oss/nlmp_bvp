#include <Eigen/Eigen>

using namespace Eigen;

int nlmp_bvp(
    int nEquations,
    VectorXd startingState,
    VectorXd tNodes,
    VectorXd (*dFunction) (double t, VectorXd x),
    VectorXd (*BCFunction) (MatrixXd BC)
    );