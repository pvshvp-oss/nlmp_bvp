#include <Eigen/Eigen>

using namespace Eigen;
using RowVectorXd = Matrix<double, 1, Dynamic>;

int nlmp_bvp(
    int nEquations,
    VectorXd startingState,
    RowVectorXd tNodes,
    VectorXd (*dFunction) (double t, VectorXd x),
    VectorXd (*BCFunction) (MatrixXd BC)
    );