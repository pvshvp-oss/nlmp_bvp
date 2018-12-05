#include <Eigen/Eigen>

using namespace Eigen;

int nlmp_bvp(
    int nEquations,
    VectorXd (*dFunction) (double t, VectorXd x, int intervalID),
    VectorXd (*BCFunction) (MatrixXd BC)
    );