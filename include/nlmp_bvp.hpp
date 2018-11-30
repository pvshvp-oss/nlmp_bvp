#include <Eigen/Dense>

using namespace Eigen;

int nlmp_bvp(
    int nEquations, int nIntervals,
    VectorXd (*dFunction) (double t, int intervalID),
    VectorXd (*BCFunction) (MatrixXd BC)
    );