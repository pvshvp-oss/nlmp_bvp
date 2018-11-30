#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

int nlmp_bvp(
    int nEquations, int nIntervals,
    VectorXd (*dFunction) (double t, int intervalID),
    VectorXd (*BCFunction) (MatrixXd BC)
    );