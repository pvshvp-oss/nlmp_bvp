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
    int nEquations,
    ArrayXd samples, ArrayXi nodeID,
    VectorXd (*dFunction) (double t, int intervalID),
    VectorXd (*BCFunction) (MatrixXd BC)
    ){  
        int nSamples = samples.size();
        int nIntervals = nodeID.size() - 1;
        int k = 0, j = 1;  
        do{

        }while(true);
        return 0;
    }



