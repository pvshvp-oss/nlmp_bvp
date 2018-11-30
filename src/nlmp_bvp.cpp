#include <nlmp_bvp.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace Eigen;
using namespace boost::numeric::odeint;

int nlmp_bvp(
    int nEquations, int nIntervals,
    VectorXd (*dFunction) (double t, int intervalID),
    VectorXd (*BCFunction) (MatrixXd BC)
    ){      
        int k = 0, j = 1;  
        do{

        }while(true);
        return 0;
    }



