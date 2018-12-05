#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>

using namespace Eigen;
using namespace boost::numeric::odeint;
typedef VectorXd state_type;

int nlmp_bvp(
    int nEquations,
    VectorXd (*dFunction) (double t, VectorXd x, int intervalID),
    VectorXd (*BCFunction) (MatrixXd BC)
    ){  
        auto dFunctionWrapper = [] (const VectorXd &x, VectorXd &dxdt, double t){
            dxdt = dFunction(t, x, 0 /*TODO: Interval ID*/);
        }

        int k = 0, j = 1;  
        do{

        }while(true);
        return 0;
    }



