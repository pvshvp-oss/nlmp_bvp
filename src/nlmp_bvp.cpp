#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>

using namespace Eigen;
using namespace boost::numeric::odeint;
typedef VectorXd state_type;

const double STEPPER_STEP = 1e-5;

int nlmp_bvp(
    int nEquations,
    VectorXd startingState,
    VectorXd tNodes,
    VectorXd (*dFunction) (double t, VectorXd x),
    VectorXd (*BCFunction) (MatrixXd BC)
    ){  
        auto dFunctionWrapper = [] (const VectorXd &x, VectorXd &dxdt, double t){
            dxdt = dFunction(t, x);
        }
        int k = 0, j = 1;  
        do{
            integrate(dFunctionWrapper, startingState, tNodes(1), tNodes(tNodes.rows() - 1), STEPPER_STEP);
        }while(true);
        return 0;
    }



