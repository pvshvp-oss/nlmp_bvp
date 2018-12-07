#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <cmath>

using namespace Eigen;
using namespace boost::numeric::odeint;
typedef std::vector<double> state_type;

const double STEPPER_STEP = 1e-5;

int nlmp_bvp(
    int nEquations,
    VectorXd startingState,
    VectorXd tNodes,
    VectorXd (*dFunction) (double t, VectorXd x),
    VectorXd (*BCFunction) (MatrixXd BC)
    ){  
        auto dFunctionWrapper = [/*&dFunction*/] (const state_type &x, state_type &dxdt, double t){
            // dxdt = dFunction(t, x);
            dxdt.push_back(1.0);
            dxdt.push_back(2.0);
        };
        int k = 0, j = 1;  

        state_type startingStateDEMO;
        startingStateDEMO.push_back(1.0);
        startingStateDEMO.push_back(2.0);

        do{
            integrate(dFunctionWrapper, /*startingState*/ startingStateDEMO, tNodes(1), tNodes(tNodes.rows() - 1), STEPPER_STEP);
        }while(true);
        return 0;
    }



