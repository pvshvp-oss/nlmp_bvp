// ========================
// Includes and definitions
// ========================
#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <cmath>
using namespace Eigen;
using namespace boost::numeric::odeint;
using StateType = VectorXd;

const double STEPPER_STEP = 1e-5;
// ========================

// ================
// The BVP function
// ================
int nlmp_bvp(
    int nEquations,
    StateType startingState,
    StateType tNodes,
    StateType dFunction(double t, StateType x),
    StateType BCFunction(MatrixXd BC)
    ){  
        
        // Capture function calls by the ODEInt library for differentials and convert it to a custom form 
        auto dFunctionWrapper = [dFunction] (const StateType &x, StateType &dxdt, double t){
            dxdt = dFunction(t, x);
        };
        
        int k = 0, j = 1;  

        runge_kutta_dopri5<StateType,double,StateType,double,vector_space_algebra> stepper;

        do{
            integrate_const(stepper, dFunctionWrapper, /*startingState*/ startingState, tNodes(1), tNodes(tNodes.rows() - 1), STEPPER_STEP);
            
        }while(true);
        return 0;

    }
// ================


