// ========================
// Includes and definitions
// ========================
#include <nlmp_bvp.hpp>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <cmath>
using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;
using RowVectorXd = Matrix<double, 1, Dynamic>;
using StateType = VectorXd;

const double STEPPER_STEP = 1e-2;
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

        MatrixXd odeXSolutions;    
        RowVectorXd odeTSolutions;
        auto odeObserver = [nEquations, &odeXSolutions, &odeTSolutions] (const StateType &x , const double t){
            odeXSolutions.conservativeResize(nEquations, odeXSolutions.cols()+1);
            odeXSolutions.col(odeXSolutions.cols()-1) = x;
            odeTSolutions.conservativeResize(1, odeTSolutions.cols()+1);
            odeTSolutions(0, odeTSolutions.cols()-1) = t;
        };               

        int k = 0, j = 1;  

        runge_kutta_dopri5<StateType,double,StateType,double,vector_space_algebra> stepper;
        do{
            cout<<"Before integration..."<<endl;
            integrate_const(stepper, dFunctionWrapper, /*startingState*/ startingState, tNodes(0), tNodes(tNodes.rows() - 1), STEPPER_STEP, odeObserver);
            cout<<"x = "<<odeXSolutions<<endl;
            cout<<"t = "<<odeTSolutions<<endl;
            cout<<"*************************"<<endl;
            cout<<"After integration..."<<endl;
            break;
        }while(true);

        return 0;

    }
// ================


