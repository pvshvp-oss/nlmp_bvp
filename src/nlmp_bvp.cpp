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
***REMOVED*** int nEquations,
***REMOVED*** StateType startingState,
***REMOVED*** StateType tNodes,
***REMOVED*** StateType dFunction(double t, StateType x),
***REMOVED*** StateType BCFunction(MatrixXd BC)
***REMOVED*** ){  
***REMOVED******REMOVED***  
***REMOVED******REMOVED***  // Capture function calls by the ODEInt library for differentials and convert it to a custom form 
***REMOVED******REMOVED***  auto dFunctionWrapper = [dFunction] (const StateType &x, StateType &dxdt, double t){
***REMOVED******REMOVED******REMOVED******REMOVED***dxdt = dFunction(t, x);
***REMOVED******REMOVED***  };

***REMOVED******REMOVED***  MatrixXd odeXSolutions;***REMOVED*** 
***REMOVED******REMOVED***  RowVectorXd odeTSolutions;
***REMOVED******REMOVED***  auto odeObserver = [nEquations, &odeXSolutions, &odeTSolutions] (const StateType &x , const double t){
***REMOVED******REMOVED******REMOVED******REMOVED***odeXSolutions.conservativeResize(nEquations, odeXSolutions.cols()+1);
***REMOVED******REMOVED******REMOVED******REMOVED***odeXSolutions.col(odeXSolutions.cols()-1) = x;
***REMOVED******REMOVED******REMOVED******REMOVED***odeTSolutions.conservativeResize(1, odeTSolutions.cols()+1);
***REMOVED******REMOVED******REMOVED******REMOVED***odeTSolutions(0, odeTSolutions.cols()-1) = t;
***REMOVED******REMOVED***  };***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***

***REMOVED******REMOVED***  int k = 0, j = 1;  

***REMOVED******REMOVED***  runge_kutta_dopri5<StateType,double,StateType,double,vector_space_algebra> stepper;
***REMOVED******REMOVED***  do{
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"Before integration..."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***integrate_const(stepper, dFunctionWrapper, /*startingState*/ startingState, tNodes(0), tNodes(tNodes.rows() - 1), STEPPER_STEP, odeObserver);
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"x = "<<odeXSolutions<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"t = "<<odeTSolutions<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"*************************"<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***cout<<"After integration..."<<endl;
***REMOVED******REMOVED******REMOVED******REMOVED***break;
***REMOVED******REMOVED***  }while(true);

***REMOVED******REMOVED***  return 0;

***REMOVED*** }
// ================


