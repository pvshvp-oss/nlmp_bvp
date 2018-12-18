// ========================
// Includes and definitions
// ========================
#include <cmath>
#include <algorithm>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <nlmp_bvp.hpp>
using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;
using RowVectorXd = Matrix<double, 1, Dynamic>;

const double h = 2e-3;
const double epsilon = 10^(-8);
// ========================

// ================
// The BVP function
// ================
int nlmp_bvp(
    int n,
    VectorXd x0,
    RowVectorXd tNodes,
    VectorXd dFunction(double t, VectorXd x),
    VectorXd BCFunction(MatrixXd BC)
    ){  
        int k = 0;     
        int IVPIColumnIndex = 0;
        int IVPPColumnIndex = 0;
        int kEpsilon = 0;
        int nSamples = 0;  
        double t0;
        double tm;      
        RowVectorXd IVPITSolutions, IVPPTSolutions;
        RowVectorXi boundaryColumns;
        VectorXd kX0;
        VectorXd pX;
        VectorXd S;         
        MatrixXd IVPIXSolutions, IVPPXSolutions; 
        runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> stepper;

        t0 = tNodes(0);
        tm = tNodes(tNodes.cols()-1);
        nSamples = (int)((tm - t0)/h) + 1;
        IVPIXSolutions.resize(n, nSamples);
        IVPITSolutions.resize(1, nSamples);
        IVPPXSolutions.resize(n, nSamples);
        IVPPTSolutions.resize(1, nSamples);
        boundaryColumns.resize(1, tNodes.cols());
        kX0.resize(n,1);
        pX.resize(n,1);
        S.resize(n,n);        
        kX0 = x0;

        // Capture function calls by the ODEInt library for differentials and convert it to a custom form 
        auto dFunctionWrapper = [dFunction] (const VectorXd &x, VectorXd &dxdt, double t){
            dxdt = dFunction(t, x);
        };

        // Observer to handle the solutions of the IVP Solver
        auto odeIObserver = [n, &IVPIColumnIndex, &IVPIXSolutions, &IVPITSolutions] (const VectorXd &x , const double t){
            IVPIXSolutions.col(IVPIColumnIndex) = x;
            IVPITSolutions(IVPIColumnIndex) = t;
            ++IVPIColumnIndex;
        };    
        auto odePObserver = [n, &IVPPColumnIndex, &IVPPXSolutions, &IVPPTSolutions] (const VectorXd &x , const double t){
            IVPPXSolutions.col(IVPPColumnIndex) = x;
            IVPPTSolutions(IVPPColumnIndex) = t;
            ++IVPPColumnIndex;
        };      

        do{            
            // Solve the initial value problem   
            integrate_const(stepper, dFunctionWrapper, x0, t0, tm, h, odeIObserver);   
            IVPIColumnIndex = 0;  
            for(int j = 0; j < n; j++){   
                // Determine the perturbation parameter
                kEpsilon = max(epsilon, abs(epsilon * kX0(j,0)));
            
                // Perturb the initial conditions            
                pX = kX0 + kEpsilon;

                // Solve the perturbed initial value problem            
                integrate_const(stepper, dFunctionWrapper, pX, t0, tm, h, odePObserver);
                IVPPColumnIndex = 0;

                cout<<"x = "<<IVPXSolutions<<endl;
                cout<<"t = "<<IVPTSolutions<<endl;
                cout<<"*************************"<<endl;
                cout<<"After integration..."<<endl;     
                break;           
                // VectorXi boundaryColumns = (tNodes/h).cast<int>();
            
                S.col(j) = (BCFunction(IVPPXSolutions) - BCFunction(IVPIXSolutions));
            }
            break;
        }while(true);    
        return 0;
    }
// ================


