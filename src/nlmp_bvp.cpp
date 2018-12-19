// ========================
// Includes and definitions
// ========================
#include <cmath>
#include <algorithm>
#include <Eigen/Eigen>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
// #include <boost/multiprecision/eigen.hpp>
#include <nlmp_bvp.hpp>
#include <mlinterp.hpp>
using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;
using namespace mlinterp;
using RowVectorXd = Matrix<double, 1, Dynamic>;
using RowVectorXi = Matrix<int, 1, Dynamic>;

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
        int m;
        int k = 0;     
        int IVPIColumnIndex = 0;
        int IVPPColumnIndex = 0;        
        int nSamples = 0;  
        double kEpsilon = 0;
        double t0;
        double tm;    
        m = tNodes.cols();  
        RowVectorXd IVPITSolutions, IVPPTSolutions;
        RowVectorXi boundaryColumns;
        VectorXd kX0;
        VectorXd pX;
        MatrixXd S;         
        MatrixXd IVPIXSolutions, IVPPXSolutions; 
        runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> stepper;

        t0 = tNodes(0);
        tm = tNodes(m-1);
        nSamples = (int)((tm - t0)/h) + 1;
        IVPIXSolutions.resize(n, nSamples);
        IVPITSolutions.resize(1, nSamples);
        IVPPXSolutions.resize(n, nSamples);
        IVPPTSolutions.resize(1, nSamples);
        boundaryColumns.resize(1, m);
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

        auto getBCs = [n, m, tNodes] (RowVectorXd IVPTSolutions, MatrixXd xSolutions) -> MatrixXd{
            RowVectorXi BCIndices, BCIndicesNext;
            MatrixXd BCs;
            BCs.resize(n,m);
            BCIndices.resize(1, m);
            BCIndices = ((tNodes-tNodes(0)*RowVectorXd::Ones(m))/h).cast<int>();
            BCIndicesNext = BCIndices+ RowVectorXi::Ones(m);
            BCIndices = ((tNodes.array() - IVPTSolutions(Eigen::all, BCIndices).array()) < (IVPTSolutions(Eigen::all, BCIndicesNext).array() - tNodes.array())).select(BCIndices, BCIndicesNext); 
            return xSolutions(Eigen::all, BCIndices);
        };

        do{            
            // Solve the initial value problem   
            integrate_const(stepper, dFunctionWrapper, x0, t0, tm, h, odeIObserver);   
            IVPIColumnIndex = 0;  

            // cout<<"Regular IVP solutions"<<endl;
            // cout<<"x = "<<IVPIXSolutions<<endl;
            // cout<<"t = "<<IVPITSolutions<<endl;
            // cout<<"*************************"<<endl;

            for(int j = 0; j < n; j++){   
                // Determine the perturbation parameter
                kEpsilon = max(epsilon, abs(epsilon * kX0(j,0)));
            
                // Perturb the initial conditions            
                pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(j);

                // Solve the perturbed initial value problem            
                integrate_const(stepper, dFunctionWrapper, pX, t0, tm, h, odePObserver);
                IVPPColumnIndex = 0;

                // cout<<"Perturbed IVP solutions"<<endl;
                // cout<<"x = "<<IVPPXSolutions<<endl;
                // cout<<"t = "<<IVPPTSolutions<<endl;
                // cout<<"*************************"<<endl;  
       
                // VectorXi boundaryColumns = (tNodes/h).cast<int>();
            
                S.col(j) = (BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))- BCFunction(getBCs(IVPITSolutions, IVPIXSolutions)));
                cout<<"S = "<<S<<endl;
                break;
            }
            break;
        }while(true);    
        return 0;
    }
// ================


