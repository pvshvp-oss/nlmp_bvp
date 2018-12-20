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

const double h =0.1;
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
            RowVectorXi BCIndices;
            BCIndices.resize(1, m);
            BCIndices = ((tNodes-tNodes(0)*RowVectorXd::Ones(m))/h).array().round().cast<int>();
            // cout<<"<Inside getBCs>"<<endl;
            // cout<<"BCIndices = "<<BCIndices<<endl;
            // cout<<"xSolutions(Eigen::all, BCIndices) = "<<xSolutions(Eigen::all, BCIndices)<<endl;
            // cout<<"</Inside getBCs>"<<endl<<endl;
            return xSolutions(Eigen::all, BCIndices);
        };

        do{            
            // Solve the initial value problem   
            integrate_const(stepper, dFunctionWrapper, x0, t0, tm, h, odeIObserver);   
            IVPIColumnIndex = 0;  

            // cout<<"<IVPI solutions>"<<endl;
            // cout<<"x = "<<IVPIXSolutions<<endl;
            // cout<<"t = "<<IVPITSolutions<<endl;
            // cout<<"</IVPI solutions>"<<endl<<endl;

            for(int j = 0; j < n; j++){   
                // Determine the perturbation parameter
                kEpsilon = max(epsilon, abs(epsilon * kX0(j,0)));
            
                // Perturb the initial conditions            
                pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(j);

                // Solve the perturbed initial value problem            
                integrate_const(stepper, dFunctionWrapper, pX, t0, tm, h, odePObserver);
                IVPPColumnIndex = 0;

                // cout<<"<IVPP solutions>"<<endl;
                // cout<<"x = "<<IVPPXSolutions<<endl;
                // cout<<"t = "<<IVPPTSolutions<<endl;
                // cout<<"</IVPP solutions>"<<endl<<endl;

                // cout<<"<Calculated BCs>"<<endl;
                // cout<<"BC for IVPI = "<<getBCs(IVPITSolutions, IVPIXSolutions)<<endl;
                // cout<<"BC for IVPP = "<<getBCs(IVPPTSolutions, IVPPXSolutions)<<endl;                
                // cout<<"</Calculated BCs>"<<endl<<endl;

                // cout<<"<Calculated residues>"<<endl;
                // cout<<"Residues for IVPI = "<<BCFunction(getBCs(IVPITSolutions, IVPIXSolutions))<<endl;
                // cout<<"Residues for IVPP = "<<BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))<<endl;                
                // cout<<"</Calculated residues>"<<endl<<endl;

                // Compute a column of the adjsuting matrix
                S.col(j) = (BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))- BCFunction(getBCs(IVPITSolutions, IVPIXSolutions)))/kEpsilon;                
            }
            cout<<"S = "<<S<<endl;
            break;
        }while(true);    
        return 0;
    }
// ================


