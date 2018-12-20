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

const double h = 4e-3;
const double epsilon = 1e-10;
const double alpha = 1;
const double sigma = 1e-14;
const double beta = 1e-3;
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
        double kAlpha = alpha;
        double t0;
        double tm;    
        double kG, kP1G;
        m = tNodes.cols();  
        RowVectorXd IVPITSolutions, IVPPTSolutions;
        RowVectorXi boundaryColumns;
        VectorXd kX0;
        VectorXd kX0Temp;
        VectorXd pX;
        VectorXd gkX0;
        MatrixXd S;         
        MatrixXd IVPIXSolutions, IVPPXSolutions; 
        runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> IVPIStepper;
        runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> IVPPStepper;

        t0 = tNodes(0);
        tm = tNodes(m-1);
        nSamples = floor((tm - t0)/h) + 1;
        IVPIXSolutions.resize(n, nSamples);
        IVPITSolutions.resize(1, nSamples);
        IVPPXSolutions.resize(n, nSamples);
        IVPPTSolutions.resize(1, nSamples);
        boundaryColumns.resize(1, m);
        kX0.resize(n,1);
        kX0Temp.resize(n,1);
        gkX0.resize(n,1);
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
            return xSolutions(Eigen::all, BCIndices);
        };        

        // Solve the initial value problem for the first time  
        kX0Temp = kX0;
        integrate_const(IVPIStepper, dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver); 
        IVPIColumnIndex = 0;  
        gkX0 = BCFunction(getBCs(IVPITSolutions, IVPIXSolutions));
        // cout<<"gkX0 = "<<endl<<gkX0<<endl;
        kP1G = gkX0.norm();
        kG = kP1G;
        while(kP1G > sigma){  
            if(kP1G < 0.1*kG) {
                cout<<"Going too fast. Changing alpha from "<<kAlpha<<" to "<<min(1.2*kAlpha, 1.0)<<"..."<<endl;
                kAlpha = min(1.2*kAlpha, 1.0);                
            } else if(kP1G >= kG){
                cout<<"Oops. Error increased. To make it go faster, changing alpha from "<<kAlpha<<" to "<<0.8*kAlpha<<"..."<<endl;
                kAlpha = 0.8*kAlpha;
            }     
            for(int j = 0; j < n; j++){   
                // Determine the perturbation parameter
                kEpsilon = max(epsilon, abs(epsilon * kX0(j)));
                //kEpsilon = epsilon;
            
                // Perturb the initial conditions            
                pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(j);
                
                // Solve the perturbed initial value problem            
                integrate_const(IVPPStepper, dFunctionWrapper, pX, t0, tm, h, odePObserver);
                IVPPColumnIndex = 0;

                // Compute a column of the adjusting matrix                
                S.col(j) = (BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))- gkX0)/kEpsilon;                
                // cout<<"gkxp = "<<endl<<BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))<<endl;
            }
            // cout<<"S = "<<endl<<S<<endl;

            // Solve the linarized adjusting equation
            kX0 = S.colPivHouseholderQr().solve(-kAlpha*gkX0) + kX0;

            // Solve the initial value problem   
            kX0Temp = kX0;
            integrate_const(IVPIStepper, dFunctionWrapper, kX0Temp, t0, tm, h, odeIObserver);   
            IVPIColumnIndex = 0; 
            gkX0 = BCFunction(getBCs(IVPITSolutions, IVPIXSolutions));

            kG = kP1G;
            kP1G = gkX0.norm()/sqrt(n);
            cout<<"kP1G = "<<kP1G<<endl;
            ++k;
        }  
        cout<<"x(0) = "<<endl<<IVPIXSolutions.col(0)<<endl;
        cout<<"x(nSamples-1) = "<<endl<<IVPIXSolutions.col(nSamples-1)<<endl;
        //x0 = x0 + epsilon*MatrixXd::Identity(n,n).col(0);
        // pX = kX0 + kEpsilon*MatrixXd::Identity(n,n).col(0);
        // integrate_const(IVPPStepper, dFunctionWrapper, x0, t0, tm, h, odePObserver);
        // IVPPColumnIndex = 0; 

        // Perturb the initial conditions            
        // x0 = x0 + kEpsilon*MatrixXd::Identity(n,n).col(1);                
        // Solve the perturbed initial value problem            
        // integrate_const(IVPPStepper, dFunctionWrapper, x0, t0, tm, h, odePObserver);
        // IVPPColumnIndex = 0;
        // cout<<"Check BCs = "<<endl<<BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))<<endl;
        return 0;
    }
// ================


