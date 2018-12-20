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
const double epsilon = 1e-8;
const double alpha = 15;
const double beta = 3;
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
        double kG, kP1G;
        m = tNodes.cols();  
        RowVectorXd IVPITSolutions, IVPPTSolutions;
        RowVectorXi boundaryColumns;
        VectorXd kX0;
        VectorXd pX;
        VectorXd gkX0;
        MatrixXd S;         
        MatrixXd IVPIXSolutions, IVPPXSolutions; 
        runge_kutta_dopri5<VectorXd,double,VectorXd,double,vector_space_algebra> stepper;

        t0 = tNodes(0);
        tm = tNodes(m-1);
        nSamples = floor((tm - t0)/h) + 1;
        IVPIXSolutions.resize(n, nSamples);
        IVPITSolutions.resize(1, nSamples);
        IVPPXSolutions.resize(n, nSamples);
        IVPPTSolutions.resize(1, nSamples);
        boundaryColumns.resize(1, m);
        kX0.resize(n,1);
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
            // cout<<"<Inside getBCs>"<<endl;
            // cout<<"BCIndices = "<<BCIndices<<endl;
            // cout<<"xSolutions(Eigen::all, BCIndices) = "<<xSolutions(Eigen::all, BCIndices)<<endl;
            // cout<<"</Inside getBCs>"<<endl<<endl;
            return xSolutions(Eigen::all, BCIndices);
        };

        // Solve the initial value problem for the first time  
        integrate_const(stepper, dFunctionWrapper, kX0, t0, tm, h, odeIObserver);  
        IVPIColumnIndex = 0;  
        gkX0 = BCFunction(getBCs(IVPITSolutions, IVPIXSolutions));
        kG = gkX0.norm();
        // cout<<"<IVPI solutions>"<<endl;
        // cout<<"x = "<<IVPIXSolutions<<endl;
        // cout<<"t = "<<IVPITSolutions<<endl;
        // cout<<"</IVPI solutions>"<<endl<<endl;
        do{          
            for(int j = 0; j < n; j++){   
                // Determine the perturbation parameter
                kEpsilon = max(epsilon, abs(epsilon * kX0(j,0)));
                //kEpsilon = epsilon;
            
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

                // Compute a column of the adjusting matrix                
                S.col(j) = (BCFunction(getBCs(IVPPTSolutions, IVPPXSolutions))- gkX0)/kEpsilon;                
            }

            // Solve the linarized adjusting equation
            kX0 = S.colPivHouseholderQr().solve(-gkX0) + kX0;

            // Solve the initial value problem   
            integrate_const(stepper, dFunctionWrapper, kX0, t0, tm, h, odeIObserver);   
            IVPIColumnIndex = 0; 
            gkX0 = BCFunction(getBCs(IVPITSolutions, IVPIXSolutions));
            // cout<<"<IVPI solutions>"<<endl;
            // cout<<"x = "<<IVPIXSolutions<<endl;
            // cout<<"t = "<<IVPITSolutions<<endl;
            // cout<<"</IVPI solutions>"<<endl<<endl;

            kP1G = gkX0.norm();

            if(kP1G <= pow(10, -alpha)){
                cout<<"Stopping criterion reached..."<<endl;
                cout<<"<IVPI solutions>"<<endl;
                cout<<"x(t0) = "<<IVPIXSolutions.col(0)<<endl;
                cout<<"x(tm) = "<<IVPIXSolutions.col(m-1)<<endl;
                cout<<"t = "<<tNodes<<endl;
                cout<<"</IVPI solutions>"<<endl<<endl;
                break;
            } else if(pow(10, -alpha) < kP1G && kP1G <= pow(10, -beta) && kP1G/kG < 1){                
                ++k;           
                cout<<"Stopping conditon = "<<kP1G - pow(10, -alpha)<<endl;
                cout<<"Proceeding for iteration "<<k<<"..."<<endl;
            } else{
                cout<<"Singular problem..."<<endl;
                cout<<"k = "<<endl;
                cout<<"pow(10, -alpha) = "<<pow(10, -alpha)<<endl;
                cout<<"kP1G = "<<kP1G<<endl;
                cout<<"pow(10, -beta) = "<<pow(10, -beta)<<endl;
                cout<<"kG = "<<kG<<endl;
                cout<<"kP1G/kG = "<<kP1G/kG<<endl;
                break;
            }

            kG = kP1G;
        }while(true);    
        return 0;
    }
// ================


