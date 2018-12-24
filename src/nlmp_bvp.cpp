// ========================================
// Author: Shivanand Pattanshetti 
// Email : shivanand.pattanshetti@gmail.com
// ========================================

// ===============================
// Includes and global definitions
// ===============================
#include <cmath>                                             // C++ analog to math.h
#include <boost/numeric/odeint.hpp>                          // For the initial value problem (IVP) solver
#include <Eigen/Eigen>                                       // For matrix and vector math
#include <Eigen/MPRealSupport>                               // For arbitrary precision computation
#include <nlmp_bvp.hpp>                                      // For the function declarations  
using namespace mpfr;                                        // For arbitrary precision computation
using namespace std;                                         // For cout
using namespace Eigen;                                       // For matrix and vector data types and operations
using namespace boost::numeric::odeint;                      // For the initial value problem (IVP) solver
using MatrixXm = Matrix<mpreal, Dynamic, Dynamic>;           // Dynamic size matrix with the 'mpreal' data-type
using VectorXm = Matrix<mpreal, Dynamic, 1>;                 // Dynamic size vector with the 'mpreal' data-type
using RowVectorXi = Matrix<int, 1, Dynamic>;                 // Dynamic size row vector with the 'integer' data-type
using RowVectorXm = Matrix<mpreal, 1, Dynamic>;              // Dynamic size row vector with the 'mpreal' data-type
using StepperType = runge_kutta_dopri5<                      // For the convenience of specifying the stepper type for the IVP solver
                                        VectorXm,            // the state vector type for the IVP solver
                                        mpreal,              // the state variable value type for the IVP solver
                                        VectorXm,            // the type for the derivative of the state vector x for the IVP solver
                                        mpreal,              // the type for independent variable t for the IVP solver
                                        vector_space_algebra // the type of algebra to be done for the IVP solver
                                        >;                   // 
const mpreal INF = std::numeric_limits<mpreal>::infinity();  // INF = infinity                                        
// ===============================

// ==================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
    int n,                                 // n              = the number of differential equations = the number of boundary conditions
    int m,                                 // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                             // nGrid          = the number of points at which the state is evaluated
    RowVectorXm tBC,                       // tBC            = row vector of values at which the boundary conditions are specified           -- (1xm)
    VectorXm oxt1,                         // oxt1           = column vector of the guessed initial state                                    -- (nx1)    
    VectorXm dxBydt(mpreal t, VectorXm x), // dxBydt         = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXm BCResidues(MatrixXm xBC),     // BCResidues     = a function that defines the boundary condition residues at state vectors xBC  -- (nx1) 
    IVAMParameters ivamParameters          // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    ){  

        mpreal::set_default_prec(256);

        // Variable declarations  
        int omega = 0.5; 
        int j;                       // j         = the inner iterating variable for IVAM                                                              -- [0,n-1]
        int k;                       // k         = the outer iterating variable for IVAM                                                              -- [0,Inf)
        int iCol;                    // iCol      = the column index of the x solution for the IVP solver                                              -- [0,nGrid-1]
        int iColP;                   // iColP     = the column index of the perturbed x solution for the IVP solver   
        mpreal h;                    // h         = the stepper step size for the IVP solver                                                           -- [0,nGrid-1] 
        mpreal t0;                   // t0        = the first boundary value of the independent variable
        mpreal tm;                   // tm        = the last boundary value of the independent variable
        mpreal kepsilonj;            // kepsilonj = the perturbation parameter for a state variable at a particular iteration k
        mpreal kalpha;               // kalpha    = the relaxation factor at a particular iteration k to scale the adjustment to the initial condition 
        mpreal kG;                   // kG        = the Root Mean Square (RMS) error of boundary residues at a particular iteration k
        mpreal kGPrev;               // kGPrev    = the Root Mean Square (RMS) error of boundary residues at the previous iteration k-1
        RowVectorXm tSol(nGrid);     // tSol      = the independent variable t over the whole grid in the solution of the IVP solver                   -- (1xnGrid)
        MatrixXm xSol(n,nGrid);      // xSol      = the state vector x integrated over the whole grid in the solution of the IVP solver                -- (nxnGrid)
        RowVectorXm tSolPert(nGrid); // tSolPert  = the independent variable t over the whole grid in the perturbed solution of the IVP solver         -- (1xnGrid)    
        MatrixXm xSolPert(n,nGrid);  // xSolPert  = the state vector x integrated over the whole grid in the perturbed solution of the IVP solver      -- (nxnGrid)
        RowVectorXi BCCols(m);       // BCCols    = the columns in the grid that correspond to boundary values                                         -- (1xm)
        VectorXm kxt1(n);            // kxt1      = the computed initial state vector in the k-th iteration                                            -- (nx1)
        VectorXm kxt1Prev(n);        // kxt1Prev  = the computed initial state vector in the previous (k-1)-th iteration                               -- (nx1)
        VectorXm kxt1P(n);           // kxt1      = the computed perturbed initial state vector in the k-th iteration                                  -- (nx1)
        VectorXm xt1(n);             // xt1       = the computed initial state vector to be input to the IVP solver                                    -- (nx1)
        VectorXm xt1P(n);            // xt1P      = the computed perturbed initial state vector to be input to the IVP solver                          -- (nx1)
        VectorXm kg(n);              // kg        = the boundary condition residues in the k-th iteration                                              -- (nx1)
        VectorXm kgj(n);             // kgj       = the j-th boundary condition perturbed system residues in the k-th iteration                        -- (nx1)
        MatrixXm kS(n,n);            // kS        = the adjusting matrix for correcting the initial condition k-th iteration                           -- (nxn) 
        BVPSolution bvpSolution;

        // Variable definitions
        t0     = tBC(0);             
        tm     = tBC(m-1);
        h      = (tm - t0)/(nGrid-1);
        BCCols = ((tBC-t0*RowVectorXm::Ones(m))/h).array().round().cast<int>();

        cout<<"Boundary nodes correspond to the below columns: "<<endl<<BCCols<<endl<<endl;
        
        // Wrapper function to be called by the IVP solver to retrieve the definitions for the differential equations
        auto dxBydtWrapper = [dxBydt] // Captured variables
                             (const VectorXm &x, VectorXm &dxdt, mpreal t){
            dxdt = dxBydt(t, x);
        };

        // Observer to store the solutions of the IVP Solver
        auto storeSol = [&iCol, &tSol, &xSol] // Captured variables
                        (const VectorXm &x , const mpreal t){
            tSol(iCol)     = t;                    
            xSol.col(iCol) = x;            
            ++iCol;            
        };

        // Observer to store the perturbed solutions of the IVP Solver
        auto storeSolP = [&iColP, &tSolPert, &xSolPert] // Captured variables
                         (const VectorXm &x , const mpreal t){
            tSolPert(iColP)     = t;                    
            xSolPert.col(iColP) = x;            
            ++iColP;            
        };       
 
        // Initialize the required variables
        kGPrev   = INF;
        kalpha   = ivamParameters.ALPHA;
        k        = 0;       
        kxt1     = oxt1; 
        kxt1Prev = kxt1;

        // Solve the initial value problem for the first time 
        xt1      = kxt1; // Assign the current initial condition state vector to a dummy variable
        iCol     = 0;    // Set the solution column index to 0 before the IVP solver starts integrating
        integrate_const(StepperType(), dxBydtWrapper, xt1, t0, tm, h, storeSol); 
        kg = BCResidues(xSol(Eigen::all, BCCols));
        kG = kg.norm()/sqrt(n);

        while(kG > ivamParameters.SIGMA /* When the error is more than the max. threshold */){      

            // Adjust the relaxation parameter to control the rate of convergence
            if(kG < 0.01*kGPrev) {
                cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Increasing alpha..."<<endl;
                kalpha = fmin(1.2*kalpha, 1.0); 
            } else if(kG >= kGPrev){
                cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Decreasing alpha..."<<endl;             
                kalpha = 0.8*kalpha;        
            } else{
                cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<endl;
            }

            // Inner loop to perturb each state variable separately and find the normalized change in the boundary condition residues
            for(j = 0; j < n; j++){    
                // Determine the perturbation parameter
                kepsilonj = fmax(ivamParameters.EPSILON, fabs(ivamParameters.EPSILON * kxt1(j)));
                // kepsilonj = ivamParameters.EPSILON;
            
                // Perturb the initial conditions   
                kxt1P = kxt1 + kepsilonj*MatrixXm::Identity(n,n).col(j);
                
                // Solve the perturbed initial value problem   
                xt1P = kxt1P; // Assign the current initial condition state vector to a dummy variable    
                iColP  = 0;   // Set the solution column index to 0 before the IVP solver starts integrating     
                integrate_const(StepperType(), dxBydtWrapper, xt1P, t0, tm, h, storeSolP); 
                kgj =  BCResidues(xSolPert(Eigen::all, BCCols));

                // Compute one column of the adjusting matrix
                // Each column corresponds to the change in boundary condition residues due to a corresponding perturbation in *one* state variable, normalized with respect to the perturbation
                kS.col(j) = (kgj - kg)*pow(kepsilonj,-1);  
            }

            kalpha = 1;
            // Solve the linarized adjusting equation            
            kxt1Prev = kxt1; 
            kxt1 = kxt1 - kS.colPivHouseholderQr().solve(kalpha*kg);

            // Start the next iteration
            kGPrev = kG;
            ++k;

            // Solve the initial value problem   
            xt1     = kxt1; // Assign the current initial condition state vector to a dummy variable
            iCol     = 0;   // Set the solution column index to 0 before the IVP solver starts integrating
            integrate_const(StepperType(), dxBydtWrapper, xt1, t0, tm, h, storeSol); 
            kg = BCResidues(xSol(Eigen::all, BCCols));
            kG = kg.norm()/sqrt(n);

            if(k >= 1000){
                cout<<"[WARNING]: The solution did not converge after 1000 iterations. Terminating the process."<<endl;
                break;
            }
        }  

        cout<<"Ran "<<k<<" iteration(s)."<<endl;

        bvpSolution.t   = tSol;
        bvpSolution.x   = xSol;
        bvpSolution.tBC = tBC;
        bvpSolution.xBC = xSol(Eigen::all, BCCols);
        return bvpSolution;
    }
// ===================
