// ========================================
// Author: Shivanand Pattanshetti 
// Email : shivanand.pattanshetti@gmail.com
// ========================================

// ===============================
// Includes and global definitions
// ===============================
#include <cmath>                                                        // C++ analog to math.h
#include <limits>                                                       // For NaN
#include <boost/numeric/odeint.hpp>                                     // For the initial value problem (IVP) solver
#include "ivamparameters.hpp"                                           // For the struct "IVAMParameters"
#include "bvpsolution.hpp"                                              // For the struct "BVPSolution"
using namespace std;                                                    // For cout
using namespace Eigen;                                                  // For matrix and vector data types and operations
using namespace boost::numeric::odeint;                                 // For the initial value problem (IVP) solver     
template <typename T> using StepperType = runge_kutta_dopri5<           // For the convenience of specifying the stepper type for the IVP solver
                                        VectorXm<T>,                    // the state vector type for the IVP solver
                                        T,                              // the state variable value type for the IVP solver
                                        VectorXm<T>,                    // the type for the derivative of the state vector x for the IVP solver
                                        T,                              // the type for independent variable t for the IVP solver
                                        vector_space_algebra            // the type of algebra to be done for the IVP solver
                                        >;                              // 
template <typename T> const T INF = std::numeric_limits<T>::infinity(); // INF = infinity                         
// ===============================

// ==================
// Function "nlmpBVP"
// ==================
template <typename T> BVPSolution<T> nlmpBVP(
    int n,                                   // n              = the number of differential equations = the number of boundary conditions
    int m,                                   // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                               // nGrid          = the number of points at which the state is evaluated
    RowVectorXm<T> tBC,                      // tBC            = row vector of values at which the boundary conditions are specified           -- (1xm)
    VectorXm<T> oxt1,                        // oxt1           = column vector of the guessed initial state                                    -- (nx1)    
    VectorXm<T> dxBydt(T t, VectorXm<T> x),  // dxBydt         = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXm<T> BCResidues(MatrixXm<T> xBC), // BCResidues     = a function that defines the boundary condition residues at state vectors xBC  -- (nx1) 
    IVAMParameters<T> ivamParameters         // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    ){  

        // Variable declarations  
        int j;                          // j         = the inner iterating variable for IVAM                                                              -- [0,n-1]
        int k;                          // k         = the outer iterating variable for IVAM                                                              -- [0,Inf)
        int col;                        // col       = the column index of the x solution for the IVP solver                                              -- [0,nGrid-1]
        int colP;                       // colP      = the column index of the perturbed x solution for the IVP solver                                    -- [0,nGrid-1]
        T h;                            // h         = the stepper step size for the IVP solver                                                            
        T t0;                           // t0        = the first boundary value of the independent variable
        T tm;                           // tm        = the last boundary value of the independent variable
        T kepsilonj;                    // kepsilonj = the perturbation parameter for a state variable at a particular iteration k
        T kalpha;                       // kalpha    = the relaxation factor at a particular iteration k to scale the adjustment to the initial condition 
        T kG;                           // kG        = the Root Mean Square (RMS) error of boundary residues at a particular iteration k
        T kGPrev;                       // kGPrev    = the Root Mean Square (RMS) error of boundary residues at the previous iteration k-1
        RowVectorXm<T> tSol(nGrid);     // tSol      = the independent variable t over the whole grid in the solution of the IVP solver                   -- (1xnGrid)
        MatrixXm<T> xSol(n,nGrid);      // xSol      = the state vector x integrated over the whole grid in the solution of the IVP solver                -- (nxnGrid)
        RowVectorXm<T> tSolP(nGrid);    // tSolP     = the independent variable t over the whole grid in the perturbed solution of the IVP solver         -- (1xnGrid)    
        MatrixXm<T> xSolP(n,nGrid);     // xSolP     = the state vector x integrated over the whole grid in the perturbed solution of the IVP solver      -- (nxnGrid)
        RowVectorXi BCCols(m);          // BCCols    = the columns in the grid that correspond to boundary values                                         -- (1xm)
        VectorXm<T> kxt1(n);            // kxt1      = the computed initial state vector in the k-th iteration                                            -- (nx1)
        VectorXm<T> kxt1P(n);           // kxt1      = the computed perturbed initial state vector in the k-th iteration                                  -- (nx1)
        VectorXm<T> xt(n);              // xt        = the computed initial state vector to be input to the IVP solver                                    -- (nx1)
        VectorXm<T> xtP(n);             // xtP       = the computed perturbed initial state vector to be input to the IVP solver                          -- (nx1)
        VectorXm<T> kg(n);              // kg        = the boundary condition residues in the k-th iteration                                              -- (nx1)
        VectorXm<T> kgj(n);             // kgj       = the j-th boundary condition perturbed system residues in the k-th iteration                        -- (nx1)
        MatrixXm<T> kS(n,n);            // kS        = the adjusting matrix for correcting the initial condition k-th iteration                           -- (nxn) 
        BVPSolution<T> bvpSolution;         

        // Variable definitions
        t0     = tBC(0);             
        tm     = tBC(m-1);
        h      = (tm - t0)/(nGrid-1);
        BCCols = ((tBC-t0*RowVectorXm<T>::Ones(m))/h).template array().template round().template cast<int>();

        if(ivamParameters.printDebug){
            cout<<"Boundary nodes correspond to the below columns: "<<endl<<BCCols<<endl<<endl;
        }
        
        // Wrapper function to be called by the IVP solver to retrieve the definitions for the differential equations
        auto dxBydtWrapper = [dxBydt] // Captured variables
                             (const VectorXm<T> &x, VectorXm<T> &dxdt, T t){
            dxdt = dxBydt(t, x);
        };

        // Observer to store the solutions of the IVP Solver
        auto storeSol = [&col, &tSol, &xSol] // Captured variables
                        (const VectorXm<T> &x , const T t){
            tSol(col)     = t;                    
            xSol.col(col) = x;            
            ++col;            
        };

        // Observer to store the perturbed solutions of the IVP Solver
        auto storeSolP = [&colP, &tSolP, &xSolP] // Captured variables
                         (const VectorXm<T> &x , const T t){
            tSolP(colP)     = t;                    
            xSolP.col(colP) = x;            
            ++colP;            
        };       
 
        // Initialize the required variables
        kGPrev   = INF<T>;
        kalpha   = ivamParameters.ALPHA;
        k        = 0;       
        kxt1     = oxt1; 

        // Solve the initial value problem for the first time 
        xt      = kxt1; // Assign the current initial condition state vector to a dummy variable
        col     = 0;    // Set the solution column index to 0 before the IVP solver starts integrating
        integrate_const(StepperType<T>(), dxBydtWrapper, xt, t0, tm, h, storeSol); 
        kg = BCResidues(xSol(Eigen::all, BCCols));
        kG = kg.norm()/sqrt(n);

        while(kG > ivamParameters.SIGMA /* When the error is more than the max. threshold */){      

            // Adjust the relaxation parameter to control the rate of convergence
            if(kG < 0.01*kGPrev) {
                if(ivamParameters.printDebug){
                    cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Increasing alpha..."<<endl;
                }
                kalpha = fmin(1.2*kalpha, 1.0); 
            } else if(kG >= kGPrev){
                if(ivamParameters.printDebug){
                    cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Decreasing alpha..."<<endl;             
                }
                kalpha = 0.8*kalpha;        
            } else{
                if(ivamParameters.printDebug){
                    cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<endl;
                }
            }

            // Inner loop to perturb each state variable separately and find the normalized change in the boundary condition residues
            for(j = 0; j < n; j++){    
                // Determine the perturbation parameter
                kepsilonj = fmax(ivamParameters.EPSILON, fabs(ivamParameters.EPSILON * kxt1(j)));
                // kepsilonj = ivamParameters.EPSILON;
            
                // Perturb the initial conditions   
                kxt1P = kxt1 + kepsilonj*MatrixXm<T>::Identity(n,n).col(j);
                
                // Solve the perturbed initial value problem   
                xtP = kxt1P; // Assign the current initial condition state vector to a dummy variable    
                colP  = 0;   // Set the solution column index to 0 before the IVP solver starts integrating     
                integrate_const(StepperType<T>(), dxBydtWrapper, xtP, t0, tm, h, storeSolP); 
                kgj =  BCResidues(xSolP(Eigen::all, BCCols));

                // Compute one column of the adjusting matrix
                // Each column corresponds to the change in boundary condition residues due to a corresponding perturbation in *one* state variable, normalized with respect to the perturbation
                kS.col(j) = (kgj - kg)*pow(kepsilonj,-1);  
            }

            kalpha = 1;
            // Solve the linarized adjusting equation            
            kxt1 = kxt1 - kS.colPivHouseholderQr().solve(kalpha*kg);

            // Start the next iteration
            kGPrev = kG;
            ++k;

            // Solve the initial value problem   
            xt      = kxt1; // Assign the current initial condition state vector to a dummy variable
            col     = 0;   // Set the solution column index to 0 before the IVP solver starts integrating
            integrate_const(StepperType<T>(), dxBydtWrapper, xt, t0, tm, h, storeSol); 
            kg = BCResidues(xSol(Eigen::all, BCCols));
            kG = kg.norm()/sqrt(n);

            if(k >= 1000){
                if(ivamParameters.printDebug){
                    cout<<"[WARNING]: The solution did not converge after 1000 iterations. Terminating the process."<<endl;
                }
                break;
            }
        }  

        if(ivamParameters.printDebug){
            cout<<"Ran "<<k<<" iteration(s)."<<endl;
        }

        bvpSolution.t   = tSol;
        bvpSolution.x   = xSol;
        bvpSolution.tBC = tBC;
        bvpSolution.xBC = xSol(Eigen::all, BCCols);
        return bvpSolution;
    }
// ==================

// ===================
// Function "nlmpBVP2"
// ===================
template <typename T> BVPSolution<T> nlmpBVP2(
    int n,                                     // n              = the number of differential equations = the number of boundary conditions
    int m,                                     // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                                 // nGrid          = the number of points at which the state can be evaluated
    RowVectorXm<T> tBC,                        // tBC            = row vector of values at which boundary conditions are specified               -- (1xm)
    VectorXm<T> oxt1,                          // oxt1           = matrix of the guessed initial state                                           -- (nx(m-1))    
    VectorXm<T> dxBydt(T t, VectorXm<T> x),    // dxBydt         = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXm<T> BCResidues(MatrixXm<T> xBCL,   // BCResidues     = a function that defines the boundary condition residues...
                           MatrixXm<T> xBCR),  //                  ...at the left and right state vectors xBCL and xBCR                          -- (n(m-1)x1)      
    IVAMParameters<T> ivamParameters           // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    ){  

        // Variable declarations  
        int i;                           // i         = the iterating variable that keeps track of the independent variable node                                -- [0,m-1]
        int j;                           // j         = the inner iterating variable for IVAM                                                                   -- [0,n-1]
        int k;                           // k         = the outer iterating variable for IVAM                                                                   -- [0,Inf)
        int l;
        int col;                         // col       = the column index of the x solution for the IVP solver                                                   -- [0,nGrid-1]
        int colP;                        // colP      = the column index of the perturbed x solution for the IVP solver                                         -- [0,nGrid-1]
        T h;                             // h         = the stepper step size for the IVP solver                                                            
        T kepsilonj;                     // kepsilonj = the perturbation parameter for a state variable at a particular iteration k
        T kalpha;                        // kalpha    = the relaxation factor at a particular iteration k to scale the adjustment to the initial condition 
        T kG;                            // kG        = the Root Mean Square (RMS) error of boundary residues at a particular iteration k
        T kGPrev;                        // kGPrev    = the Root Mean Square (RMS) error of boundary residues at the previous iteration k-1
        RowVectorXm<T> tSol(nGrid+m-2);  // tSol      = the independent variable t over the whole grid in the solution of the IVP solver                        -- (1x(nGrid+m-2))
        MatrixXm<T> xSol(n,nGrid+m-2);   // xSol      = the state vector x integrated over the whole grid in the solution of the IVP solver                     -- (nx(nGrid+m-2)m)
        RowVectorXm<T> tSolP(nGrid);     // tSolP     = the independent variable t over the whole grid in the perturbed solution of the IVP solver              -- (1xnGrid)    
        MatrixXm<T> xSolP(n,nGrid);      // xSolP     = the state vector x integrated over the whole grid in the perturbed solution of the IVP solver           -- (nxnGrid)
        RowVectorXi BCCols(m);           // BCCols    = the columns in the grid that correspond to boundary values                                         -- (1xm)
        MatrixXm<T> kxt1(n,m-1);         // kxt1      = the computed initial state vector in the k-th iteration at the left side of every interval              -- (nx(m-1))
        MatrixXm<T> kxtm(n,m-1);         // kxtm      = the computed final state vector in the k-th iteration at the left side of every interval                -- (nx(m-1))
        MatrixXm<T> kxt1P(n,m-1);        // kxt1      = the computed perturbed initial state vector in the k-th iteration                                       -- (nx(m-1))
        MatrixXm<T> kxtmP(n,m-1);        // kxtm      = the computed perturbed initial state vector in the k-th iteration                                       -- (nx(m-1))
        VectorXm<T> xt(n);               // xt        = the computed state vector to be input to the IVP solver                                                 -- (nx1)
        VectorXm<T> xtP(n);              // xtP       = the computed perturbed state vector to be input to the IVP solver                                       -- (nx1)
        VectorXm<T> kg(n*(m-1));         // kg        = the boundary condition residues in the k-th iteration                                                   -- ((n(m-1))x1)
        VectorXm<T> kgj(n*(m-1));       // kgj       = the j-th boundary condition perturbed system residues in the k-th iteration                             -- ((n(m-1))x1)
        MatrixXm<T> kS(n*(m-1),n*(m-1)); // kS        = the adjusting matrix for correcting the initial condition k-th iteration                                -- ((n(m-1))x(n(m-1))) 
        BVPSolution<T> bvpSolution;

        // Variable definitions
        h = (tBC(tBC.cols()-1) - tBC(0))/(nGrid-1);
        BCCols = ((tBC-tBC(0)*RowVectorXm<T>::Ones(m))/h).template array().template round().template cast<int>();

        if(ivamParameters.printDebug){
            cout<<"Boundary nodes correspond to the below columns: "<<endl<<BCCols<<endl<<endl;
        }
        
        // Wrapper function to be called by the IVP solver to retrieve the definitions for the differential equations
        auto dxBydtWrapper = [dxBydt] // Captured variables
                             (const VectorXm<T> &x, VectorXm<T> &dxdt, T t){
            dxdt = dxBydt(t, x);
        };

        // Observer to store the solutions of the IVP Solver
        auto storeSol = [&col, &tSol, &xSol] // Captured variables
                        (const VectorXm<T> &x , const T t){
            tSol(col)     = t;                    
            xSol.col(col) = x;            
            ++col;            
        };

        // Observer to store the perturbed solutions of the IVP Solver
        auto storeSolP = [&colP, &tSolP, &xSolP] // Captured variables
                         (const VectorXm<T> &x , const T t){
            tSolP(colP)     = t;                    
            xSolP.col(colP) = x;            
            ++colP;            
        };       
 
        // Initialize the required variables
        kGPrev         = INF<T>;
        kalpha         = ivamParameters.ALPHA;
        k              = 0;       
        kxt1           = oxt1;

        col = 0;
        for(i=0; i<(m-1); i++){
            xt  = kxt1.col(i);
            integrate_const(StepperType<T>(), dxBydtWrapper, xt, tBC(i), tBC(i+1), h, storeSol); 
            kxtm.col(i) = xt;
        }
        kg = BCResidues(kxt1,kxtm);
        kG = kg.norm()/sqrt(n*(m-1));

        while(kG > ivamParameters.SIGMA /* When the error is more than the max. threshold */){      

            // Adjust the relaxation parameter to control the rate of convergence
            if(kG < 0.01*kGPrev) {
                if(ivamParameters.printDebug){
                    cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Increasing alpha..."<<endl;
                }
                kalpha = fmin(1.2*kalpha, 1.0); 
            } else if(kG >= kGPrev){
                if(ivamParameters.printDebug){
                    cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<" Decreasing alpha..."<<endl;             
                }
                kalpha = 0.8*kalpha;        
            } else{
                if(ivamParameters.printDebug){
                    cout<<"k = "<<k<<"... kG = "<<kG<<"... kGPrev = "<<kGPrev<<" alpha = "<<kalpha<<endl;
                }
            }

            for(l = 0; l<(m-1); l++){  
                kxt1P = kxt1;
                for(j=0; j<n; i++){   
                    kepsilonj = fmax(ivamParameters.EPSILON, fabs(ivamParameters.EPSILON * kxt1(j,i)));
                    // kepsilonj = ivamParameters.EPSILON;        
                    kxt1P.col(l) = kxt1.col(l) + kepsilonj*MatrixXm<T>::Identity(n,n).col(j);     
                    colP = 0;               
                    for(i=0; i<(m-1); i++){
                        xtP  = kxt1P.col(i);
                        integrate_const(StepperType<T>(), dxBydtWrapper, xtP, tBC(i), tBC(i+1), h, storeSolP); 
                        kxtmP.col(i) = xtP;
                    }             
                    kgj = BCResidues(kxt1P,kxtmP);
                    kS.col(l*n+j) = (kgj - kg)/kepsilonj;
                }              
            }

            kalpha = 1;
            // Solve the linarized adjusting equation            
            kxt1 = kxt1 - Map<MatrixXm<T>>(kS.colPivHouseholderQr().solve(kalpha*kg).data(),n,m-1);

            // Start the next iteration
            kGPrev = kG;
            ++k;

            // Solve the initial value problems
            col = 0;
            for(i=0; i<(m-1); i++){
                xt  = kxt1.col(i);
                integrate_const(StepperType<T>(), dxBydtWrapper, xt, tBC(i), tBC(i+1), h, storeSol); 
                kxtm.col(i) = xt;
            }
            kg = BCResidues(kxt1,kxtm);
            kG = kg.norm()/sqrt(n*(m-1));

            if(k >= 1000){
                if(ivamParameters.printDebug){
                    cout<<"[WARNING]: The solution did not converge after 1000 iterations. Terminating the process."<<endl;
                }
                break;
            }
        }  

        if(ivamParameters.printDebug){
            cout<<"Ran "<<k<<" iteration(s)."<<endl;
        }

        bvpSolution.t   = tSol;
        bvpSolution.x   = xSol;
        bvpSolution.tBC = tBC;
        bvpSolution.xBC = xSol(Eigen::all, BCCols);
        return bvpSolution;
    }
// ==================