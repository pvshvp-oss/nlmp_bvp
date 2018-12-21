// ===============================
// Includes and global definitions
// ===============================
#include <cmath>                                             // C++ analog to math.h
#include <Eigen/Eigen>                                       // For matrix and vector math
#include <boost/numeric/odeint.hpp>                          // For the initial value problem (IVP) solver
#include <nlmp_bvp.hpp>                                      // For the function declarations   
using namespace std;                                         //
using namespace Eigen;                                       //
using namespace boost::numeric::odeint;                      //
using RowVectorXi = Matrix<int, 1, Dynamic>;                 // For the convenience of declaring row vectors
using RowVectorXd = Matrix<double, 1, Dynamic>;              // For the convenience of declaring row vectors
using StepperType = runge_kutta_dopri5<                      // For the convenience of specifying the stepper type for the IVP solver
                                        VectorXd,            // the state vector type for the IVP solver
                                        double,              // the state variable value type for the IVP solver
                                        VectorXd,            // the type for the derivative of the state vector x for the IVP solver
                                        double,              // the type for independent variable t for the IVP solver
                                        vector_space_algebra // the type of algebra to be done for the IVP solver
                                        >;                   // 
const double INF = std::numeric_limits<double>::infinity();  // INF = infinity                                        
// ===============================

// ==================
// Function "nlmp_bvp"
// ===================
BVPSolution nlmp_bvp(
    int n,                                 // n              = the number of differential equations = the number of boundary conditions
    int m,                                 // m              = the number of nodes at which boundary conditions are specified
    int nGrid,                             // nGrid          = the number of points at which the state is evaluated
    RowVectorXd t_BC,                      // t_BC           = row vector of values at which the boundary conditions are specified           -- (1xm)
    VectorXd _0_x_t1,                      // _0_x_t1        = column vector of the guessed initial state                                    -- (nx1)    
    VectorXd dxBydt(double t, VectorXd x), // dxBydt         = a function that defines the derivative of a state vector x at t               -- (nx1)
    VectorXd BCResidues(MatrixXd x_BC),    // BCResidues     = a function that defines the boundary condition residues at state vectors x_BC -- (nx1) 
    IVAMParameters ivamParameters          // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)
    ){  

        // Variable declarations        
        int j;                      // j            = the inner iterating variable for IVAM                                                              -- [0,n-1]
        int k;                      // k            = the outer iterating variable for IVAM                                                              -- [0,Inf)
        int iCol;                   // iCol         = the column index of the x solution for the IVP solver                                              -- [0,nGrid-1]
        int iColP;                  // iColP        = the column index of the perturbed x solution for the IVP solver   
        double h;                   // h            = the stepper step size for the IVP solver                                 -- [0,nGrid-1] 
        double t0;                  // t0           = the first boundary value of the independent variable
        double tm;                  // tm           = the last boundary value of the independent variable
        double _k_epsilon_j;        // _k_epsilon_j = the perturbation parameter for a state variable at a particular iteration k
        double _k_alpha;            // _k_alpha     = the relaxation factor at a particular iteration k to scale the adjustment to the initial condition 
        double _k_G;                // _k_G         = the Root Mean Square (RMS) error of boundary residues at a particular iteration k
        double _k_GPrev;            // _k_GPrev     = the Root Mean Square (RMS) error of boundary residues at the previous iteration k-1
        RowVectorXd tSol(nGrid);    // tSol         = the independent variable t over the whole grid in the solution of the IVP solver                   -- (1xnGrid)
        MatrixXd xSol(n,nGrid);     // xSol         = the state vector x integrated over the whole grid in the solution of the IVP solver                -- (nxnGrid)
        RowVectorXd tsolP(nGrid);   // tsolP        = the independent variable t over the whole grid in the perturbed solution of the IVP solver         -- (1xnGrid)    
        MatrixXd xSolP(n,nGrid);    // xSolP        = the state vector x integrated over the whole grid in the perturbed solution of the IVP solver      -- (nxnGrid)
        RowVectorXi BCCols(m);      // BCCols       = the columns in the grid that correspond to boundary values                                         -- (1xm)
        VectorXd _k_x_t1(n);        // _k_x_t1      = the computed initial state vector in the k-th iteration                                            -- (nx1)
        VectorXd _k_x_t1P(n);       // _k_x_t1      = the computed perturbed initial state vector in the k-th iteration                                  -- (nx1)
        VectorXd x_t1(n);           // x_t1         = the computed initial state vector to be input to the IVP solver                                    -- (nx1)
        VectorXd x_t1P(n);          // x_t1P        = the computed perturbed initial state vector to be input to the IVP solver                          -- (nx1)
        VectorXd _k_g(n);           // _k_g         = the boundary condition residues in the k-th iteration                                              -- (nx1)
        VectorXd _k_g_j(n);         // _k_g_j       = the j-th boundary condition perturbed system residues in the k-th iteration                        -- (nx1)
        MatrixXd _k_S(n,n);         // _k_S         = the adjusting matrix for correcting the initial condition k-th iteration                           -- (nxn) 
        BVPSolution bvpSolution;

        // Variable definitions
        t0     = t_BC(0);
        tm     = t_BC(m-1);
        h      = (tm - t0)/(nGrid-1);
        BCCols = ((t_BC-t0*RowVectorXd::Ones(m))/h).array().round().cast<int>();

        cout<<"BCCols = "<<BCCols<<endl;

        // Wrapper function to be called by the IVP solver to retrieve the definitions for the differential equations
        auto dxBydtWrapper = [dxBydt] // Captured variables
                             (const VectorXd &x, VectorXd &dxdt, double t){
            dxdt = dxBydt(t, x);
        };

        // Observer to store the solutions of the IVP Solver
        auto storeSol = [&iCol, &tSol, &xSol] // Captured variables
                        (const VectorXd &x , const double t){
            tSol(iCol)     = t;                    
            xSol.col(iCol) = x;            
            ++iCol;            
        };

        // Observer to store the perturbed solutions of the IVP Solver
        auto storeSolP = [&iColP, &tsolP, &xSolP] // Captured variables
                         (const VectorXd &x , const double t){
            tsolP(iColP)     = t;                    
            xSolP.col(iColP) = x;            
            ++iColP;            
        };       
 
        _k_GPrev = INF;
        _k_alpha = ivamParameters.ALPHA;
        k        = 0;       // Set k for the first time
        _k_x_t1  = _0_x_t1; // Assign the initial condition state vector for the 0-th iteration

        // Solve the initial value problem for the first time 
        x_t1     = _k_x_t1; // Assign the current initial condition state vector to a dummy variable
        iCol     = 0;       // Set the solution column index to 0 before the IVP solver starts integrating
        integrate_const(StepperType(), dxBydtWrapper, x_t1, t0, tm, h, storeSol); 
        _k_g = BCResidues(xSol(Eigen::all, BCCols));
        _k_G = _k_g.norm()/sqrt(n);
        while(_k_G > ivamParameters.SIGMA){ 

            cout<<"k = "<<k<<"... _k_G = "<<_k_G<<"... _k_GPrev = "<<_k_GPrev<<endl;

            if(_k_G < 0.1*_k_GPrev) {
                cout<<"Going too fast. Changing alpha from "<<_k_alpha<<" to "<<fmin(1.2*_k_alpha, 1.0)<<"..."<<endl;
                _k_alpha = fmin(1.2*_k_alpha, 1.0);                
            } else if(_k_G >= _k_GPrev){
                cout<<"Error increased. Changing alpha from "<<_k_alpha<<" to "<<0.8*_k_alpha<<"..."<<endl;
                _k_alpha = 0.8*_k_alpha;
            }     
            for(j = 0; j < n; j++){   
                // Determine the perturbation parameter
                _k_epsilon_j = fmax(ivamParameters.EPSILON, fabs(ivamParameters.EPSILON * _k_x_t1(j)));
            
                // Perturb the initial conditions            
                _k_x_t1P = _k_x_t1 + _k_epsilon_j*MatrixXd::Identity(n,n).col(j);
                
                // Solve the perturbed initial value problem   
                x_t1P = _k_x_t1P; // Assign the current initial condition state vector to a dummy variable    
                iColP  = 0;        // Set the solution column index to 0 before the IVP solver starts integrating     
                integrate_const(StepperType(), dxBydtWrapper, x_t1P, t0, tm, h, storeSolP); 
                _k_g_j =  BCResidues(xSolP(Eigen::all, BCCols));

                // Compute a column of the adjusting matrix                
                _k_S.col(j) = (_k_g_j- _k_g)/_k_epsilon_j;                
            }
                        
            _k_alpha = 1;
            VectorXd xPrev = _k_x_t1;
            VectorXd gPrev = _k_g;

            // Solve the linarized adjusting equation
            _k_x_t1 = _k_S.colPivHouseholderQr().solve(-_k_alpha*_k_g) + _k_x_t1;
            _k_GPrev = _k_G;
            ++k;

            // Solve the initial value problem   
            x_t1     = _k_x_t1; // Assign the current initial condition state vector to a dummy variable
            iCol     = 0;       // Set the solution column index to 0 before the IVP solver starts integrating
            integrate_const(StepperType(), dxBydtWrapper, x_t1, t0, tm, h, storeSol); 
            _k_g = BCResidues(xSol(Eigen::all, BCCols));
            _k_G = _k_g.norm()/sqrt(n);      

            cout<<"dRes/dx = "<<endl<<_k_S<<endl;
            cout<<"dx = "<<endl<<_k_x_t1 - xPrev<<endl;
            cout<<"Res = "<<endl<<_k_g<<endl;

            if(k >= 9){
                break;
            }
        }  
        bvpSolution.t    = tSol;
        bvpSolution.x    = xSol;
        bvpSolution.t_BC = t_BC;
        bvpSolution.x_BC = xSol(Eigen::all, BCCols);
        return bvpSolution;
    }
// ===================


