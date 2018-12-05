#include<iostream>
#include<nlmp_bvp.hpp>

using namespace std;

VectorXd dFunction(double t, VectorXd x, int intervalID = 0){
***REMOVED*** VectorXd dxdt();
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = 0.01*(1 - (x(0))^2)*x(1) - x(0);
***REMOVED*** return dxdt;
}

VectorXd BCFunction(MatrixXd BC){
***REMOVED*** VectorXd residues();
***REMOVED*** residues(0) = 4*BC(0,0) + (BC(1,1))^2 - 8;
***REMOVED*** residues(1) = (BC(0,0))^2 + 2*BC(1,1) - 5;
***REMOVED*** return residues;
}

int main(int argc, char **argv)
{
***REMOVED***nlmp_bvp(1, 1, dFunction, BCFunction);
***REMOVED***cout << "Hello" << endl;
***REMOVED***return 0;
}