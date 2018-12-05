#include<iostream>
#include<nlmp_bvp.hpp>

using namespace std;

VectorXd dFunction(double t, VectorXd x, int intervalID = 0){
    VectorXd dxdt();
    dxdt(0) = x(1);
    dxdt(1) = 0.01*(1 - (x(0))^2)*x(1) - x(0);
    return dxdt;
}

VectorXd BCFunction(MatrixXd BC){
    VectorXd residues();
    residues(0) = 4*BC(0,0) + (BC(1,1))^2 - 8;
    residues(1) = (BC(0,0))^2 + 2*BC(1,1) - 5;
    return residues;
}

int main(int argc, char **argv)
{
   nlmp_bvp(1, 1, dFunction, BCFunction);
   cout << "Hello" << endl;
   return 0;
}