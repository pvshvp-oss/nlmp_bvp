#include<iostream>
#include<nlmp_bvp.hpp>

using namespace std;

VectorXd dFunction(double t, int intervalID = 0){
    return VectorXd();
}
VectorXd BCFunction(MatrixXd BC){
    return VectorXd();
}

int main(int argc, char **argv)
{
   nlmp_bvp(1, 1, dFunction, BCFunction);
   cout << "Hello" << endl;
   return 0;
}