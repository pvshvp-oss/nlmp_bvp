#include<iostream>
#include<nlmp_bvp.hpp>

using namespace std;

VectorXd dFunction(double t, int intervalID = 0){
***REMOVED*** return VectorXd();
}
VectorXd BCFunction(MatrixXd BC){
***REMOVED*** return VectorXd();
}

int main(int argc, char **argv)
{
***REMOVED***nlmp_bvp(1, 1, dFunction, BCFunction);
***REMOVED***cout << "Hello" << endl;
***REMOVED***return 0;
}