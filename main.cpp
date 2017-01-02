#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "structs.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include "mikadoclass.h"
#include <fstream>

const double pi=4.0*atan(1.0);

using namespace std;
using namespace Eigen;

double ftest(double *x){
	return (x[0]-1)*(x[0]-1)+(1+x[1])*(1+x[1]);
}

void grad(double *x, double *xi){
	xi[0]=2*(x[0]-1);
	xi[1]=2*(1+x[1]);
}





int main(){

	network testnetwork;
	testnetwork.import_params();
	
	cout<<"Stretch-factor = "<<testnetwork.parameters.stretchfactor<<endl;
	cout<<"Bending on = "<<testnetwork.parameters.Bendon<<endl;
	cout<<"Lenght of sticks = "<<testnetwork.parameters.Lstick<<endl;
	cout<<"The number of sticks = "<<testnetwork.parameters.NumberMikado<<endl;

	testnetwork.initiate_randomnetwork();
	testnetwork.write_sticks();
	testnetwork.write_springs();
	
	double etest=testnetwork.energy();
	cout<<etest<<endl;
	
	return 0;
}
