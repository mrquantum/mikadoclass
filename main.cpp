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
#include "random.h"

const double pi=4.0*atan(1.0);

using namespace std;
using namespace Eigen;


int network::set_affine(vector<int> &aff){
	for(int i=0;i<aff.size();i++){
		affine.push_back(aff[i]);
	}

	return 1;
}




int main(int argc,char **argv){

	int SEED;
    if(argc==0) SEED=0;
    if(argc>1){
        SEED=stoi(argv[1]);
        my_random::set_seed(SEED);
    }
	
	
	
	network testnetwork;
	testnetwork.import_params();
	
	cout<<"Stretch-factor = "<<testnetwork.parameters.stretchfactor<<endl;
	cout<<"Bending on = "<<testnetwork.parameters.Bendon<<endl;
	cout<<"Lenght of sticks = "<<testnetwork.parameters.Lstick<<endl;
	cout<<"The number of sticks = "<<testnetwork.parameters.NumberMikado<<endl;

	testnetwork.initiate_randomnetwork();
	testnetwork.write_sticks();
	//testnetwork.write_springs();
	
	vector<int> affinenodes;
	affinenodes.push_back(1);
	affinenodes.push_back(4);

	testnetwork.set_affine(affinenodes);

	double xi[2*testnetwork.numbernodes];
	
	cout<<"X	Y"<<endl;
	cout<<testnetwork.xy[3]<<"\t"<<testnetwork.xy[3+testnetwork.numbernodes]<<endl;
	testnetwork.xy[3]=2;
	testnetwork.xy[3+testnetwork.numbernodes]=1;
	cout<<testnetwork.xy[3]<<"\t"<<testnetwork.xy[3+testnetwork.numbernodes]<<endl;
	
	double fret;
	testnetwork.frprmn(testnetwork.xy,2*testnetwork.numbernodes,&fret);
	testnetwork.write_init_nodes();
	
	cout<<testnetwork.xy[3]<<"\t"<<testnetwork.xy[3+testnetwork.numbernodes]<<endl;
	
	cout<<fret<<"\t"<<testnetwork.cg_iter<<endl;


	
	return 0;
}
