#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "mikadoclass.h"

using namespace std;

int network::box2phys_shear(double xbox,double ybox, double &X_PHYS, double &Y_PHYS){
	X_PHYS=xbox+gamma*ybox;
	Y_PHYS=ybox;
	
	return 1;
}
