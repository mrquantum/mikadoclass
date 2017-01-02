#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>

#include "mikadoclass.h"
#include "structs.h"

double network::energy()
{
	int one,two,wlr,wud;
	double x1,x2,y1,y2;
	double X1,X2,Y1,Y2;
	double k,rlen;
	double L;
	double E=0.0;
	
	for(int i=0;i<springlist.size();i++){
		one=springlist[i].one;
		two=springlist[i].two;
		wlr=springlist[i].wlr;
		wud=springlist[i].wud;
		k=springlist[i].k;
		rlen=springlist[i].rlen;
		
		x1=xy[one];
		x2=xy[two]+wlr;
		y1=xy[one+numbernodes];
		y2=xy[two+numbernodes]+wud;
		
		box2phys_shear(x1,y1,X1,Y1);
		box2phys_shear(x2,y2,X2,Y2);
		
		L=sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2));
		E+=k*((L-rlen)*(L-rlen));
		
	}
	return E;
	
}
