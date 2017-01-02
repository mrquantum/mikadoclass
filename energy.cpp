#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>

#include "mikadoclass.h"
#include "structs.h"

double network::energy(double *xy)
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

void network::grad(double *xy,double *xi)
{
		for(int i=0;i<numbernodes*2;i++){
			xi[i]=0.0;
		}
	
	
		int one,two,wlr,wud;
		double x1,x2,y1,y2,X1,X2,Y1,Y2;
		double gradx1,gradx2,grady1,grady2;
		double L,rlen,k;
		double gradx,grady;
		
		double dX,dY,foverr;
		
		for(int i=0;i<springlist.size();i++){
			one=springlist[i].one;
			two=springlist[i].two;
			dY=xy[one+numbernodes]-xy[two+numbernodes]-springlist[i].wud; 
			dX=xy[one]+sheardeformation*dY-xy[two]-springlist[i].wlr;
			
			foverr=2*springlist[i].k*(1.0-springlist[i].rlen/sqrt(dX*dX+dY*dY));
			gradx=foverr*dX;
			grady=foverr*dY+sheardeformation*gradx;
			
			xi[one]+=gradx;
			xi[two]-=gradx;
			xi[one+numbernodes]+=grady;
			xi[two+numbernodes]-=grady;
		}	
		
		for(int i=0;i<affine.size();i++){
			xi[affine[i]]=0.0;
			xi[affine[i]+numbernodes]=0.0;
		}
		


}



//~ double network::energy(double *x){
	//~ return (x[0]-1)*(x[0]-1)+(1+x[1])*(1+x[1]);
//~ }

//~ void network::grad(double *x, double *xi){
	//~ xi[0]=2*(x[0]-1);
	//~ xi[1]=2*(1+x[1]);
//~ }
