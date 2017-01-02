#ifndef MIKADOCLASS_H
#define MIKADOCLASS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/LU"

#include "importparam.h"
#include "structs.h"





class network{
	public:
		int NumberMikado;
		double LStick;
		int Bendon;

		double springk=1;
		double b_rigid;
		double sheardeformation=0;
		double gamma=0;

		
		std::vector<stick> mikadosticks;
		std::vector<node> nodes;
		std::vector<spring> springlist;
		std::vector<vector<int>> springpairs;
		std::vector<int> affine;
		
		Eigen::VectorXd XY;
		double *xy;
		
		int numbernodes;
		
		param parameters;
	
		int import_params();
		int initiate_randomnetwork();
		int set_affine(std::vector<int> &aff);
		int write_sticks();
		int write_sticks(char *filename);
		int write_springs();
		int write_springs(char *filename);
		int write_connectivity_hist(int nr2,int nr3,int nr4);
		int write_init_nodes(char *filename);
		int write_init_nodes();
		
		double Connectivity();
		int connectivity_hist();
		double energy(double *x);
		void grad(double *x, double *xi);
		
		//~ void frprmn(double p[],int n, double *fret,double (*func)(double []), void (*dfunc)(double [],double []));
		void frprmn(double p[],int n, double *fret);
		int cg_iter=0;

	private:
	
		vector<connected> Connection;
		vector<elonstick> ELONSTICK; 
		vector<int> order;
		vector<stick> mikadosticks_original; //The original copy of the sticks
		
	
		int make_sticks();
		int periodic_images();
		std::vector<stick> make_ghost_lr(std::vector<stick> &mikado);
		std::vector<stick> make_ghost_ud(std::vector<stick> &mikado);
		int make_connections();
		int sortELEMENTSperMIKADO();
		int combineElementsOnMikado(vector<connected> &connection,vector<elonstick> &elements);
		int orderElonstick();
		int makeSpringsAndNodes();
		int ordernodes();
		spring makespring(int node1,int node2,double x1,double x2, double y1, double y2,int stick,double k,double stretchf);
		int trim_network();
		int set_parameters();
		
		int box2phys_shear(double xbox,double ybox, double &X_PHYS, double &Y_PHYS);

		
		//functions for minimisation
		
		void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
		//void linmin(double p[],double xi[], int n,double *fret);
		void dlinmin(double p[], double xi[], int n, double *fret);
		double f1dim(double alpha);
		double df1dim(double x);
		double dbrent(double ax, double bx, double cx,double tol, double *xmin);
		double brent(double ax,double bx,double cx,double (*f)(double),double tol,
            double *xmin);
		double *dvector( int size );
		
		double cg_ftol=1e-36;
		double cg_lengrad;

		//function to pop elements from an eigen-style vector
		int pop_eigen(Eigen::VectorXd &X,int element);
		
		
				

};




	


#endif
