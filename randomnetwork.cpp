#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/LU"
#include <math.h>
#include "mikadoclass.h"
#include "random.h"
#include "combineElementsOnMikado.h"
const double pi=4.0*atan(1.0);


using namespace std;
using namespace Eigen;


bool operator<(const node& first,const node& second){
	return first.number < second.number;
}

bool operator<(const elonstick& first, const elonstick& second){
	return first.sticki < second.sticki;
}

bool operator<(const stick &first, const stick &second)
{
    return first.nr < second.nr;
}

int network::import_params(){
	parameters.init("params.txt");
	return 1;
}

int network::initiate_randomnetwork()
/*This is the function that calls all functions to make a
random network*/
{ 
	make_sticks();	
	mikadosticks_original=mikadosticks;
	periodic_images();

	make_connections();
	sortELEMENTSperMIKADO();
	orderElonstick();
		
	makeSpringsAndNodes();
	
	set_parameters();
	
	
	}

int network::make_sticks()
/*Throw some sticks in a box
*/
{
	stick tempstick;
	my_random::get_gre(2);
	
	for(int i=0;i<parameters.NumberMikado;++i){
		tempstick.x=randf();
		tempstick.y=randf();
		tempstick.th=2*pi*randf();
		tempstick.nr=i;
		tempstick.wlr=0;
		tempstick.wud=0;
		tempstick.length=parameters.Lstick;
	
		mikadosticks.push_back(tempstick);
	}
	return 1;
}

int network::periodic_images()
{
	//make periodic images w. left-right wall
	//make_ghost_lr(mikadosticks);
	vector<stick> ghostlr=make_ghost_lr(mikadosticks);
	mikadosticks.insert(mikadosticks.end(),ghostlr.begin(),ghostlr.end());
	
	//make periodic images w. up-down wall
	vector<stick> ghostud=make_ghost_ud(mikadosticks);
	if(ghostud.size()!=0){
		vector<stick> ghostlr2=make_ghost_lr(ghostud);
		mikadosticks.insert(mikadosticks.end(),ghostud.begin(),ghostud.end());
		if(ghostlr2.size()!=0){
			mikadosticks.insert(mikadosticks.end(),ghostlr2.begin(),ghostlr2.end());
		}
	}
	std::sort(mikadosticks.begin(),mikadosticks.end());
	return 1;
}



vector<stick> network::make_ghost_lr(vector<stick> &mikado)
//duplictes the mikadosticks over the left-right periodic wall
{
  vector<stick> GhostLR;
  for (int i=0;i<mikado.size();i++){
    if(mikado[i].th!=pi/2 || mikado[i].th!=3*pi/2){ // Check here for the left wall 
		stick ghostRowLR;
		Matrix2d A; //The matrix to solve with
		Vector2d b,b2; // A*st=b
		Vector2d st,st2;
 
		A<<-cos(mikado[i].th),0,-sin(mikado[i].th),1; 
		b<<mikado[i].x,mikado[i].y;
		b2<<mikado[i].x-1,mikado[i].y;
		st=A.lu().solve(b);
		st2=A.lu().solve(b2);
		
		if((st(0)>0.0 && st(0)<parameters.Lstick) && (st(1)>0.0&&st(1)<1.0)){ //If mikado passes throug wall make ghost mikado 
			ghostRowLR=mikado[i];
			ghostRowLR.x=ghostRowLR.x+1;
			ghostRowLR.wlr=ghostRowLR.wlr-1;
			GhostLR.push_back(ghostRowLR); 

		}
		if((st2(0)>0.0&&st2(0)<parameters.Lstick)&&(st2(1)>0.0&&st2(1)<1.0)){ //Now do the same for the upper wall
			ghostRowLR=mikado[i];
			ghostRowLR.x=ghostRowLR.x-1;
			ghostRowLR.wlr=ghostRowLR.wlr+1;
			GhostLR.push_back(ghostRowLR);
		}
	}
  }
	return GhostLR;
}

vector<stick> network::make_ghost_ud(vector<stick> &mikado)
//suplicates the mikadosticks over the upper-down periodic wall
{
  vector<stick> GhostUD;
  for(int i=0;i<mikado.size();i++){
    if(mikado[i].th!=0||mikado[i].th!=pi){
		stick ghostRowUD;
		Matrix2d A;
		Vector2d b,b2;
		Vector2d st, st2;
      
		A<<-cos(mikado[i].th),1,-sin(mikado[i].th),0;
		b<<mikado[i].x,mikado[i].y;
		b2<<mikado[i].x,mikado[i].y-1;
		st=A.lu().solve(b);
		st2=A.lu().solve(b2);
		if((st(0)>0.0 && st(0)<parameters.Lstick)&&(st(1)>0.0 &&st(1)<1.0)){
			ghostRowUD=mikado[i];
			ghostRowUD.y=ghostRowUD.y+1;
			ghostRowUD.wud=ghostRowUD.wud-1;
			GhostUD.push_back(ghostRowUD);
		}
     
		if((st2(0)>0&&st2(0)<parameters.Lstick)&&(st2(1)>0.0&&st2(1)<1.0)){
			ghostRowUD=mikadosticks[i];
			ghostRowUD.y=ghostRowUD.y-1;
			ghostRowUD.wud=ghostRowUD.wud+1;
			GhostUD.push_back(ghostRowUD);
		}
     }
    
  }
	return GhostUD;
}

int network::make_connections()
                   {
  Matrix2d A; //The matrix to solve with
  Vector2d b,st; // A*st=b
  int nrcon=0;
  connected xtrarow, xtrarow2;
  //Loop over all sticks to find coordinates
    for(int i=0;i<mikadosticks.size()-1;i++){
        for(int j=i+1;j<mikadosticks.size();j++){
			
            if(mikadosticks[i].th!=mikadosticks[j].th || mikadosticks[i].nr!=mikadosticks[j].nr){
                A<<-cos(mikadosticks[i].th),cos(mikadosticks[j].th),-sin(mikadosticks[i].th),sin(mikadosticks[j].th);
                b<<mikadosticks[i].x-mikadosticks[j].x,mikadosticks[i].y-mikadosticks[j].y;
                st=A.lu().solve(b);
                if ((st(0)>0.0 && st(0)<parameters.Lstick)&&(st(1)>0.0 && st(1)<parameters.Lstick)){
                    xtrarow.first=mikadosticks[i].nr;	//[stick i stick j sij sji] 
                    xtrarow.second=mikadosticks[j].nr;
                    xtrarow.s1=st(0);
                    xtrarow.s2=st(1);
                    xtrarow.nrCon=nrcon;
                    xtrarow.recur=0;
                    xtrarow.type=0;
                    xtrarow.backgroundspring[0]=-1;
                    xtrarow.backgroundspring[1]=-1;
                    
                    xtrarow2.first=mikadosticks[j].nr;	//[stick j sticki sji sij]
                    xtrarow2.second=mikadosticks[i].nr;
                    xtrarow2.s1=st(1);
                    xtrarow2.s2=st(0);
                    xtrarow2.nrCon=nrcon;
                    xtrarow2.recur=0;
                    xtrarow2.type=0;
                    xtrarow2.backgroundspring[0]=-1;
                    xtrarow2.backgroundspring[1]=-1;
                    Connection.push_back(xtrarow);
                    Connection.push_back(xtrarow2);
                    nrcon++;
                }
            }
        }
    } 
    
    
    return 1;
}

int network::sortELEMENTSperMIKADO()
{
  combineElementsOnMikado(Connection,ELONSTICK);
  
// now sort extrarow on descending order per stick;
for(int j=0; j<ELONSTICK.size();j++){
    vector<double> distances=ELONSTICK[j].S;
    vector<int> numbers=ELONSTICK[j].nr;
    vector<int> type=ELONSTICK[j].type;
    vector<array<int,2>> backgroundspring=ELONSTICK[j].backgroundspringn;
// now sort extrarow on descending order per stick;
    
    int swapped=0;
      do{
        int k=0;
        swapped=0; //this is the control parameter, checks 1 if elements are swapped
            if(distances.size()>1){
                for(int i=0;i<distances.size()-1;i++){ //loop through the list thill the end-k-1 th element;
                    if(distances[i]>distances[i+1]){ //checks if neighbours are in right order, if not then swap en change swap parameter
                        double d1=distances[i];
                        double d2=distances[i+1];
                        distances[i]=d2;
                        distances[i+1]=d1;
                        int n1=numbers[i]; 
                        int n2=numbers[i+1];
                        numbers[i]=n2; 
                        numbers[i+1]=n1;
                        int t1=type[i];
                        int t2=type[i+1];
                        type[i]=t2;
                        type[i+1]=t1;
                        int backspring1 [2]={backgroundspring[i][0],backgroundspring[i][1]};
                        int backspring2 [2]={backgroundspring[i+1][0],backgroundspring[i+1][1]};
                        backgroundspring[i][0]=backspring2[0];
                        backgroundspring[i][1]=backspring2[1];
                        backgroundspring[i+1][0]=backspring1[0];
                        backgroundspring[i+1][1]=backspring1[1];
                        
                        swapped=1;
                        k++;
                    }
                }
            }
    } while(swapped==1);
        ELONSTICK[j].S=distances; //Put the new data back into the original vectors
        ELONSTICK[j].nr=numbers;
		ELONSTICK[j].type=type;
        ELONSTICK[j].backgroundspringn=backgroundspring;
    }
	std::sort(ELONSTICK.begin(),ELONSTICK.end());

	return 1;
}

int network::orderElonstick()
/*Function that ordens all the elements per stick on an ascending way of 
lenght
*/
{
    order.push_back(0);
    for(std::size_t i=0;i<ELONSTICK.size();i++){
        vector<int> numbervec=ELONSTICK[i].nr;
        for(std::size_t j=0;j<numbervec.size();j++){
            int flag=0;
            for(std::size_t k=0;k<order.size();k++){
                if(numbervec[j]==order[k]) flag=1;
            }
            if(flag==0) order.push_back(numbervec[j]);
        }
    }
    std::sort(order.begin(),order.end());
    for(std::size_t i=0;i<ELONSTICK.size();i++){
        for(std::size_t j=0;j<ELONSTICK[i].nr.size();j++){
            for(std::size_t k=0;k<order.size();k++){
                if(ELONSTICK[i].nr[j]==order[k]) ELONSTICK[i].nr[j]=k;
            }
        }
    }
    return 1;
}

int network::makeSpringsAndNodes()
{
    for(int i=0;i<ELONSTICK.size();i++){
        int sticknr=ELONSTICK[i].sticki;
        vector<int> nodesonsticki=ELONSTICK[i].nr;
        vector<double> posonsticki=ELONSTICK[i].S;
        
        stick CURRENTSTICK=mikadosticks_original[sticknr];
        int node1, node2;
        double lenspring;
        double springconstant;
        double k1=1;
        
        
        if(nodesonsticki.size()>1){
            for(int j=0;j<nodesonsticki.size()-1;j++){
                spring newspring;
                node1=nodesonsticki[j];
                node2=nodesonsticki[j+1];
                lenspring=posonsticki[j+1]-posonsticki[j];
                springconstant=lenspring*k1/CURRENTSTICK.length;
                double x1, x2, y1, y2;
                x1=CURRENTSTICK.x+posonsticki[j]*cos(CURRENTSTICK.th); //calculate the position of the node
                x2=CURRENTSTICK.x+posonsticki[j+1]*cos(CURRENTSTICK.th);//and the position of the adjacent one
                y1=CURRENTSTICK.y+posonsticki[j]*sin(CURRENTSTICK.th);
                y2=CURRENTSTICK.y+posonsticki[j+1]*sin(CURRENTSTICK.th);
                newspring=makespring(node1,node2,x1,x2,y1,y2,sticknr,springconstant,parameters.stretchfactor);
                node nodetemp1, nodetemp2;
                nodetemp1.number=newspring.one;
                nodetemp1.x=x1-floor(x1);
                nodetemp1.y=y1-floor(y1);
                
                nodetemp2.number=newspring.two;
                nodetemp2.x=x2-floor(x2);
                nodetemp2.y=y2-floor(y2);
                nodes.push_back(nodetemp1);
                nodes.push_back(nodetemp2);
                springlist.push_back(newspring);              
            }
        }
    }
    //sort the nodes and if the background absent remove double nodes.
    //in the other case this happens in makeanddeletebondsonbackground
    std::sort(nodes.begin(),nodes.end());
    ordernodes();

 
    int mi=0;
    int ma=0;
    for(int i=0;i<springlist.size();i++){
        springlist[i].one>ma ? ma=springlist[i].one : ma;
        springlist[i].two>ma ? ma=springlist[i].two : ma;

    }
    
    //Put all the nodes in an vector XY=<x0 x1 x2 ... x_N;y0,y0 y1 y2 ... y_N>
	VectorXd X(nodes.size());
	VectorXd Y(nodes.size());
	for(int i=0;i<X.size();i++){
		X[i]=nodes[i].x;
	}
	for(int i=0;i<Y.size();i++){
		Y[i]=nodes[i].y;
	}
	XY.resize(2*X.size());
	XY<<X,Y;
	

	return 1;
	
	
}


spring network::makespring(int node1,int node2,double x1,double x2, double y1, double y2,int stick,double k,double stretchf)
{
    spring newspring;
    newspring.one=node1;
    newspring.two=node2;
    newspring.rlen=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))/stretchf;
    newspring.k=k;
    newspring.sticki=stick;
    
    if((x1<1 && x1>0)&&(x2>0&&x2<1)){ //Check if crossed wlr or wud wall.
        newspring.wlr=0;
    }
    else if((x1>0&&x1<1)&&x2>1){
        newspring.wlr=1;       
    }
    else if((x1>0&&x1<1)&&x2<0){
        newspring.wlr=-1;
    }
    else if(x1<0&&x2<0){
        newspring.wlr=0;
    }
    else if(x2>0 && x2<1.0 && x1<0.0){
        newspring.wlr=1;
    }
    else if(x2>0 && x2<1.0 && x1>1.0){
          newspring.wlr=-1;
    }
    else if(x1>1&&x2>1){
        newspring.wlr=0;
    }
    if((y1<1 && y1>0)&&(y2>0&&y2<1)){
        newspring.wud=0;
    }
    else if((y1>0&&y1<1)&&y2>1){
        newspring.wud=1;
    }
    else if((y1>0&&y1<1)&&y2<0){
        newspring.wud=-1;
    }
    else if(y1<0&&y2<0){
        newspring.wud=0;
    }
    else if(y1>1&&y2>1){
        newspring.wud=0;
    }
    else if(y2>0 && y2<1.0 && y1<0.0){
        newspring.wud=1;
    }
    else if(y2>0 && y2<1.0 && y1>1.0){
        newspring.wud=-1;
    }

    return newspring;
}


int network::ordernodes(){
    
    std::sort(nodes.begin(),nodes.end());

    vector<int> remove(0);
    
    for(int i=0;        i<nodes.size()-1;       i++){
        if(nodes[i].number==nodes[i+1].number){
            remove.push_back(i+1);
        }
    }
    
    //now erase the elements from the remove vector
    int del;
    while(remove.size()>0){
        del=remove[remove.size()-1];
        remove.pop_back();
        nodes.erase(nodes.begin()+del);
    }
    
    // next step: renumber and do the same to the springs
    int num;
    for(int i=0;        i<nodes.size(); i++){
        if(nodes[i].number!=i){
           num=nodes[i].number; 
           //correct in the springlist
           for(int j=0; j<springlist.size();    j++){
               if(springlist[j].one==num){
                   springlist[j].one=i;
               } 
               if(springlist[j].two==num){
                   springlist[j].two=i;
               }
            }
           nodes[i].number=i;
        }
    }
    
 
    return 1;

}

int network::set_parameters()
{
	xy=XY.data();
	numbernodes=XY.size()/2;
	
	return 1;
	
}



int network::connectivity_hist()
{
	int nr_2=0;
	int nr_3=0;
	int nr_4=0;
	
	int maxele=0;
	for(int i=0;i<springlist.size();i++){
		if(springlist[i].one>maxele){
			maxele=springlist[i].one;
		}
		if(springlist[i].two>maxele){
			maxele=springlist[i].two;
	}
}
	
	int count=0;
	for(int i=0;i<=maxele;i++){
		count=0;
		for(int j=0;j<springlist.size();j++){
			if(springlist[j].one==i){
				count++;
			}
			if(springlist[j].two==i){
				count++;
			}
		}
		if(count==2) nr_2++;
		if(count==3) nr_3++;
		if(count==4) nr_4++;
	}

	write_connectivity_hist(nr_2,nr_3,nr_4);

	
	return 1;
	
}







int network::write_sticks(){
	ofstream mikadofile("mikado.txt");
    for(int i=0;i<mikadosticks.size();i++){
        mikadofile<<mikadosticks[i].nr<<"\t"
        <<mikadosticks[i].x<<"\t"
        <<mikadosticks[i].y<<"\t"
        <<mikadosticks[i].th<<"\t"
        <<mikadosticks[i].wlr<<"\t"
        <<mikadosticks[i].wud<<endl;
    } 
    mikadofile.close();
    return 1;
};

int network::write_sticks(char *filename){
	ofstream mikadofile(filename);
    for(int i=0;i<mikadosticks.size();i++){
        mikadofile<<mikadosticks[i].nr<<"\t"
        <<mikadosticks[i].x<<"\t"
        <<mikadosticks[i].y<<"\t"
        <<mikadosticks[i].th<<"\t"
        <<mikadosticks[i].wlr<<"\t"
        <<mikadosticks[i].wud<<endl;
    } 
    mikadofile.close();
    return 1;
};







int network::write_springs()
{
	ofstream springfile("springs.txt");
    for(int i=0;i<springlist.size();i++){
        springfile<<springlist[i].one<<"\t"
                  <<springlist[i].two<<"\t"
                <<springlist[i].wlr<<"\t"
                <<springlist[i].wud<<"\t"
                <<springlist[i].rlen<<"\t"
                <<springlist[i].k<<"\t"
                <<springlist[i].sticki<<endl;       
    }
    springfile.close();
    return 1;
};



int network::write_springs(char *filename)
{
	ofstream springfile(filename);
    for(int i=0;i<springlist.size();i++){
        springfile<<springlist[i].one<<"\t"
                  <<springlist[i].two<<"\t"
                <<springlist[i].wlr<<"\t"
                <<springlist[i].wud<<"\t"
                <<springlist[i].rlen<<"\t"
                <<springlist[i].k<<"\t"
                <<springlist[i].sticki<<endl;       
    }
    springfile.close();
    return 1;
};


int network::write_connectivity_hist(int nr2,int nr3,int nr4)
{
	ofstream c_hist("chist.txt");
	c_hist<<nr2<<"\t"<<nr3<<"\t"<<nr4<<"\t"<<parameters.NumberMikado<<"\t"<<parameters.Lstick<<endl;
	c_hist.close();
	return 1;
	
}





