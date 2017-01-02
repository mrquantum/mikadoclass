#ifndef IMPORTPARAM_H
#define IMPORTPARAM_H

#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include <stdio.h>
#include <string.h>

using namespace std;

class param
{
public: 
	int NumberMikado,Bendon;
    double Lstick,stretchfactor; //Number of mikado, max number of conjugate iterations
    param()
    {
     Bendon=0;
     stretchfactor=1.0;
     Lstick=.0;
     NumberMikado=0;    
    }  
    
    ~param(){         
    }
    
    int validate() //Checks the parameters 
    {
       return 0;
    }
    
    int parse_param_assignment(char *str)
    {
        const int buffsize=101;
        char lvalue[buffsize], rvalue[buffsize];
        int rv=sscanf(str,"%100[^=]=%s",lvalue,rvalue);
		if(rv!=2 && strcmp(lvalue,"end")!=0){
            cout<<"error in param assignment"<<endl;

        }
        // rv=sscanf(lvalue,"%s",lvalue);
        // cout<<lvalue<<endl;
        
        if(strcmp(lvalue,"LStick")==0){
            Lstick=atof(rvalue);
        }
		else if(strcmp(lvalue,"NumberMikado")==0){
		NumberMikado=atoi(rvalue);
		}
		else if(strcmp(lvalue,"Stretchfactor")==0){
		stretchfactor=atof(rvalue);
		}
		else if(strcmp(lvalue,"Bend_on")==0){
			Bendon=atoi(rvalue);
		}
		else if(strcmp(lvalue,"end")==0){
		return 1;
		}
        else{
            cout<<"Did not recognise param"<<endl;
            return -1;            
            
        }
        return 1;   
        
    }   

	int get_next_script_tag(FILE *in, char *buf);
	int init(const char *FILENAME);
    
};












#endif // IMPORTPARAM_H
