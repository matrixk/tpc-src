#ifndef _SIGNAL_H		
#define _SIGNAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>



#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TComplex.h>


void signal_1(double *re,double *im,int length){
	
	for(int i=0;i<length;i++){
		if(i>length/3&&i<length/2)re[i]=1;else re[i]=0;
		im[i]=0;
	}
}
void signal_2(double *re,double *im,int length){
	double tao1=10;
	double tao2=20;
	for(int i=0;i<length;i++){
		if(i<length/3-5*tao1)
			re[i]=0;
		else if(i>=length/3-5*tao1&&i<length/3) 
			re[i]=1-exp(-(i-length/3+5*tao1)/tao1); 
		else
			re[i]=exp(-i/tao2)/exp(-length/3/tao2);
		im[i]=0;
	}
}

#endif	
