/*This Program will create a root file. If the input root file is abc.root, then the output file is abc_signal.root
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <string>
#include <vector>
#include "../lib/lib.h"
#include <TRandom.h>
#include <time.h>
using std::string;
using std::vector;

using namespace std;

int help=0;
int gargc;char** gargv;int para(char *a);int para(char *a){if(help==1)cout<<a<<" ";int i;for(i=0;i<gargc;i++){if(strcmp(a,gargv[i])==0) return i;}return -1;}

int main(int argc, char **argv)
{

	
	help=0;
	gargc=argc;
	gargv=argv;
	int pid;
	pid=para("-help");if(pid>0){help=1;};

	char input[200]="";
	pid=para("-input");
	if(pid>0){
	sprintf(input,"%s",argv[pid+1]);
	}


	
	cout<<"input  "<<input<<endl;
	
	
//////////////////////
	TFile *tf, *tof;
	TTree *t1, *t2;

	EventData *d;
	d = new EventData();

	if(help==1){cout<<endl; return 0;}


	tf = new TFile(input);
	tf->ls();
	t1 = (TTree*)(tf->Get("t1"));
	t1->ls();

	set_tree(t1,d);
	int event_number=t1->GetEntries();	
		cout<<"event_number  "<<event_number<<endl;
	for(int i=0; i<event_number; i++) {  //for each event

        t1->GetEntry(i);	
	   cout<<d->m_nbHits<<endl;

	for(int j=0; j<d->m_nbHits; j++) {	
	

         cout<<(*d->m_xp)[j]<<endl;	
	}
	
	}
    return EXIT_SUCCESS;
}
