#include <TCanvas.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TPaletteAxis.h>
#include <vector>

#include <string>
#include <sstream>

#include "../lib/lib.h"

#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>
using namespace std;










int help=0;
int gargc;char** gargv;int para(char *a);int para(char *a){if(help==1)cout<<a<<" ";int i;for(i=0;i<gargc;i++){if(strcmp(a,gargv[i])==0) return i;}return -1;}




// void find_cluster(TString &filename,int &eventnumber=0,int &signal=0)
int main(int argc, char* argv[])
{

	help=0;
	gargc=argc;
	gargv=argv;
	int pid;

	pid=para("-help");if(pid>0){help=1;};


	char input[200]="";
	char input_name[200]="";
	pid=para("-input");
	if(pid>0){
	sprintf(input,"%s",argv[pid+1]);
	
	}else if(para("-input_name")>0){
	
	pid=para("-input_name");
	sprintf(input_name,"%s",argv[pid+1]);

    ifstream infile;
	infile.open(input_name);
	infile>>input;
	infile.close();
	}else{
	cout<<"No input file"<<endl;
	return 0;
	}
	
	char output[200]="rf.root";
	pid=para("-output");
	if(pid>0){
	sprintf(output,"%s",argv[pid+1]);
	}

	int photo_on=0;
	pid=para("-photo_on");
	if(pid>0){
	photo_on=atoi(argv[pid+1]);
	}	
	
	if(help==1) {
	cout<<endl;
	return 0;}
int eventnumber=1;
int signal=1;



	TFile *tf=new TFile(input);
	TTree *t1;
	t1=(TTree*)tf->Get("tree");
	digitization_tree_data *m;
    m = new digitization_tree_data();	
    m->set_digitization_tree(t1);
	
	cout<<"event number "<<t1->GetEntries()<<endl;
	

    
	
    TFile *aarf=new TFile(output,"RECREATE");
	TH2D *h=new TH2D("h","h",625,0,625,125,0,125);
	
	


	for(int i=0;i<t1->GetEntries();i++){
	// for(int i=98;i<99;i++){
	t1->GetEntry(i);		
	
	if(i==photo_on){
	for(int step=0;step<m->time1->size();step++)
	{
	if((*m->Digit_adc)[step]>0)h->SetBinContent((*m->Z_bin)[step],(*m->X_bin)[step],(*m->Digit_adc)[step]);
	}	}	

m->load_to_data();

m->finding(5);
m->set_max();
m->build_chain();	
cout<<"cluster number= "<<m->clu.size()<<endl;
cout<<"chain size= "<<m->chain.size()<<endl;
	}


	


	aarf->Write();
    delete m;

	
}
