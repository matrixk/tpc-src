/*This Program Calculate The Steps Of Gamma Interaction which is Compton
Scatterig Or Photelectric Effect, And Will Create A Root FIle Named _t2.root.
This File Can Be Run By:./GammaMultiScatter InputFile Energy
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
#include <iostream>
#include <string>
#include <vector>

using namespace std;
TFile *tf, *tof;
TTree *t1, *t2;

// for input
int nbSteps;
vector<int> *trackId;
vector<string> *edProc;
vector<double> *xp, *yp, *zp, *ed, *ek;
double etot;
// for output
int t2NbSteps;
vector<int> *gammaStepIndex;
double Ee;//Electron energy in the first step created by compton scattering
double Ee_cal;//Electron energy calculated by incident gamma energy direction and outgoing gamma direction and energy
double Ee_sim;//simulated first electron energy
double Ee_cal_theta;//calculated electron outing theta with xaxis
double Ee_sim_theta;//simulated electron outing theta with xaxis
double Eg;//total energy deposited in LXenon
double theta;//incident gamma direction calculated
double E_incident=0;

int help=0;
int gargc;char** gargv;int para(char *a);int para(char *a){if(help==1)cout<<a<<" ";int i;for(i=0;i<gargc;i++){if(strcmp(a,gargv[i])==0) return i;}return -1;}

int main(int argc, char **argv)
{

help=0;
gargc=argc;
gargv=argv;
int pid;
pid=para("-help");if(pid>0){help=1;};

char input_de[200]="";
pid=para("-input_de");
if(pid>0){
sprintf(input_de,"%s",argv[pid+1]);
}

char input_sigma[200]="";
pid=para("-input_sigma");
if(pid>0){
sprintf(input_sigma,"%s",argv[pid+1]);
}


double start=0;
pid=para("-start");
if(pid>0){
start=atof(argv[pid+1]);
}

double stop=10;
pid=para("-stop");
if(pid>0){
stop=atof(argv[pid+1]);
}

double step=1;
pid=para("-step");
if(pid>0){
step=atof(argv[pid+1]);
}


char output[200]="";
pid=para("-output");
if(pid>0){
sprintf(output,"%s",argv[pid+1]);
}

char output_opt[200]="RECREATE";
pid=para("-output_opt");
if(pid>0){
sprintf(output_opt,"%s",argv[pid+1]);
}


char gr_name[200]="gr";
pid=para("-gr_name");
if(pid>0){
sprintf(gr_name,"%s",argv[pid+1]);
}

char sel_dif[200]="gr_inner";
pid=para("-sel_dif");
if(pid>0){
sprintf(sel_dif,"%s",argv[pid+1]);
}


////////////////////////
TFile *tf, *tof;
TTree *t1, *t2;
//For reading TGraph
TFile *de_gr,*sigma_gr;
TGraphErrors *gr1;
TGraphErrors *gr2;



if(help==1){cout<<endl; return 0;}
	cout<<"input_de  "<<input_de<<endl;
	cout<<"input_sigma  "<<input_sigma<<endl;

    de_gr = new TFile(input_de);
    // tf->ls();
    gr1 = (TGraphErrors*)(de_gr->Get("ratio"));
	
	sigma_gr = new TFile(input_sigma);
    // tf->ls();
    gr2 = (TGraphErrors*)(sigma_gr->Get(sel_dif));
	
	
    // t1->ls();

    // t1->SetBranchAddress("nbsteps", &nbSteps);
    // t1->SetBranchAddress("trackid", &trackId);
    // t1->SetBranchAddress("edproc", &edProc);
    // t1->SetBranchAddress("etot", &etot);
    // t1->SetBranchAddress("ed", &ed);
    // t1->SetBranchAddress("ek", &ek);
    // t1->SetBranchAddress("xp", &xp);
    // t1->SetBranchAddress("yp", &yp);
    // t1->SetBranchAddress("zp", &zp);

    // for output
    tof = new TFile(output, output_opt);
	
int nstep=(stop-start)/step;
	double *x_de=new double[nstep];
	double *x_sigma=new double[nstep];

	double *y_de=new double[nstep];
	double *y_sigma=new double[nstep];

	double *error_de=new double[nstep];
	double *error_sigma=new double[nstep];
	
	double *de_sigma=new double[nstep];
	double *xe=new double[nstep];
	double *ye=new double[nstep];

	double *sum=new double[nstep];
	
    for(int i=0; i<nstep; i++) {
	x_de[i]=0;
	y_de[i]=0;
	x_sigma[i]=0;
	y_sigma[i]=0;
	de_sigma[i]=0;
	xe[i]=0;
	ye[i]=0;

	}
	
	for(int i=0; i<nstep; i++) {
	gr1->GetPoint(i,x_de[i],y_de[i]);
	error_de[i]=gr1->GetErrorY(i);
	gr2->GetPoint(i,x_sigma[i],y_sigma[i]);
	error_sigma[i]=gr2->GetErrorY(i);
	de_sigma[i]=y_sigma[i]/sqrt(y_de[i]);
	
	xe[i]=0;
	ye[i]=sqrt(pow(error_de[i],2)*pow(y_sigma[i],2)/(pow(y_de[i],3)*4)+pow(error_sigma[i],2)/y_de[i]);
	// cout<<de_sigma[i]<<endl;
	}
	cout<<"done"<<endl;

	


    TGraphErrors *gr=new TGraphErrors(nstep,x_de,de_sigma,xe,ye);
	gr->SetMarkerStyle(20);
    // td->Draw("ap");
    gr->SetName(gr_name);
    gr->SetTitle("resolution per photon (mm) vs CdZnTe thickness");
    gr->GetXaxis()->SetTitle("CdZnTe thick mm");
    gr->GetYaxis()->SetTitle("resolution per photon mm");

	


    gr->Write();
    tof->Write();
	
    delete tf;
	delete []x_de;
	delete []y_de;
	delete []error_de;
	delete []x_sigma;
	delete []y_sigma;
	delete []error_sigma;
	delete []de_sigma;
	delete []xe;
	delete []ye;
	delete []sum;

    return EXIT_SUCCESS;
}
