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
#include "../lib/lib.h"

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

char input[200]="";
pid=para("-input");
if(pid>0){
sprintf(input,"%s",argv[pid+1]);
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

char ratio_name[200]="ratio";
pid=para("-ratio_name");
if(pid>0){
sprintf(ratio_name,"%s",argv[pid+1]);
}

double energy=100;
pid=para("-energy");
if(pid>0){
energy=atof(argv[pid+1]);
}

int energy_bin=100;
pid=para("-energy_bin");
if(pid>0){
energy_bin=atof(argv[pid+1]);
}

int plot_number=10;
pid=para("-plot_number");
if(pid>0){
plot_number=atof(argv[pid+1]);
}

////////////////////////
TFile *tf, *tof;
TTree *t1, *t2;



if(help==1){cout<<endl; return 0;}
	cout<<"input  "<<input<<endl;

    tf = new TFile(input);
    tf->ls();
    t1 = (TTree*)(tf->Get("t1"));
    t1->ls();

    t1->SetBranchAddress("nbsteps", &nbSteps);
    t1->SetBranchAddress("trackid", &trackId);
    t1->SetBranchAddress("edproc", &edProc);
    t1->SetBranchAddress("etot", &etot);
    t1->SetBranchAddress("ed", &ed);
    t1->SetBranchAddress("ek", &ek);
    t1->SetBranchAddress("xp", &xp);
    t1->SetBranchAddress("yp", &yp);
    t1->SetBranchAddress("zp", &zp);

    // for output
    tof = new TFile(output, output_opt);
	
int nstep=(stop-start)/step;
	double *x=new double[nstep];
	double *y=new double[nstep];
	double *xe=new double[nstep];
	double *ye=new double[nstep];
	double *de=new double[nstep];
	double *de2=new double[nstep];
	double *dedx=new double[nstep];
	double *dedx2=new double[nstep];
    for(int i=0; i<nstep; i++) {
	de[i]=0;
	de2[i]=0;
	dedx[i]=0;
	dedx2[i]=0;
	}
	TH1D *energy_deposit[plot_number+1];
	char th1d_name[200];
    for(int i=0; i<=plot_number; i++) {
	sprintf(th1d_name,"energy_deposit_%f-%f",start,start+(nstep/plot_number)*step*i);
	energy_deposit[i]=new TH1D(th1d_name,th1d_name,energy_bin,0,energy);
	}	
	// double sum=0;	
	cout<<"GetEntries  "<<t1->GetEntries()<<endl;
    for(int i=0; i<t1->GetEntries(); i++) {
	// if(i>10) break;
        t1->GetEntry(i);
    // cout<<"-------------------trackId->size()  "<<trackId->size()<<"  "<<xp->size()<<endl;		

	// for(unsigned step=0;step<trackId->size();step++){
     // cout<<(*trackId)[step]<<" ek "<<(*ek)[step]<<" ed "<<(*ed)[step]<<endl;
	// sum+=(*ed)[step];
	// }


	
    for(int j=0; j<nstep; j++) {
	double sum=0;
	double sumdx=0;
    for(int k=0; k<trackId->size(); k++) {
	
	if((*zp)[k]>start+step*(j-1)&&(*zp)[k]<=start+step*j){
	sumdx+=(*ed)[k];
	}
	
	if((*zp)[k]>start+step*j)break;
	sum+=(*ed)[k];  //add the deposit energy
	}

	if(j%(nstep/plot_number)==0&&j/(nstep/plot_number)<=plot_number){
	energy_deposit[j/(nstep/plot_number)]->Fill(sum);
	}
	
	de[j]+=sum;
	de2[j]+=sum*sum;
	
	dedx[j]+=sumdx;
	dedx2[j]+=sumdx*sumdx;
	
    }
	
	}
// energy filled into hist 
double n=t1->GetEntries();

    for(int j=0; j<nstep; j++) {
	x[j]=start+step*j;
	xe[j]=0;
	y[j]=de[j]/n;
	ye[j]=sqrt(de2[j]/n-pow(y[j],2));
	//ye[j]=y[j]/sqrt(n);
	}
	
char mat[200];
lib_matter(input,mat);
// cout<<"matter  is "<<mat<<endl;
char par[200];
lib_particle(input,par);
// cout<<"particle  is "<<par<<endl;
char ene[200];
lib_energy(input,ene);
// cout<<"energy  is "<<ene<<endl;

char title_buf[2000];

    TGraphErrors *gr=new TGraphErrors(nstep,x,y,xe,ye);
	gr->SetMarkerStyle(20);
    // td->Draw("ap");
    gr->SetName(gr_name);
	sprintf(title_buf,"%s KeV %s energy loss vs %s thickness",ene,par,mat);
    gr->SetTitle(title_buf);
	sprintf(title_buf,"%s thickness (mm)",mat);
    gr->GetXaxis()->SetTitle(title_buf);
    gr->GetYaxis()->SetTitle("dE (KeV)");
    gr->Write();
	
    for(int j=0; j<nstep; j++) {
	y[j]=y[j]/energy;
	ye[j]=ye[j]/energy;
	}
	

    TGraphErrors *ratio=new TGraphErrors(nstep,x,y,xe,ye);
	ratio->SetMarkerStyle(20);
    // td->Draw("ap");
    ratio->SetName(ratio_name);
	sprintf(title_buf,"%s KeV %s energy loss vs %s thickness",ene,par,mat);
    ratio->SetTitle(title_buf);
	sprintf(title_buf,"%s thickness (mm)",mat);
    ratio->GetXaxis()->SetTitle(title_buf);
    ratio->GetYaxis()->SetTitle("dE/E");
    ratio->Write();	


    for(int j=0; j<nstep; j++) {
	x[j]=start+step*j;
	xe[j]=0;
	y[j]=dedx[j]/n/step;
	ye[j]=sqrt(dedx2[j]/n-pow(y[j]*step,2))/step;
	//ye[j]=y[j]/sqrt(n);
	}
	

    TGraphErrors *grdx=new TGraphErrors(nstep,x,y,xe,ye);
	grdx->SetMarkerStyle(20);
    // td->Draw("ap");
	grdx->SetName("dedx");
	sprintf(title_buf,"%s KeV %s dedx vs %s thickness",ene,par,mat);
    grdx->SetTitle(title_buf);
	sprintf(title_buf,"%s thickness (mm)",mat);
    grdx->GetXaxis()->SetTitle(title_buf);
    grdx->GetYaxis()->SetTitle("dedx (KeV/mm)");	
    grdx->Write();
	
	
	
    tof->Write();
    tf->Close();
    delete tf;
	delete []x;
	delete []y;
	delete []xe;
	delete []ye;
	delete []de;
	delete []de2;
    return EXIT_SUCCESS;
}
