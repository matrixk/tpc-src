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
#include <TH2D.h>
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


char gr_name[200]="diffusion";
pid=para("-gr_name");
if(pid>0){
sprintf(gr_name,"%s",argv[pid+1]);
}

char gr_inner_name[200]="gr_inner";
pid=para("-gr_inner_name");
if(pid>0){
sprintf(gr_inner_name,"%s",argv[pid+1]);
}

char gr_outer_name[200]="gr_outer";
pid=para("-gr_outer_name");
if(pid>0){
sprintf(gr_outer_name,"%s",argv[pid+1]);
}


double diffusion_range=100;
pid=para("-diffusion_range");
if(pid>0){
diffusion_range=atof(argv[pid+1]);
}

int diffusion_bin=100;
pid=para("-diffusion_bin");
if(pid>0){
diffusion_bin=atof(argv[pid+1]);
}
int diffusion_xybin=100;
pid=para("-diffusion_xybin");
if(pid>0){
diffusion_xybin=atof(argv[pid+1]);
}
int diffusion_zbin=100;
pid=para("-diffusion_zbin");
if(pid>0){
diffusion_zbin=atof(argv[pid+1]);
}
double ed_range=100;
pid=para("-ed_range");
if(pid>0){
ed_range=atof(argv[pid+1]);
}

int ed_bin=1000;
pid=para("-ed_bin");
if(pid>0){
ed_bin=atof(argv[pid+1]);
}


int plot_number=10;
pid=para("-plot_number");
if(pid>0){
plot_number=atoi(argv[pid+1]);
}


double diff_radius=1;
pid=para("-diff_radius");
if(pid>0){
diff_radius=atof(argv[pid+1]);
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
	double *rms_sum=new double[nstep];
	double *weight_sum=new double[nstep];
	double *rms_sum_outer=new double[nstep];
	double *weight_sum_outer=new double[nstep];
	double *rms_sum_inner=new double[nstep];
	double *weight_sum_inner=new double[nstep];

    for(int i=0; i<nstep; i++) {
	ye[i]=0;
	rms_sum[i]=0;
	weight_sum[i]=0;
	rms_sum_inner[i]=0;
	weight_sum_inner[i]=0;
	rms_sum_outer[i]=0;
	weight_sum_outer[i]=0;
	
	}
	// TH1D *diffusion[plot_number+1];
	// TH2D *diffusion2D[plot_number+1];
	TH1D hed("energy deposit","",ed_bin,0,ed_range);
	// char th1d_name[200];
	// char th2d_name[200];
    // for(int i=0; i<=plot_number; i++) {
	// sprintf(th1d_name,"diffusion_%f-%f",start,start+(nstep/plot_number)*step*i);
	// sprintf(th2d_name,"diffusion_2D_%f-%f",start,start+(nstep/plot_number)*step*i);
	// diffusion[i]=new TH1D(th1d_name,th1d_name,diffusion_bin,-diffusion_range,diffusion_range);
	// diffusion2D[i]=new TH2D(th2d_name,th2d_name,diffusion_bin,-diffusion_range,diffusion_range,diffusion_bin,-diffusion_range,diffusion_range);
	// }		
TH2D * diffusionXY=new TH2D("diffusionXY","diffusionXY",diffusion_xybin,-diffusion_range,diffusion_range,diffusion_xybin,-diffusion_range,diffusion_range);
TH2D * diffusionXZ=new TH2D("diffusionXZ","diffusionXZ",diffusion_xybin,-diffusion_range,diffusion_range,diffusion_zbin,start,stop);	
TH2D * diffusionYZ=new TH2D("diffusionYZ","diffusionYZ",diffusion_xybin,-diffusion_range,diffusion_range,diffusion_zbin,start,stop);	
	
	double dr=diff_radius*diff_radius;	
	cout<<"GetEntries  "<<t1->GetEntries()<<endl;
    for(int i=0; i<t1->GetEntries(); i++) {
	 // if(i>1000) break;
        t1->GetEntry(i);
    // cout<<"-------------------trackId->size()  "<<trackId->size()<<"  "<<xp->size()<<endl;		

	// for(unsigned step=0;step<trackId->size();step++){
     // cout<<(*trackId)[step]<<" ek "<<(*ek)[step]<<" ed "<<(*ed)[step]<<endl;
	// sum+=(*ed)[step];
	// }

    for(int k=0; k<trackId->size(); k++) {
	hed.Fill((*ed)[k]);	

	double r2=(pow((*xp)[k],2)+pow((*yp)[k],2));
	
    // for(int j=0; j<=plot_number; j++) {
    // int jp=j*(nstep/plot_number);
    // if((*zp)[k]>start+step*jp) continue;
	// diffusion[j]->Fill((*xp)[k],(*ed)[k]);    
	// diffusion2D[j]->Fill((*xp)[k],(*yp)[k],(*ed)[k]);    
    // }
	diffusionXY->Fill((*xp)[k],(*yp)[k],(*ed)[k]);
	diffusionXZ->Fill((*xp)[k],(*zp)[k],(*ed)[k]);
	diffusionYZ->Fill((*yp)[k],(*zp)[k],(*ed)[k]);
	
	int jm;
	if(((*zp)[k]-start)<0) jm=0;
	else jm=((*zp)[k]-start)/step+1;
	
    if(jm>=0&&jm<nstep){

	if(r2>dr){
	rms_sum_outer[jm]+=(pow((*xp)[k],2))*(*ed)[k];
	weight_sum_outer[jm]+=(*ed)[k];	
	}else{
	rms_sum_inner[jm]+=(pow((*xp)[k],2))*(*ed)[k];
	weight_sum_inner[jm]+=(*ed)[k];	
	}
	
	// rms_sum[jm]+=(pow((*xp)[k],2))*(*ed)[k];
	// weight_sum[jm]+=(*ed)[k];		
	
    }

    }
	}

    for(int i=0; i<nstep; i++) {
	
	// cout<<rms_sum[i]<<"  "<<rms_sum_inner[i]+rms_sum_outer[i]<<endl;
	rms_sum[i]=rms_sum_inner[i]+rms_sum_outer[i];
	weight_sum[i]=weight_sum_inner[i]+weight_sum_outer[i];
    }		
	
    for(int i=nstep-1; i>=0; i--) {
    for(int j=0; j<nstep; j++) {
	if(i<=j) continue;
    rms_sum[i]+=rms_sum[j];
    weight_sum[i]+=weight_sum[j];
    rms_sum_inner[i]+=rms_sum_inner[j];
    weight_sum_inner[i]+=weight_sum_inner[j];
    rms_sum_outer[i]+=rms_sum_outer[j];
    weight_sum_outer[i]+=weight_sum_outer[j];
	}}	

	
	
    for(int j=0; j<nstep; j++) {
	x[j]=start+step*j;
	xe[j]=0;
	y[j]=sqrt(rms_sum[j]/weight_sum[j]);
    ye[j]=pow(2*pow(y[j],4)/weight_sum[j],0.25);
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
	
    // TGraphErrors *gr=new TGraphErrors(nstep,x,y,xe,ye);
	// gr->SetMarkerStyle(20);
    // gr->SetName(gr_name);
	// sprintf(title_buf,"%s KeV %s energy diffusion vs %s thickness",ene,par,mat);
    // gr->SetTitle(title_buf);	
	// sprintf(title_buf,"%s thickness (mm)",mat);
    // gr->GetXaxis()->SetTitle(title_buf);
    // gr->GetYaxis()->SetTitle("energy diffusion (mm)");
    // gr->Write();

	sprintf(title_buf,"XY plane of %s KeV %s energy distribution in %s ",ene,par,mat);
    diffusionXY->SetTitle(title_buf);
	sprintf(title_buf,"%s (mm)",mat);
    diffusionXY->GetXaxis()->SetTitle(title_buf);
    diffusionXY->GetYaxis()->SetTitle(title_buf);

	sprintf(title_buf,"XZ plane of %s KeV %s energy distribution in %s ",ene,par,mat);
    diffusionXZ->SetTitle(title_buf);
	sprintf(title_buf,"%s (mm)",mat);
    diffusionXZ->GetXaxis()->SetTitle(title_buf);
    diffusionXZ->GetYaxis()->SetTitle(title_buf);	

	sprintf(title_buf,"YZ plane of %s KeV %s energy distribution in %s ",ene,par,mat);
    diffusionYZ->SetTitle(title_buf);
	sprintf(title_buf,"%s (mm)",mat);
    diffusionYZ->GetXaxis()->SetTitle(title_buf);
    diffusionYZ->GetYaxis()->SetTitle(title_buf);	
	
	sprintf(title_buf,"%s KeV %s energy deposit in %s each simulation step",ene,par,mat);
    hed.SetTitle(title_buf);
    hed.GetXaxis()->SetTitle("energy deposit (KeV)");	
	
    for(int j=0; j<nstep; j++) {
	x[j]=start+step*j;
	xe[j]=0;
	y[j]=weight_sum_outer[j]/weight_sum[j];
    ye[j]=pow(2*y[j]/t1->GetEntries(),0.5);
	}
    // TGraphErrors *gr_outer=new TGraphErrors(nstep,x,y,xe,ye);
	// gr_outer->SetMarkerStyle(20);

    // gr_outer->SetName(gr_outer_name);
	// char outer_title[200];
	// sprintf(outer_title,"energy deposition ratio out of radius %.2f VS thickness", diff_radius);	
    // gr_outer->SetTitle(outer_title);
    // gr_outer->GetXaxis()->SetTitle("CdZnTe thick mm");
    // gr_outer->GetYaxis()->SetTitle("ratio");
    // gr_outer->Write();
	
	
	

	
    for(int j=0; j<nstep; j++) {
	x[j]=start+step*j;
	xe[j]=0;
	y[j]=sqrt(rms_sum_inner[j]/weight_sum_inner[j]);
    ye[j]=pow(2*pow(y[j],4)/weight_sum_inner[j],0.25);
	}
    // TGraphErrors *gr_inner=new TGraphErrors(nstep,x,y,xe,ye);
	// gr_inner->SetMarkerStyle(20);
    // gr_inner->SetName(gr_inner_name);
	// char inner_title[200];
	// sprintf(inner_title,"energy diffusion vs CdZnTe thickness within radius %.2f", diff_radius);
    // gr_inner->SetTitle(inner_title);
    // gr_inner->GetXaxis()->SetTitle("CdZnTe thick mm");
    // gr_inner->GetYaxis()->SetTitle("diffusion mm");
    // gr_inner->Write();



	
	
	
    tof->Write();
    tf->Close();
    delete tf;
	delete []x;
	delete []y;
	delete []xe;
	delete []ye;
	delete []rms_sum;
	delete []weight_sum;
	delete []rms_sum_inner;
	delete []weight_sum_inner;
	delete []rms_sum_outer;
	delete []weight_sum_outer;
    return EXIT_SUCCESS;
}
