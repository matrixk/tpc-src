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

	char output_opt[200]="RECREATE";
	pid=para("-output_opt");
	if(pid>0){
	sprintf(output_opt,"%s",argv[pid+1]);
	}
	double Amplification=200.0;
	pid=para("-Amplification");
	if(pid>0){
	Amplification=atof(argv[pid+1]);
	}
	double Threshold=100.0;
	pid=para("-Threshold");
	if(pid>0){
	Threshold=atof(argv[pid+1]);
	}
	double Energy_per_electron=0.03;
	pid=para("-Energy_per_electron");
	if(pid>0){
	Energy_per_electron=atof(argv[pid+1]);
	}
	double Electron_convert_volt=200.0;
	pid=para("-Electron_convert_volt");
	if(pid>0){
	Electron_convert_volt=atof(argv[pid+1]);
	}
	
	double Detector_Y=-10.0;
	pid=para("-Detector_Y");
	if(pid>0){
	Detector_Y=atof(argv[pid+1]);
	}
	double Detector_X_start=-5.0;
	pid=para("-Detector_X_start");
	if(pid>0){
	Detector_X_start=atof(argv[pid+1]);
	}
	double Detector_X_stop=5.0;
	pid=para("-Detector_X_stop");
	if(pid>0){
	Detector_X_stop=atof(argv[pid+1]);
	}
	
	double Detector_Z_start=0.0;
	pid=para("-Detector_Z_start");
	if(pid>0){
	Detector_Z_start=atof(argv[pid+1]);
	}
	double Detector_Z_stop=50.0;
	pid=para("-Detector_Z_stop");
	if(pid>0){
	Detector_Z_stop=atof(argv[pid+1]);
	}
	double Pixel_size=0.08;
	pid=para("-Pixel_size");
	if(pid>0){
	Pixel_size=atof(argv[pid+1]);
	}
	
	
	double Initial_cluster_size=0.08;
	pid=para("-Initial_cluster_size");
	if(pid>0){
	Initial_cluster_size=atof(argv[pid+1]);
	}
	double Diffusion=0.1;
	pid=para("-Diffusion");
	if(pid>0){
	Diffusion=atof(argv[pid+1]);
	}
	double Gauss_mean=0.0;
	pid=para("-Gauss_mean");
	if(pid>0){
	Gauss_mean=atof(argv[pid+1]);
	}
	double Gauss_sigma=20.0;
	pid=para("-Gauss_sigma");
	if(pid>0){
	Gauss_sigma=atof(argv[pid+1]);
	}
	double Theta=0.0;
	pid=para("-Theta");
	if(pid>0){
	Theta=atof(argv[pid+1]);
	}
	double Adc_volt_range=5000;//mV
	pid=para("-Adc_volt_range");
	if(pid>0){
	Adc_volt_range=atof(argv[pid+1]);
	}
	int Adc_bit=14;
	pid=para("-Adc_bit");
	if(pid>0){
	Adc_bit=atoi(argv[pid+1]);
	}
	
	double Drift_readout=1;//mm
	pid=para("-Drift_readout");
	if(pid>0){
	Drift_readout=atof(argv[pid+1]);
	}	
////////////////////////
	TFile *tf, *tof;
	TTree *t1, *t2;

	EventData *d;
	d = new EventData();

	if(help==1){cout<<endl; return 0;}
	cout<<"input  "<<input<<endl;

	tf = new TFile(input);
	tf->ls();
	t1 = (TTree*)(tf->Get("t1"));
	t1->ls();

	set_tree(t1,d);
	string out;
	out.assign(input);
	size_t size=out.size();
	out.resize(size-5);
	out.replace(size-5,12,"_signal.root");
	cout<<"output file is\t"<<out<<endl;
    // for output
	tof = new TFile(out.c_str(), output_opt);
	t2=new TTree("tree","TTree Containing Event Data From Geant4 That Electron Is Diffused");
// detector plane is parallel to Z axis and in X plane, cloud drift along Y axis
//variable for input
	double sita=Theta;	//counter-clockwise turn angle, change is in Z and y 
	double detector_Y=Detector_Y; //Y coordinate of detector plane
	double detector_X_start=Detector_X_start;
	double detector_X_stop=Detector_X_stop;
	double detector_Z_start=Detector_Z_start;
	double detector_Z_stop=Detector_Z_stop;
	double pixel_size=Pixel_size;
	double energy_per_electron=Energy_per_electron;
	double electron_convert_volt=Electron_convert_volt;
	double r0=Initial_cluster_size; //initial cluster size
	double r1=Diffusion;  //cluster size increase with distance^1/2
	double gauss_mean=Gauss_mean;
	double gauss_sigma=Gauss_sigma;

//variable for detector
	double amplification=100.0;
	int adc_bit=14;
	int digit_range=pow(2,adc_bit);
	double adc_volt_range=Adc_volt_range;
	double adc_precision=adc_volt_range*2.0/(double)digit_range;
	double threshold=Threshold;
//cout constant
	cout<<"sita\t"<<sita<<"\tamplification\t"<<amplification<<"\telectron_convert_volt\t"<<electron_convert_volt<<endl;
	cout<<"threshold\t"<<threshold<<"\tgauss_mean\t"<<gauss_mean<<"\tgauss_sigma\t"<<gauss_sigma<<endl;
//////////////////////////////////	
	double yp,zp,ed;	
	double x,y,z;	
	double he=0.25*pixel_size; //offset to avoid data in the hist edge	
	double pixel_number_z=(detector_Z_stop-detector_Z_start)/pixel_size;
	double pixel_number_x=(detector_X_stop-detector_X_start)/pixel_size;
	int event_number=t1->GetEntries();
	cout<<"The total number of event\t "<<event_number<<endl;	
	cout<<"pixel_number_z\t"<<pixel_number_z<<"\tpixel_number_x\t"<<pixel_number_x<<endl;
	int Eventid=0;
	vector<int> *Digit_z=new vector<int>;
	vector<int> *Digit_x=new vector<int>;
	vector<int> *Digit_adc=new vector<int>;
	vector<int> *Time=new vector<int>;
	t2->Branch("Eventid",&Eventid,"Eventid/I");
	t2->Branch("Digit_z","vector<int>",&Digit_z);
	t2->Branch("Digit_x","vector<int>",&Digit_x);
	t2->Branch("Digit_adc","vector<double>",&Digit_adc);
	t2->Branch("Time","vector<double>",&Time);

	TH2D * timeZX;// time information
	TH2D * signalZX;// electron cluster information
	char signal_name[100];
	char time_name[100];
	sprintf(time_name,"timeZX");
	sprintf(signal_name,"signalZX");	
	timeZX=new TH2D(time_name,time_name,pixel_number_z,detector_Z_start,detector_Z_stop,pixel_number_x,detector_X_start,detector_X_stop);
	signalZX=new TH2D(time_name,time_name,pixel_number_z,detector_Z_start,detector_Z_stop,pixel_number_x,detector_X_start,detector_X_stop);
	for(int i=0; i<event_number/100; i++) {  //for each event
		Digit_z->clear();
		Digit_x->clear();
		Digit_adc->clear();
		Time->clear();
        t1->GetEntry(i);
		Eventid=i;
	
	for(int ix=0;ix<=(detector_X_stop-detector_X_start)/pixel_size+1;ix++){
	for(int iz=0;iz<=(detector_Z_stop-detector_Z_start)/pixel_size+1;iz++){	
	timeZX->SetBinContent(iz,ix,0);	
	signalZX->SetBinContent(iz,ix,0);
	}}	

	if(i%100==0)cout<<i<<"---------------"<<d->m_nbHits<<endl;
			
	if(d->m_nbHits!=d->m_nbSteps||d->m_nbHits!=d->m_trackId->size()){
	   cout<<"---------------"<<i<<endl;
	   cout<<d->m_nbHits<<endl;
	   cout<<d->m_nbSteps<<endl;
	   cout<<d->m_trackId->size()<<endl;
			}
	//for each step,diffuse electron along y axis
	for(int j=0; j<d->m_nbHits; j++) { //for each step
	x=(*d->m_xp)[j];
	yp=(*d->m_yp)[j];
	zp=(*d->m_zp)[j];
	ed=(*d->m_energyDeposited)[j];
	z=zp*cos(sita)+yp*sin(sita);  //rotate
	y=yp*cos(sita)-zp*sin(sita);
	double size=ed/energy_per_electron;
	double distance=fabs(y-detector_Y);
	double r=r0+sqrt(distance)*r1;
	double x_start=x-r;
	double x_stop=x+r;
	double z_start=z-r;
	double z_stop=z+r;
	int x1=x_start/pixel_size;
	if(x_start<0) x1--;
	int x2=x_stop/pixel_size;
	if(x_stop>=0) x2++;
	int z1=z_start/pixel_size;
	if(z_start<0) z1--;
	int z2=z_stop/pixel_size;
	if(z_stop>=0) z2++;
	double dy_sum=0;
	for(int ix=x1;ix<=x2;ix++){
	for(int iz=z1;iz<=z2;iz++){
	double dx=ix*pixel_size+0.5*pixel_size-x;
	double dz=iz*pixel_size+0.5*pixel_size-z;
	double dy=r*r-dx*dx-dz*dz;
	if(dy<=0) continue;
	dy_sum+=sqrt(dy);
	}}
	for(int ix=x1;ix<=x2;ix++){
	for(int iz=z1;iz<=z2;iz++){
	double cx=ix*pixel_size+0.5*pixel_size;
	double cz=iz*pixel_size+0.5*pixel_size;
	double dx=cx-x;
	double dz=cz-z;
	double dy=r*r-dx*dx-dz*dz;
	int inside=0;
	if((ix*pixel_size<=x)&&(ix*pixel_size+pixel_size>x)&&(iz*pixel_size<=z)&&(iz*pixel_size+pixel_size>z))
	{inside=1;}
	if(dy<=0&&inside==0) continue;
	
	if(dy<=0&&inside==1) {
	signalZX->Fill(cz+he,cx+he,size);
		cout<<"size="<<size<<"----------------"<<endl;
	double  time=timeZX->GetBinContent((cz+he-detector_Z_start)/pixel_size+1,(cx+he-detector_X_start)/pixel_size+1);
	if(time==0){ timeZX->SetBinContent((cz+he-detector_Z_start)/pixel_size+1,(cx+he-detector_X_start)/pixel_size+1,distance-r);}
	else if(time>distance-r){
	timeZX->SetBinContent((cz+he-detector_Z_start)/pixel_size+1,(cx+he-detector_X_start)/pixel_size+1,distance-r);
	}
	break;}
	
	dy=sqrt(dy);
	signalZX->Fill(cz+he,cx+he,size*dy/dy_sum);
//			cout<<"size*dy/dy_sum="<<size*dy/dy_sum<<"----------------"<<endl;

	double  t_time=timeZX->GetBinContent((cz+he-detector_Z_start)/pixel_size+1,(cx+he-detector_X_start)/pixel_size+1);
	if(t_time==0){ timeZX->SetBinContent((cz+he-detector_Z_start)/pixel_size+1,(cx+he-detector_X_start)/pixel_size+1,distance-dy);
	}
	else if(t_time>distance-dy){
	timeZX->SetBinContent((cz+he-detector_Z_start)/pixel_size+1,(cx+he-detector_X_start)/pixel_size+1,distance-dy);
	}
	}}		
	}
	
	TRandom random;// To Create Gauss Noise
	time_t t=time(0);
	random.SetSeed(t);
	for(int ix=1;ix<=pixel_number_x;ix++)
	{
		for(int iz=1;iz<=pixel_number_z;iz++)
		{	
			double a=signalZX->GetBinContent(iz,ix)*amplification;
			double b=random.Gaus(gauss_mean,gauss_sigma);
			int adc_value=((a+b)/electron_convert_volt+adc_volt_range)/adc_precision;
			int time_value=timeZX->GetBinContent(iz,ix)/Drift_readout;
			if(adc_value>threshold)
			{
			Digit_z->push_back(iz);
			Digit_x->push_back(ix);
			Digit_adc->push_back(adc_value);
			Time->push_back(time_value);
			}
		}
	}	
	
	t2->Fill(); 		
	}
    d->Clear();

	tof->Write();
	delete Digit_z;
	delete Digit_x;
	delete Digit_adc;
	delete Time;

    return EXIT_SUCCESS;
}
