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
#include <fstream>
#include <string>
#include <vector>
#include "../lib/lib.h"
#include <TRandom.h>
#include <time.h>
using std::string;
using std::vector;

using namespace std;



int main(int argc, char **argv)
{
	
 paras par;
 if(argc==1){
 par.run("para.txt"); 
 }
 else if(argv[1][0]!='-'){
 par.run(argv[1]);
 }
 else{
	 par.gargc=argc;
	 par.gargv=argv;
 
 }

 
 int pid=-1; 
	pid=par.para("-help");
	if(pid>=0){
		par.help=1;
	}
 
	char input[200]="";
	pid=par.para("-input");
	if(pid>0){
	sprintf(input,"%s",par.v(1));
	}

	pid=par.para("-inputFile");
			
	if(pid>0){
		ifstream in;
		in.open(par.v(1));
		in>>input;
		in.close();
	
	}	

 	
	char output_opt[200]="RECREATE";
	pid=par.para("-output_opt");
	if(pid>0){
	sprintf(output_opt,"%s",par.v(1));
	}
	
	char output_name[200]="output_name.txt";
	pid=par.para("-output_name");
	if(pid>0){
	sprintf(output_name,"%s",par.v(1));
	}	
	
	pid=par.para("-outputFile");
	if(pid>0){
		ifstream in;
		in.open(par.v(1));
		in>>output_name;
		in.close();
	}	
	
	
	double amplification=1.0;
	pid=par.para("-amplification");
	if(pid>0){
	amplification=atof(par.v(1));
	}
	double threshold=100.0;  //in ADC count
	pid=par.para("-threshold");
	if(pid>0){
	threshold=atof(par.v(1));
	}
	double energyPerElectron=0.03;
	pid=par.para("-energyPerElectron");
	if(pid>0){
	energyPerElectron=atof(par.v(1));
	}
	double electronPerAdcCnt=1.0;   //electron number per ADC count
	pid=par.para("-Electron_per_adc_cnt");
	if(pid>0){
	electronPerAdcCnt=atof(par.v(1));
	}
	
	double detectorY=-5.0;
	pid=par.para("-detectorY");
	if(pid>0){
	detectorY=atof(par.v(1));
	}
	double detectorXStart=-5.0;
	pid=par.para("-detectorXStart");
	if(pid>0){
	detectorXStart=atof(par.v(1));
	}
	double detectorXStop=5.0;
	pid=par.para("-detectorXStop");
	if(pid>0){
	detectorXStop=atof(par.v(1));
	}

	double detectorYStart=-5.0;
	pid=par.para("-detectorYStart");
	if(pid>0){
	detectorYStart=atof(par.v(1));
	}
	double detectorYStop=5.0;
	pid=par.para("-detectorYStop");
	if(pid>0){
	detectorYStop=atof(par.v(1));
	}
	
	double detectorZStart=0.0;
	pid=par.para("-detectorZStart");
	if(pid>0){
	detectorZStart=atof(par.v(1));
	}
	double detectorZStop=50.0;
	pid=par.para("-detectorZStop");
	if(pid>0){
	detectorZStop=atof(par.v(1));
	}
	double pixelSize=0.08;
	pid=par.para("-pixelSize");
	if(pid>0){
	pixelSize=atof(par.v(1));
	}
	
	
	double initialClusterSize=0.08;
	pid=par.para("-initialClusterSize");
	if(pid>0){
	initialClusterSize=atof(par.v(1));
	}
	double diffusion=0.1;
	pid=par.para("-diffusion");
	if(pid>0){
	diffusion=atof(par.v(1));
	}
	double gaussMean=0.0;
	pid=par.para("-gaussMean");
	if(pid>0){
	gaussMean=atof(par.v(1));
	}
	double gaussSigma=20.0;
	pid=par.para("-gaussSigma");
	if(pid>0){
	gaussSigma=atof(par.v(1));
	}
	double sita=0.0;  //counter-clockwise turn angle, change is in Z and y 
	pid=par.para("-sita");
	if(pid>0){
	sita=atof(par.v(1));
	}
	double adcVoltRange=3000;//mV
	pid=par.para("-adcVoltRange");
	if(pid>0){
	adcVoltRange=atof(par.v(1));
	}
	int adcBits=14;
	pid=par.para("-adcBits");
	if(pid>0){
	adcBits=atoi(par.v(1));
	}
	

	double driftSpeed=1;
	pid=par.para("-driftSpeed");
	if(pid>0){
	driftSpeed=atof(par.v(1));
	}	


	double noiseOnPixel=10;
	pid=par.para("-noiseOnPixel");
	if(pid>0){
	noiseOnPixel=atof(par.v(1));
	}	

	double adcSamplingRate=1;
	pid=par.para("-adcSamplingRate");
	if(pid>0){
	adcSamplingRate=atof(par.v(1));
	}	
	


////////////////////////
	TFile *tf, *tof;
	TTree *t1, *t2;

	EventData *d;
	d = new EventData();

	if(par.help==1){cout<<endl; return 0;}
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
 
	ofstream output_name_file;
	output_name_file.open(output_name);
	output_name_file<<out<<endl;
	output_name_file.close();
	
	cout<<"output file is\t"<<out<<endl;
    // for output
	tof = new TFile(out.c_str(), output_opt);
	t2=new TTree("tree","TTree Containing Event Data From Geant4 That Electron Is Diffused");
// detector plane is parallel to Z axis and in X plane, cloud drift along Y axis
//variable for input

	double r0=initialClusterSize; //initial cluster electron_num
	double r1=diffusion;  //cluster electron_num increase with distance^1/2


//cout constant
	cout<<"sita\t"<<sita<<"\tamplification\t"<<amplification<<"\telectron_per_adc_cnt\t"<<electronPerAdcCnt<<endl;
	cout<<"threshold\t"<<threshold<<"\tgauss_mean\t"<<gaussMean<<"\tgauss_sigma\t"<<gaussSigma<<endl;
//////////////////////////////////	
	double yp,zp,ed;	
	double x,y,z;	
	double he=0.25*pixelSize; //offset to avoid data in the hist edge	
	double pixelNumberZ=(detectorZStop-detectorZStart)/pixelSize;
	double pixelNumberX=(detectorXStop-detectorXStart)/pixelSize;
	double pixelNumberY=(detectorYStop-detectorYStart)*adcSamplingRate/driftSpeed;
	int event_number=t1->GetEntries();
	cout<<"The total number of event\t "<<event_number<<endl;	
	cout<<"pixelNumberZ\t"<<pixelNumberZ<<"\tpixel_number_x\t"<<pixelNumberX<<endl;
	int Eventid=0;
	vector<int> *hitIsNoise=new vector<int>;
	vector<int> *Digit_z=new vector<int>;
	vector<int> *Digit_x=new vector<int>;
	vector<int> *Digit_adc=new vector<int>;
	vector<int> *Time=new vector<int>;
	
	vector<int> *hitIsNoiseT=new vector<int>;
	vector<int> *Digit_zT=new vector<int>;
	vector<int> *Digit_xT=new vector<int>;
	vector<int> *Digit_adcT=new vector<int>;
	vector<int> *TimeT=new vector<int>;
	
	
	t2->Branch("Eventid",&Eventid,"Eventid/I");
	t2->Branch("hitIsNoise","vector<int>",&hitIsNoise);
	t2->Branch("Digit_z","vector<int>",&Digit_z);
	t2->Branch("Digit_x","vector<int>",&Digit_x);
	t2->Branch("Digit_adc","vector<int>",&Digit_adc);
	t2->Branch("Time","vector<int>",&Time);

	TH2D * timeZX;// time information
	TH2D * signalZX;// electron cluster information
	TH1D * nosiePassThreshold;// electron cluster information
	char signal_name[100];
	char time_name[100];
	char noisePT[100];
	sprintf(time_name,"timeZX");
	sprintf(signal_name,"signalZX");	
	sprintf(noisePT,"nosiePassThreshold");	
	timeZX=new TH2D(time_name,time_name,pixelNumberZ,detectorZStart,detectorZStop,pixelNumberX,detectorXStart,detectorXStop);
	signalZX=new TH2D(signal_name,signal_name,pixelNumberZ,detectorZStart,detectorZStop,pixelNumberX,detectorXStart,detectorXStop);
	nosiePassThreshold=new TH1D(noisePT,noisePT,10000,threshold,threshold+10*noiseOnPixel);
	
	for(int i=0;i<10000;i++){
		nosiePassThreshold->Fill((double(i)+0.0001)*(noiseOnPixel/1000.0)+threshold,exp(-pow((double(i)*(noiseOnPixel/1000.0)+threshold),2)/pow(noiseOnPixel,2)));
	}
	
	
	
	
	
	TRandom random;// To Create Gauss Noise
	time_t t=time(0);
	random.SetSeed(t);	
	
	
int fakeHitNum=(1-erf(threshold/noiseOnPixel))*pixelNumberZ*pixelNumberX*pixelNumberY;
	cout<<"threshold "<<threshold<<endl;
	cout<<"noiseOnPixel "<<noiseOnPixel<<endl;
	cout<<"erf(threshold/noiseOnPixel) "<<erf(threshold/noiseOnPixel)<<endl;
	cout<<"fakeHitNum "<<fakeHitNum<<endl;
	
	for(int i=0; i<event_number; i++) {  //for each event
	
		hitIsNoise->clear();
		Digit_z->clear();
		Digit_x->clear();
		Digit_adc->clear();
		Time->clear();

		hitIsNoiseT->clear();
		Digit_zT->clear();
		Digit_xT->clear();
		Digit_adcT->clear();
		TimeT->clear();
		
		
        t1->GetEntry(i);
		Eventid=i;
	
	for(int ix=0;ix<=(detectorXStop-detectorXStart)/pixelSize+1;ix++){
	for(int iz=0;iz<=(detectorZStop-detectorZStart)/pixelSize+1;iz++){	
	timeZX->SetBinContent(iz,ix,0);	
	signalZX->SetBinContent(iz,ix,0);
	}}	

	//if(i%10==0)
		cout<<i<<"---------------"<<d->m_nbHits<<endl;
			
	if(d->m_nbHits!=d->m_nbSteps||d->m_nbHits!=d->m_trackId->size()){
	   cout<<"---------------"<<i<<endl;
	   cout<<d->m_nbHits<<endl;
	   cout<<d->m_nbSteps<<endl;
	   cout<<d->m_trackId->size()<<endl;
			}
			
			

int hits=0;			
			
			
	//for each step,diffuse electron along y axis
	for(int j=0; j<d->m_nbHits; j++) { //for each step
	x=(*d->m_xp)[j];
	yp=(*d->m_yp)[j];
	zp=(*d->m_zp)[j];
	ed=(*d->m_energyDeposited)[j];
	z=zp*cos(sita)+yp*sin(sita);  //rotate
	y=yp*cos(sita)-zp*sin(sita);
	double electron_num=ed/energyPerElectron;//electron_num is the electron number
	double distance=fabs(y-detectorY);
	double r=r0+sqrt(distance)*r1;
	double x_start=x-r;
	double x_stop=x+r;
	double z_start=z-r;
	double z_stop=z+r;
	int x1=x_start/pixelSize;
	if(x_start<0) x1--;
	int x2=x_stop/pixelSize;
	if(x_stop>=0) x2++;
	int z1=z_start/pixelSize;
	if(z_start<0) z1--;
	int z2=z_stop/pixelSize;
	if(z_stop>=0) z2++;
	double dy_sum=0;
	for(int ix=x1;ix<=x2;ix++){
	for(int iz=z1;iz<=z2;iz++){
	double dx=ix*pixelSize+0.5*pixelSize-x;
	double dz=iz*pixelSize+0.5*pixelSize-z;
	double dy=r*r-dx*dx-dz*dz;
	if(dy<=0) continue;
	dy_sum+=sqrt(dy);
	}}
	for(int ix=x1;ix<=x2;ix++){
	for(int iz=z1;iz<=z2;iz++){
	double cx=ix*pixelSize+0.5*pixelSize;
	double cz=iz*pixelSize+0.5*pixelSize;
	double dx=cx-x;
	double dz=cz-z;
	double dy=r*r-dx*dx-dz*dz;
	int inside=0;
	if((ix*pixelSize<=x)&&(ix*pixelSize+pixelSize>x)&&(iz*pixelSize<=z)&&(iz*pixelSize+pixelSize>z))
	{inside=1;}
	if(dy<=0&&inside==0) continue; //inside=1 means all charge cloud is in one pixel cell
	
	if(dy<=0&&inside==1) {
	signalZX->Fill(cz+he,cx+he,electron_num);
		// cout<<"electron_num="<<electron_num<<"----------------"<<endl;
		


			double a=electron_num*amplification;
			double b=random.Gaus(gaussMean,gaussSigma);
			int adc_value=(a+b)/electronPerAdcCnt;
			int time_value=adcSamplingRate*(distance-r)/driftSpeed; 
			if(adc_value>threshold)
			{
            hits++;
			hitIsNoiseT->push_back(0);
			Digit_zT->push_back(iz);
			Digit_xT->push_back(ix);
			Digit_adcT->push_back(adc_value);
			TimeT->push_back(time_value);
			}

	
	
	break;}
	
	dy=sqrt(dy);
	signalZX->Fill(cz+he,cx+he,electron_num*dy/dy_sum);
			// cout<<"electron_num*dy/dy_sum="<<electron_num*dy/dy_sum<<"----------------"<<endl;

			double a=electron_num*(dy/dy_sum)*amplification;
			double b=random.Gaus(gaussMean,gaussSigma);
			int adc_value=(a+b)/electronPerAdcCnt;
			int time_value=adcSamplingRate*(distance-dy)/driftSpeed; 
			if(adc_value>threshold)
			{
            hits++;
			hitIsNoiseT->push_back(0);
			Digit_zT->push_back(iz);
			Digit_xT->push_back(ix);
			Digit_adcT->push_back(adc_value);
			TimeT->push_back(time_value);
			}

	
	}}	//loop z direction and x direction

	
	}//for each step
	
	
//add noise hit
for(int i=0;i<fakeHitNum;i++){
            hits++;
			hitIsNoiseT->push_back(0);
			Digit_zT->push_back(random.Uniform(pixelNumberZ));
			Digit_xT->push_back(random.Uniform(pixelNumberX));
			Digit_adcT->push_back(nosiePassThreshold->GetRandom());
			TimeT->push_back(random.Uniform(pixelNumberY));	
	
}	



//merge hits
int merge=0;
int hitAfterMerge=0;
for(int i=0;i<hits;i++){
	merge=0;
for(int j=0;j<hitIsNoise->size();j++){
	if(
	(*Digit_zT)[i]==(*Digit_z)[j]&&
	(*Digit_xT)[i]==(*Digit_x)[j]&&
	   (*TimeT)[i]==(*Time)[j]	
	){
	(*Digit_adc)[j]+=(*Digit_adcT)[i];
	if((*hitIsNoise)[j]!=(*hitIsNoiseT)[i]){
		(*hitIsNoise)[j]=2;
	}
	merge=1;
	break;
	}
}

if(merge==0){
	hitIsNoise->push_back((*hitIsNoiseT)[i]);
	Digit_z->push_back((*Digit_zT)[i]);
	Digit_x->push_back((*Digit_xT)[i]);
	Digit_adc->push_back((*Digit_adcT)[i]);
	Time->push_back((*TimeT)[i]);	
	hitAfterMerge++;
}	
	
}	
	
	
//if(i%10==0)
	cout<<"hits= "<<hits<<endl;
	cout<<"hitAfterMerge= "<<hitAfterMerge<<endl;
	t2->Fill(); 		
	}
    d->Clear();

	tof->Write();
	delete hitIsNoise;
	delete Digit_z;
	delete Digit_x;
	delete Digit_adc;
	delete Time;
	
	delete hitIsNoiseT;
	delete Digit_zT;
	delete Digit_xT;
	delete Digit_adcT;
	delete TimeT;
	
    return EXIT_SUCCESS;
}
