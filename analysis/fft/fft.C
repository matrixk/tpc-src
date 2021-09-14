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
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <math.h>
#include <TVirtualFFT.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TMath.h>

#include <iostream>
#include <string>
#include <vector>
#include <TComplex.h>
#include "signal.h"

using namespace std;

int help=0;
int gargc;
char** gargv;
int para(char *a);
int para(char *a)
{
    if(help==1)cout<<a<<" ";
    int i;
    for(i=0;i<gargc;i++) {
        if(strcmp(a,gargv[i])==0) return i;
    }
    return -1;
}

int main(int argc, char **argv)
{
    help=0;
    gargc=argc;
    gargv=argv;
    int pid;
    pid=para("-help");if(pid>0){help=1;};

    int M=1024;
    double *re_d=new double[M];
    double *im_d=new double[M];
    TFile *rf=new TFile("rf.root", "RECREATE");

    TH1D *  signal_re=new TH1D("signal_re_name","signal_re_title",M,0,M);
    TH1D *  signal_im=new TH1D("signal_im_name","signal_im_title",M,0,M);

    TH1D *  hist_re=new TH1D("hist_re_name","hist_re_title",M,0,M);
    TH1D *  hist_im=new TH1D("hist_im_name","hist_im_title",M,0,M);
    TH1D *  hist_mag=new TH1D("hist_mag_name","hist_mag_title",M,0,M);
    TH1D *  hist_phase=new TH1D("hist_phase_phase","hist_mag_title",M,0,M);

    signal_2(re_d,im_d,M);
 	//loop pixel(N samples per pixel)
    TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &M, "C2CFORWARD");
    fftr2c->SetPointsComplex(re_d,im_d);
    fftr2c->Transform();

    double re, im;
    double mag=0;
    double phase=0;
    double sum_mag=0;  
    for (int i=0;i<M;i++)
    {
        cout<<re_d[i]<<endl;
        fftr2c->GetPointComplex(i,re,im);
	    mag=sqrt(re*re+im*im);
		phase=atan2(im,re);
		signal_re->SetBinContent(i+1,re_d[i]);
		signal_im->SetBinContent(i+1,im_d[i]);
		hist_re->SetBinContent(i+1,re);
		hist_im->SetBinContent(i+1,im);
		hist_mag->SetBinContent(i+1,mag);
		hist_phase->SetBinContent(i+1,phase);
    }   
   
    rf->Write();

    delete [] re_d;
    delete [] im_d;

    return EXIT_SUCCESS;
}
