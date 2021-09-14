/** \file
 * Simulate the digitized CSA output waveform.
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH2I.h>
#include <TH3I.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TUnuran.h>
#include <TUnuranMultiContDist.h>
#include <TMath.h>
using namespace std;

extern "C" {
#include "spectrarand.h"
#include "hexlib.h"
}

#if defined(__MAKECINT__) || defined(__CINT__)
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<long> >+;
#endif

#ifndef MIN
#define MIN(a,b) ((a)>(b)?(b):(a))
#endif

typedef std::vector< std::vector<float> > vecvec;

typedef struct param 
{
	size_t n;	   // number of samples
	double fs;	  // sampling frequency
	int l;		  // spiral coordinate
	double t0;	  // CSA parameters
	double tr;
	double tau;
	double a;
	double b;

} param_t;

param_t param_default = {
	.n  = 1000,
	.fs = 1.0e7,
	.l  = 0,
	.t0 = 0.5,   // CSA parameters
	.tr = 1e-05,
	.tau = 5,
	.a  = 0.7,
	.b  = 1.0
};

TRandom3 *tr=0;
TH1D *responFuncHist;
double drift_mu0_neg[2] = {0.466, 0.503};
double electric_field = 100.; //V·cm −1
double temperature = 293.1500; //K
double gas_pressure = 10.; //atm
double Boltzmann_const = 1.380e-23; //J/K
double elementary_charge = 1.602e-19; //C
double L = 90.; //cm
double noiseLevel = 0.783;

extern "C" double rand0_1(void)
{
	if(tr)
		return tr->Uniform();
	else
		return drand48();
}

void AddNoise(double *originIn, double *addOut, double *arrayIn, double *arrayOut, int binyNumber, int adcDataSize)
{
	fftw_complex *outCpx;

	fftw_plan fft;
	fftw_plan ifft;
	outCpx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * adcDataSize);
	arrayOut = (double *) malloc((binyNumber*2)*sizeof(double));

	fft = fftw_plan_dft_r2c_1d(adcDataSize, arrayIn, outCpx, FFTW_ESTIMATE);  //Setup fftw plan for fft
	ifft = fftw_plan_dft_c2r_1d(binyNumber*2, outCpx, arrayOut, FFTW_ESTIMATE);   //Setup fftw plan for ifft
	fftw_execute(fft);
	fftw_execute(ifft);

	for(int i = 0; i < (binyNumber*2); i++)
        {
                addOut[i] = originIn[i] + arrayOut[i] / adcDataSize;
                //std::cout << "addOut = "<< addOut[i] << std::endl;
        }
	fftw_destroy_plan(fft);
	fftw_destroy_plan(ifft);
	fftw_free(outCpx);
}
#ifndef __CINT__

int main(int argc, char **argv)
{
	//param_t pm;
	char *sourceFileName;
	std::string changeFileName;
	std::string noiseFileName;
	ssize_t changeFileNameSize;

	TFile *afterPix;
	TTree *afterPixTree;
	TH2D *ArrvingTime_neg = 0;
	TH2D *ArrvingTime_pos = 0;

	TH2D *convHist;
	TH2D *responNoiseHist;
	TFile *responFile;
	TTree *responTree;

	TH1D *responFuncHist;

	//double *originSignalNeg;
	double *convolveResultsNeg, *fftResults, *addNoiseResultsNeg;

	vector<vector<double>> adcNoiseData;
	string adcNoiseDataLine;

	int eventNumber, pixelNumber, i, seq, y;
	int startBinSinglePixel;
	double binValue, signalValue, xVal;
	int xbin;
	int adcNoiseChannelNumber;


	sourceFileName = argv[1];
	changeFileName.assign(sourceFileName);
	changeFileNameSize = changeFileName.size();
	changeFileName.resize(changeFileNameSize+3);
	changeFileName.replace(changeFileNameSize-5, 8, "_r1.root");

	noiseFileName = argv[2];

	std::cout << changeFileNameSize << std::endl;
	std::cout << sourceFileName <<" "<< changeFileName <<" "<< noiseFileName << std::endl;

	//read original data
	afterPix = new TFile(sourceFileName, "READ");
	if(afterPix->IsZombie()) {
		std::cerr << "Error opening input file " << afterPix << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "******************** Signal ************************" << std::endl;

	afterPix->ls();
	afterPixTree = (TTree*)(afterPix->Get("p1"));
	afterPixTree->ls();
	std::cout << "p1 has " << afterPixTree->GetEntries() << " entries." << std::endl;

	afterPixTree->SetBranchAddress("ArrvingTime_neg", &ArrvingTime_neg);
	afterPixTree->SetBranchAddress("ArrvingTime_pos", &ArrvingTime_pos);

	afterPixTree->GetEntry(0);
	
	//Get Negtive Bin Range
        TAxis *xaxisNeg = ArrvingTime_neg->GetXaxis();
        TAxis *yaxisNeg = ArrvingTime_neg->GetYaxis();
        Int_t binxNumberNeg = ArrvingTime_neg->GetNbinsX();
        Int_t binyNumberNeg = ArrvingTime_neg->GetNbinsY();
	Double_t binxUpNeg = xaxisNeg->GetBinUpEdge(binxNumberNeg);
	Double_t binyUpNeg = yaxisNeg->GetBinUpEdge(binyNumberNeg);
        Double_t binxWidthNeg = xaxisNeg->GetBinWidth(0);
        Double_t binyWidthNeg = yaxisNeg->GetBinWidth(0);
        std::cout << "[NEG] X axis(pixel) has " << binxNumberNeg << " bins, binwidth "<< binxWidthNeg << ", Up Edge " << binxUpNeg << std::endl;
        std::cout << "[NEG] Y axis(Time) has " << binyNumberNeg << " bins, binwidth "<< binyWidthNeg << ", Up Edge " << binyUpNeg << std::endl;

	//create new root file
	convHist = new TH2D("Convolve","Convolve;Pixel number;Time(s)", binxNumberNeg, 0, binxUpNeg, binyNumberNeg*2, 0, binyUpNeg*2);
        responNoiseHist = new TH2D("Response","Response;Pixel number;Time(s)", binxNumberNeg, 0, binxUpNeg, binyNumberNeg*2, 0, binyUpNeg*2);
	responFile = new TFile("wavsim_r1.root", "RECREATE");
	responFile->cd();
	responTree = new TTree("r1", "Response");
	responTree->Branch("Convolute", "TH2D", &convHist);
	responTree->Branch("Response", "TH2D", &responNoiseHist);
	responTree->Branch("NoiseChannel", &adcNoiseChannelNumber, "adcNoiseChannel/I");
	std::cout << "******************** Ready ************************" << std::endl;
	std::cout << "X bins = " << convHist->GetNbinsX() << "; Y bins = " << convHist->GetNbinsY() << std::endl;

	//Create response function hist
        responFuncHist = new TH1D("Response Function Hist","Response Function Hist", 30000, 0, 40);
        double responFuncHistx;
        double responFuncHisty;
        for (int rfnbin = 0; rfnbin < 30000; rfnbin++) {
                responFuncHistx = responFuncHist->GetXaxis()->GetBinCenter(rfnbin);
                //std::cout << "responFuncHistx = "<< responFuncHistx << std::endl;
                responFuncHisty = 1e-06 * (1.0 - 1.0 / (exp((responFuncHistx-0.5)/1e-05)+1.0)) * (exp(-(responFuncHistx-0.5)/5.0) + (1.0 - 1.0) * exp(-(responFuncHistx-0.5)/5.0));
                responFuncHist->SetBinContent(rfnbin, responFuncHisty);
		//std::cout << "ResponseFunc: "<< responFuncHist->GetBinContent(rfnbin) << std::endl;
        }
	std::cout << "Response Function Read In Successful" << std::endl;

	//read noise adc data to 2D array
        ifstream adcNoiseFile(noiseFileName);
        while (getline(adcNoiseFile, adcNoiseDataLine))
        {
                adcNoiseData.push_back(vector<double>());
                stringstream ss(adcNoiseDataLine);
                double adcNoiseValue;
                while (ss >> adcNoiseValue)
                {
                        adcNoiseData.back().push_back(adcNoiseValue);
                }
        }
	std::cout << "Noise ADC Data  Read In Successful" << std::endl;
	double adcNoiseSingleChannel[adcNoiseData.size()];

	//Start to process
	for (eventNumber = 0; eventNumber <= afterPixTree->GetEntries(); eventNumber++) {
	//for (eventNumber = 1; eventNumber <= 1; eventNumber++) {
		afterPixTree->GetEntry(eventNumber);
		convHist->Reset();
		responNoiseHist->Reset();
		std::cout << "Now Process Event "<< eventNumber << std::endl;
		for (pixelNumber = 0; pixelNumber <= binxNumberNeg; pixelNumber++) {
		//for (pixelNumber = 0; pixelNumber <= 1; pixelNumber++) {
			//Set Start Bin
			startBinSinglePixel = binyNumberNeg;
			for(i=0; i<=binyNumberNeg; i++) {
				if (ArrvingTime_neg->GetBinContent(ArrvingTime_neg->GetBin(pixelNumber, i)) !=0) {
					startBinSinglePixel = i;
					break;
				}
			}
			//std::cout << "startBinSinglePixel = "<< startBinSinglePixel << std::endl;

			//Allocates array memory
			convolveResultsNeg = (double*)calloc(binyNumberNeg*2, sizeof(double));
			
			//Convolute
			if (startBinSinglePixel != binyNumberNeg) {
				for(i=startBinSinglePixel; i<=binyNumberNeg*2; i++) {
					//std::cout << "x = "<< responseXVal << std::endl;
					binValue = 0;
					for (seq = startBinSinglePixel; seq <= i; seq++) {
						signalValue = ArrvingTime_neg->GetBinContent(pixelNumber, seq);
						xVal = ArrvingTime_neg->GetYaxis()->GetBinCenter(i) - ArrvingTime_neg->GetYaxis()->GetBinCenter(seq);
						xbin = responFuncHist->FindBin(xVal);
						binValue += signalValue * responFuncHist->GetBinContent(xbin);
					}
				//std::cout << "Bin = "<< i << ";Signal = " << ArrvingTime_neg->GetBinContent(ArrvingTime_neg->GetBin(pixelNumber, i)) << "; Value = " << binValue << std::endl;
				convHist->SetBinContent(pixelNumber, i, binValue);
				convolveResultsNeg[i] = binValue;
				}
			}
			
			//std::cout << "Convoluted - Pixel "<< pixelNumber << std::endl;

			//Random adc Noise Channel
			TRandom3 *adcNoiseChannelRm = new TRandom3();
			adcNoiseChannelNumber = (int)floor(adcNoiseChannelRm->Rndm() * 20);

			//double adcNoiseSingleChannel[adcNoiseData.size()];
        		for (y = 0; y < (int)adcNoiseData.size(); y++) {
                		adcNoiseSingleChannel[y] = adcNoiseData[y][adcNoiseChannelNumber];
        		}

			//Add Noise
			fftResults = (double*)calloc(binyNumberNeg*2, sizeof(double));
			addNoiseResultsNeg = (double*)calloc(binyNumberNeg*2, sizeof(double));
			AddNoise(convolveResultsNeg, addNoiseResultsNeg, adcNoiseSingleChannel, fftResults, binyNumberNeg, adcNoiseData.size());

			//Set Content
			for (i=0; i<=binyNumberNeg*2; i++) {
				//std::cout << "Final Value = "<< addNoiseResultsNeg[i]  << std::endl;
				responNoiseHist->SetBinContent(pixelNumber, i, addNoiseResultsNeg[i]);
			}
			delete adcNoiseChannelRm;
			free(convolveResultsNeg);
			convolveResultsNeg = NULL;
			free(fftResults);
			fftResults = NULL;
			free(addNoiseResultsNeg);
			addNoiseResultsNeg = NULL;
		}
		responTree->Fill();
	}
	responTree->Write();
	responFile->CurrentFile()->Close();
	delete responFile;

	delete convHist;
	delete responNoiseHist;
	
	afterPix->Close();
	delete afterPix;

	return EXIT_SUCCESS;
}

#endif
