/** \file
 * Generate pixel signals according to specific pixelation schemes.
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include <math.h>
#include <algorithm>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH2I.h>
#include <TH3I.h>
#include <TMinuit.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFitter.h>
using namespace std;

extern "C" {
#include "hexlib.h"
}


typedef std::vector< std::vector<float> > vecvec;

typedef struct param 
{
	double electric_field;
	double temperature;
	double ionization_energy;
	double gas_pressure;
	double Boltzmann_const;
	double elementary_charge;
	double Fano;
	double L;
	double np;
	double time_binwidth;
	double time_min;
	double time_max;
	double p;       // pitch [mm]
	double rp;      // randomize initial location
	double noise;
} param_t;

param_t param_default = {
	.electric_field = 100, //V·cm −1
	.temperature = 293.1500, //K
	.ionization_energy = 24.8, //eV
	.gas_pressure = 10, //atm
	.Boltzmann_const = 1.380e-23, //J/K
	.elementary_charge = 1.602e-19, //C
	.Fano = 0.19,   // for Se at 10bar
	.L = 90, //cm
	.np = 300,
	.time_binwidth = 15, //ms
	.time_min = 0,
	.time_max = 30, //s
	.p = 1.0,
	.rp = 0.0,
	.noise = 10.0,
};

TH1D *currentSinglePixel;
double drift_v0_neg[2] = {5.4012619, 5.0039524};  // SeF5-,SeF6-
TH1D *responFuncHist;

Double_t fitFunction(Double_t *x, Double_t *par) //Charge(Total), Occurrence Time, Drift Length, Pedestal
{
        Double_t xx = x[0];
        Double_t binValue, doubGaus;
        int binCount = currentSinglePixel->FindBin(xx);
        binValue = 0;
        for (int seq = 0; seq < binCount; seq++)
        {

                double xseq = currentSinglePixel->GetXaxis()->GetBinCenter(seq);
		doubGaus = (1 / (sqrt(2 * TMath::Pi()) * (sqrt(0.00047033708 * par[2]) / drift_v0_neg[0]))) * (0.333 * par[0]) * TMath::Gaus(xseq-par[1],(par[2] / drift_v0_neg[0]),(sqrt(0.00047033708 * par[2]) / drift_v0_neg[0])) + (1 / (sqrt(2 * TMath::Pi()) * (sqrt(0.00047033708 * par[2]) / drift_v0_neg[1]))) * (0.666 *  par[0]) * TMath::Gaus(xseq-par[1],(par[2] / drift_v0_neg[1]),(sqrt(0.00047033708 * par[2]) / drift_v0_neg[1]));
                if (doubGaus == 0 && xseq < (par[2] / drift_v0_neg[1])) {
                        continue;
                }
                double xVal = xx - xseq;
                int xbin =  responFuncHist->FindBin(xVal);
                binValue += doubGaus * responFuncHist->GetBinContent(xbin);
        }
        return binValue + par[3];
}

int IoniFitting(const param_t *pm, TTree *t1, const char *ofname, ssize_t iStart=0, ssize_t iStop=-1)
{
	TFile *fitFile=0;
	TTree *fitTree=0;

	TFitter::SetPrecision (1e-3);
	TH2D *Convolute = 0;
	TH2D *Response = 0;
	std::string noiseFileName;
	vector<vector<double>> adcNoiseData;
        string adcNoiseDataLine;
	int adcNoiseChannelNumber;
	Double_t *xpeaks;
	Double_t *ypeaks;
	double nSignal;
	int binxNumber;
	double pedestal;
	std::vector<double> validPixelVal;
	std::vector<double> xFitPositionVal;
	std::vector<double> yFitPositionVal;
	std::vector<double> driftLengthFitVal;
	std::vector<double> chargeFitVal;
	std::vector<double> timeFitVal;
	int qPosition, rPosition;
	double xPosition, yPosition;
	double diffBin;
	auto canvas = new TCanvas("c", "c", 2000, 1000);

        //Create response function hist
        responFuncHist = new TH1D("Response Function Hist","Response Function Hist", 30000, 0, 40);
        double responFuncHistx;
        double responFuncHisty;
        for (int rfnbin = 0; rfnbin < 30000; rfnbin++) {
                responFuncHistx = responFuncHist->GetXaxis()->GetBinCenter(rfnbin);
                responFuncHisty = 1e-06 * (1.0 - 1.0 / (exp((responFuncHistx-0.5)/1e-05)+1.0)) * (exp(-(responFuncHistx-0.5)/5.0) + (1.0 - 1.0) * exp(-(responFuncHistx-0.5)/5.0));
                responFuncHist->SetBinContent(rfnbin, responFuncHisty);
        }

	//read original root file
	t1->SetBranchAddress("Convolute", &Convolute);
	t1->SetBranchAddress("Response", &Response);
	t1->SetBranchAddress("NoiseChannel", &adcNoiseChannelNumber);

	//read noise adc data to 2D array
	noiseFileName = "inside_flat3_adc.dat";
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

	//create new root file
	fitFile = new TFile(ofname, "RECREATE");
	fitFile->cd();
        fitTree = new TTree("f1", "FittingResults");
        fitTree->Branch("n_signal", &nSignal);
        fitTree->Branch("Vaild_Pixel", &validPixelVal);
        fitTree->Branch("x_position", &xFitPositionVal);
        fitTree->Branch("y_position", &yFitPositionVal);
        fitTree->Branch("Drift_Length", &driftLengthFitVal);
        fitTree->Branch("Charge", &chargeFitVal);
        fitTree->Branch("Occur_Time", &timeFitVal);

        std::cout << "******************** Ready ************************" << std::endl;

	TSpectrum *s = new TSpectrum();
	for (int eventNumber = 20; eventNumber < 30; eventNumber++) {
		t1->GetEntry(eventNumber);
		nSignal = 0;
		pedestal = adcNoiseData[10][adcNoiseChannelNumber];

		binxNumber = Response->GetNbinsX();

		for (int pixelNumber = 1; pixelNumber <= binxNumber; pixelNumber++) {
                        currentSinglePixel = Response->ProjectionY("Single Pixel", pixelNumber, pixelNumber);
			Int_t xbinCountCurrent = currentSinglePixel->GetNbinsX();
			double xbinWidthCurrent = currentSinglePixel->GetBinWidth(1);

			//Difference
			TH1D *diffCurrent = (TH1D*)currentSinglePixel->Clone("diffCurrent");
			diffCurrent->Reset();
			int diffStep = 8;
			for (int xbinLoop = diffStep; xbinLoop < (xbinCountCurrent - diffStep - 1); xbinLoop++) {
				diffBin = (currentSinglePixel->GetBinContent(xbinLoop) - currentSinglePixel->GetBinContent(xbinLoop - diffStep)) / (xbinWidthCurrent * diffStep);
				diffCurrent->SetBinContent(xbinLoop, diffBin);
			}

			//Finding peaks
			TSpectrum *s = new TSpectrum();
                        Int_t nfound = s->Search(diffCurrent,5,"",0.4);

			//Skip?
                        if (nfound != 2) {
                                validPixelVal.push_back(0);
                                xFitPositionVal.push_back(0);
                                yFitPositionVal.push_back(0);
                                driftLengthFitVal.push_back(0);
                                chargeFitVal.push_back(0);
                                timeFitVal.push_back(0);
                                continue;
                        }

                        xpeaks = s->GetPositionX();
                        ypeaks = s->GetPositionY();
                        std::sort(xpeaks,xpeaks+nfound);
                        std::sort(ypeaks,ypeaks+nfound);
		
                        //Calculate X/Y Position of Pixel
                        hex_l2qr(pixelNumber, &qPosition, &rPosition);
                        hex_qr2xy(1, qPosition, rPosition, &xPosition, &yPosition);
                        xFitPositionVal.push_back(xPosition);
                        yFitPositionVal.push_back(yPosition);

			//Calculate initial value
                        double chargePar, driftLengthPar, timePar;
			chargePar = 2250 * ((ypeaks[0] - pedestal) +  (ypeaks[1] - pedestal));
                        driftLengthPar = (xpeaks[nfound-1] - xpeaks[nfound-2]) / ((1 / drift_v0_neg[1])-(1 / drift_v0_neg[0]));
			timePar = 0.5 * ((driftLengthPar / drift_v0_neg[0]) - xpeaks[nfound-2] + (driftLengthPar / drift_v0_neg[1]) - xpeaks[nfound-1]);

			//SetBinError
                        for(int i=0; i<4000; i++){
                                currentSinglePixel->SetBinError(i+1,0.0002764);
                        }

			//Draw Initial Plot
			TF1 *fitFunc = new TF1("fitFunc", fitFunction, 0, 60, 4);
                        fitFunc->SetParNames("Charge(Total)", "Occurrence Time", "Drift Length", "Pedestal");
			fitFunc->SetParameters(chargePar, timePar , driftLengthPar, pedestal);
			canvas->cd();
			fitFunc->Draw();
			currentSinglePixel->Draw("same");
			canvas->SaveAs("canvas.png");

			//Fitting
			fitFunc->SetParLimits(0, 0.5*chargePar, 2*chargePar);
			fitFunc->SetParLimits(1, -0.5, 0.5);
			fitFunc->SetParLimits(2, driftLengthPar-10, driftLengthPar+10);
			fitFunc->SetParLimits(3, pedestal-0.001, pedestal+0.001);
                        fitFunc->SetNpx(6000);
			std::cout << "****************** Start Fitting ******************" << std::endl;
                        currentSinglePixel->Fit("fitFunc");
			canvas->SaveAs("canvas_func_re.png");

			//Get Final Parameters
                        double chargeVal, timeVal, driftLengthVal, pedestalVal;
                        chargeVal = fitFunc->GetParameter(0);
                        timeVal = fitFunc->GetParameter(1);
			driftLengthVal = fitFunc->GetParameter(2);
                        pedestalVal = fitFunc->GetParameter(3);

			//Failed?
			TString fitStat;
			fitStat = gMinuit->fCstatu;

			if (fitStat != "CONVERGED ") {
				validPixelVal.push_back(0);
                                xFitPositionVal.push_back(0);
                                yFitPositionVal.push_back(0);
                                driftLengthFitVal.push_back(0);
                                chargeFitVal.push_back(0);
                                timeFitVal.push_back(0);
				std::cout << "FIT FAILED, SKIP CURRENT PIXEL" << std::endl;
				s = NULL;
				continue;
			}

                        //Store Fitting Value
			nSignal += 1;
			validPixelVal.push_back(1);
			xFitPositionVal.push_back(xPosition);
                        yFitPositionVal.push_back(yPosition);
                        driftLengthFitVal.push_back(driftLengthVal);
                        chargeFitVal.push_back(chargeVal);
                        timeFitVal.push_back(timeVal);
			std::cout << "EVENT\tPIXEL\tX POSI\tY POSI\tChargeFitVal\tOccurTimeFitVal\tDriftLengthFitVal\tPedestalVal" << std::endl;
			std::cout << eventNumber << "\t" << pixelNumber << "\t" << xPosition << "\t" << yPosition << "\t" << chargeVal << "\t" << timeVal << "\t" << driftLengthVal << "\t" << pedestalVal << std::endl;
			s = NULL;
		}
		fitTree->Fill();
		validPixelVal.clear();
                xFitPositionVal.clear();
                yFitPositionVal.clear();
                driftLengthFitVal.clear();
                chargeFitVal.clear();
                timeFitVal.clear();
	}
	fitTree->Write();
	fitFile->CurrentFile()->Close();
}

#ifndef __CINT__

void print_usage(const param_t *pm)
{
	printf("Usage:\n");
}

/** Main entry if this file is compiled outside of root */
int main(int argc, char **argv)
{
	param_t pm;
	
	char *sRootFname;
	std::string ofname;
	std::stringstream ss;
	ssize_t iStart=0, iStop=-1, sz;
	
	TFile *tfData;
	TTree *tfDataTree;

	memcpy(&pm, &param_default, sizeof(pm));
	argc -= optind;
	argv += optind;
	if(argc<1 || argc>=2) {
		print_usage(&pm);
		return EXIT_FAILURE;
	}

	sRootFname = argv[0];
	ofname.assign(sRootFname);
	sz = ofname.size();
	ofname.resize(sz+3);
	ofname.replace(sz-5, 8, "_f1.root");

	// signal			
	tfData = new TFile(sRootFname, "READ");
	if(tfData->IsZombie()) {
		std::cerr << "Error opening input file " << sRootFname << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "******************** Signal ************************" << std::endl;
	tfData->ls();
	tfDataTree = (TTree*)(tfData->Get("r1"));
	tfDataTree->ls();
	std::cout << "r1 has " << tfDataTree->GetEntries() << " entries." << std::endl;

	IoniFitting(&pm, tfDataTree, ofname.c_str(), iStart, iStop);

	tfData->Close();

	return EXIT_SUCCESS;
}
#endif
