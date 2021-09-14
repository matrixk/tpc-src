/** \file
 * Generate pixel signals according to specific pixelation schemes.
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include <getopt.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

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
	.noise = 10.0
};

int IoniImage(const param_t *pm, TTree *t1, const char *ofname, ssize_t iStart=0, ssize_t iStop=-1)
{
	TFile *tfp=0;
	TTree *tp1=0;

	std::vector<int> *parentId=0;
	std::vector<double> *xp=0, *yp=0, *zp=0, *ed=0; // must be initialized to NULL
	int nIonTot, nPixTot;
	double sigTot;
	TRandom3 *tr = new TRandom3;
	TH2D *ArrvingTime_neg;
	TH2D *ArrvingTime_pos;
	TH1D *lH_neg;
	TH1D *lH_pos;
	TH1D *diffusion_neg_x;
	TH1D *diffusion_neg_y;
	TH1D *diffusion_pos_x;
	TH1D *diffusion_pos_y;
	int nIon;
	double mIon, sIon;
	double zbin;
	double drift_mu0_neg[2] = {0.466, 0.503};
	double drift_mu0_pos[6] = {0.466, 0.503, 0.546, 0.599, 0.663, 0.746};
	double Dt_neg[2], Dl_neg[2], distance_trans_neg[2], distance_long_neg[2];
	double Dt_pos[6], Dl_pos[6], distance_trans_pos[6], distance_long_pos[6];
	double drift_v_neg[2];
	double drift_v_pos[6];
	double distance_trans_neg_x_rm, distance_trans_neg_y_rm, distance_long_neg_rm;
	double distance_trans_pos_x_rm, distance_trans_pos_y_rm, distance_long_pos_rm;
	int particle_type_rm_pos, particle_type_rm_neg;
	double particle_select_pos, particle_select_neg;
	double drift_time;
	int spiral_number;
	int i, j, k;
	int q, r;

	//count zbin
	zbin = (pm->time_max - pm->time_min) * 1e03 / pm->time_binwidth;
	std::cout << "zbin = "<< zbin << std::endl;

	ArrvingTime_neg = new TH2D("ArrvingTime_neg","Arrving Time(Negative ion, 20cm);Pixel number;Time(s)", pm->np, 0, pm->np, zbin, pm->time_min, pm->time_max);
	ArrvingTime_pos = new TH2D("ArrvingTime_pos","Arrving Time(Positive ion, 20cm);Pixel number;Time(s)", pm->np, 0, pm->np, zbin, pm->time_min, pm->time_max);
	diffusion_neg_x = new TH1D("diffusion_neg_x","diffusion of negtive ion (x);dx;count", 1000, -6, 6);
	diffusion_neg_y = new TH1D("diffusion_neg_y","diffusion of negtive ion (y);dy;count", 1000, -6, 6);
	diffusion_pos_x = new TH1D("diffusion_pos_x","diffusion of positive ion (x);dx;count", 1000, -6, 6);
	diffusion_pos_y = new TH1D("diffusion_pos_y","diffusion of positive ion (y);dy;count", 1000, -6, 6);

	//read original root file
	t1->SetBranchAddress("parentid", &parentId);
	t1->SetBranchAddress("ed", &ed);
	t1->SetBranchAddress("xp", &xp);
	t1->SetBranchAddress("yp", &yp);
	t1->SetBranchAddress("zp", &zp);

	//create new root file
	tfp = new TFile(ofname, "RECREATE");
	tfp->cd();
	tp1 = new TTree("p1", "IoniImage");
	tp1->Branch("nIonTot", &nIonTot, "nIonTot/I");
	tp1->Branch("nPixTot", &nPixTot, "nPixTot/I");
	tp1->Branch("sigTot", &sigTot, "sigTot/D");
	tp1->Branch("ArrvingTime_neg", "TH2D", &ArrvingTime_neg);
	tp1->Branch("ArrvingTime_pos", "TH2D", &ArrvingTime_pos);
	tp1->Branch("diffusion_neg_x", "TH1D", &diffusion_neg_x);
	tp1->Branch("diffusion_neg_y", "TH1D", &diffusion_neg_y);
	tp1->Branch("diffusion_pos_x", "TH1D", &diffusion_pos_x);
	tp1->Branch("diffusion_pos_y", "TH1D", &diffusion_pos_y);

	if(iStart < 0) iStart = 0;
	if(iStart >= t1->GetEntries()) iStart = t1->GetEntries()-1;
	if(iStop < 0) iStop = t1->GetEntries();
	if(iStop > t1->GetEntries()) iStop = t1->GetEntries();
	std::cout << "iStart=" << iStart <<" iStop=" << iStop << std::endl;

	for (i=0; i<2; i++) {
		//std::cout << "drift mu0 ="<< drift_mu0_neg[i] << std::endl;
		//std::cout << "gas_pressure = "<< pm->gas_pressure << " electric_field = " << pm->electric_field << std::endl;
		drift_v_neg[i] = (1. / pm->gas_pressure) * (pm->temperature / 273) * drift_mu0_neg[i] * pm->electric_field;
		//std::cout << "drift velocity(neg) = "<< drift_v_neg[i] << std::endl;

		//diffusion
		Dt_neg[i] = (drift_mu0_neg[i] * pm->Boltzmann_const * pm->temperature) / (pm->elementary_charge * pm->gas_pressure);
		Dl_neg[i] = (drift_mu0_neg[i] * pm->Boltzmann_const * pm->temperature) / (pm->elementary_charge * pm->gas_pressure);
		distance_trans_neg[i] = sqrt((2. * Dt_neg[i] * pm->L) / drift_v_neg[i]);
		//std::cout << "distance_trans_neg = "<< distance_trans_neg[i] << " Dt_neg = "<< Dt_neg[i] << " pm->L = " << pm->L << " drift_v_neg = " << drift_v_neg[i] << std::endl;
		distance_long_neg[i] = sqrt((2. * Dl_neg[i] * pm->L) / drift_v_neg[i]);
		//std::cout << "distance_trans_neg = "<< distance_trans_neg[i] << " Dl_neg = "<< Dl_neg[i] << " Tran = " << distance_trans_neg[i] << " Longt = " << distance_long_neg[i] << std::endl;
	}

	for (i=0; i<6; i++) {
		//std::cout << "drift mu0 ="<< drift_mu0_pos[i] << std::endl;
		//std::cout << "gas_pressure = "<< pm->gas_pressure << " electric_field = " << pm->electric_field << std::endl;
		drift_v_pos[i] = (1. / pm->gas_pressure) * (pm->temperature / 273) * drift_mu0_pos[i] * pm->electric_field;
		//std::cout << "drift velocity(pos) = "<< drift_v_pos[i] << std::endl;

		//diffusion
		Dt_pos[i] = (drift_mu0_pos[i] * pm->Boltzmann_const * pm->temperature) / (pm->elementary_charge * pm->gas_pressure);
		Dl_pos[i] = (drift_mu0_pos[i] * pm->Boltzmann_const * pm->temperature) / (pm->elementary_charge * pm->gas_pressure);
		distance_trans_pos[i] = sqrt(2. * Dt_pos[i] * pm->L / drift_v_pos[i]);
		//std::cout << "distance_trans_pos = "<< distance_trans_pos[i] << " Dt_pos = "<< Dt_pos[i] << " pm->L = " << pm->L << " drift_v_pos = " << drift_v_pos[i] << std::endl;
		distance_long_pos[i] = sqrt(2. * Dl_pos[i] * pm->L / drift_v_pos[i]);
		//std::cout << "distance_long_pos = "<< distance_long_pos[i] << " Dl_pos = "<< Dl_pos[i] << " Tran = " << distance_trans_pos[i] << " Longt = " << distance_long_pos[i] << std::endl;
	}

	for(i=0; i<iStop; i++) {
		t1->GetEntry(i);
		std::cout << "Get entry "<< i << std::endl;
		ArrvingTime_pos->Reset();
		ArrvingTime_neg->Reset();
		std::cout << "Start to process event "<< i << std::endl;

		for(j=0; j<(ssize_t)ed->size(); j++) {
			mIon = (*ed)[j] * 1000.0 / pm->ionization_energy;
			sIon = std::sqrt(pm->Fano * mIon);
			nIon = (int)tr->Gaus(mIon, sIon);
			//std::cout << "nIon = " << nIon << std::endl;
			if(nIon<0) nIon = 0;
			nIonTot += nIon;
			if(nIon == 0) continue;
			for (k=0; k<nIon; k++) {
			
				//particle random
				particle_select_pos = tr->Rndm();
				if(particle_select_pos<0.510){
					particle_type_rm_pos = 0;
				}
				else if(particle_select_pos>=0.510 && particle_select_pos<0.541){
					particle_type_rm_pos = 1;
				}
				else if(particle_select_pos>=0.541 && particle_select_pos<0.668){
					particle_type_rm_pos = 2;
				}
				else if(particle_select_pos>=0.668 && particle_select_pos<0.719){
					particle_type_rm_pos = 3;
				}
				else if(particle_select_pos>=0.719 && particle_select_pos<0.780){
					particle_type_rm_pos = 4;
				}
				else {
					particle_type_rm_pos = 5;
				}
				//std::cout << "particle_type_pos = "<< particle_type_rm_pos << std::endl;

				particle_select_neg = tr->Rndm();
				if(particle_select_neg<0.666){
					particle_type_rm_neg = 0;
				}
				else {
					particle_type_rm_neg = 1;
				}
				//std::cout << "particle_type_neg = "<< particle_type_rm_neg << std::endl;

				//positive Ion
				//gauss random
				distance_trans_pos_x_rm = tr->Gaus(0,distance_trans_pos[particle_type_rm_pos]);
				distance_trans_pos_y_rm = tr->Gaus(0,distance_trans_pos[particle_type_rm_pos]);
				distance_long_pos_rm = tr->Gaus(0,distance_long_pos[particle_type_rm_pos]);

				//drifting
				drift_time = (pm->L + distance_long_pos_rm) / drift_v_pos[particle_type_rm_pos];
				spiral_number = hex_xy2qr(pm->p, (*xp)[j] + distance_trans_pos_x_rm, (*yp)[j] + distance_trans_pos_y_rm, &q, &r);

				//std::cout << "drift time = " << drift_time << " spiral number = " << spiral_number << std::endl;
				//std::cout << "x = " << (*xp)[j] << " diff = " << distance_trans_neg[particle_type_rm_neg] << std::endl;
				//std::cout << "y = " << (*xp)[j] << " diff = " << distance_trans_neg[particle_type_rm_neg] << std::endl;
	
				//filling
				ArrvingTime_pos->Fill(spiral_number, drift_time);
				diffusion_pos_x->Fill(distance_trans_pos_x_rm);
				diffusion_pos_y->Fill(distance_trans_pos_y_rm);

				//nagetive Ion
				//gauss random
				distance_trans_neg_x_rm = tr->Gaus(0,distance_trans_neg[particle_type_rm_neg]);
				distance_trans_neg_y_rm = tr->Gaus(0,distance_trans_neg[particle_type_rm_neg]);
				//std::cout << "distance_trans_neg_y_rm = "<< distance_trans_neg[particle_type_rm_neg] << std::endl;
				distance_long_neg_rm = tr->Gaus(0,distance_long_neg[particle_type_rm_neg]);

				//drifting
				drift_time = (pm->L + distance_long_neg_rm) / drift_v_neg[particle_type_rm_neg];
				spiral_number = hex_xy2qr(pm->p, (*xp)[j] + distance_trans_neg_x_rm, (*yp)[j] + distance_trans_neg_y_rm, &q, &r);
				//std::cout << "drift time = " << drift_time << " spiral number = " << spiral_number << std::endl;
				//std::cout << "x = " << (*xp)[j] << " diff = " << distance_trans_neg[particle_type_rm_neg] << std::endl;
				//std::cout << "y = " << (*xp)[j] << " diff = " << distance_trans_neg[particle_type_rm_neg] << std::endl;

				//filling
				ArrvingTime_neg->Fill(spiral_number, drift_time);
				diffusion_neg_x->Fill(distance_trans_neg_x_rm);
				diffusion_neg_y->Fill(distance_trans_neg_y_rm);
			}
		}
		std::cout << "Filling done" << std::endl;
		lH_pos = ArrvingTime_pos->ProjectionX();
		lH_neg = ArrvingTime_neg->ProjectionX();
		nPixTot = 0; sigTot = 0.0;
		for(j=0; j<lH_pos->GetNbinsX(); j++) {
			double c = lH_pos->GetBinContent(j+1);
			if(c > 0.5) {
			nPixTot++;
			sigTot += tr->Gaus(c, 10);
			}
		}
		for(j=0; j<lH_neg->GetNbinsX(); j++) {
			double d = lH_neg->GetBinContent(j+1);
			if(d > 0.5) {
			nPixTot++;
			sigTot += tr->Gaus(d, 10);
			}
		}
		delete lH_pos;
		delete lH_neg;
		tp1->Fill();
	}
	tp1->Write();
	std::cout << "Write done" << std::endl;
	tfp->CurrentFile()->Close();
	delete tfp;
	
	delete ArrvingTime_neg;
	delete ArrvingTime_pos;
	delete tr;
	return i;
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
	
	TFile *tfs;
	TTree *t1s;

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
	ofname.replace(sz-5, 8, "_p1.root");

	// signal			
	tfs = new TFile(sRootFname, "READ");
	if(tfs->IsZombie()) {
		std::cerr << "Error opening input file " << sRootFname << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "******************** Signal ************************" << std::endl;
	tfs->ls();
	t1s = (TTree*)(tfs->Get("t1"));
	t1s->ls();
	std::cout << "t1 has " << t1s->GetEntries() << " entries." << std::endl;

	IoniImage(&pm, t1s, ofname.c_str(), iStart, iStop);

	tfs->Close();
	delete tfs;

	return EXIT_SUCCESS;
}
#endif
