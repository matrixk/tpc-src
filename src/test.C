#include <vector>
#include <string>
#include <iostream>
#include <sstream>

TFile *tfData;
TTree *tfDataTree;
TH2D *Convolute;
TH2D *Response;

void test()
{
tfData = new TFile("wavsim_r1.root", "READ");
tfDataTree = (TTree*)(tfData->Get("r1"));
tfDataTree->ls();
tfDataTree->SetBranchAddress("Convolute", &Convolute);
tfDataTree->SetBranchAddress("Response", &Response);
tfDataTree->GetEntry(1);
cout <<"test done"<<endl;
}
