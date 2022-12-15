#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooBernstein.h"
#include "RooAbsPdf.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include <time.h>
#include "TBox.h"


using namespace RooFit;


//int AllFits(std::string filename,std::string path, double mass,double ctau,int channel, std::string wd,std::string sel,double lumiRatio, bool binned, bool blinded ){


void Make_x2(std::string ToysFile, std::string WorkspaceFile, std::string FitFile, int NToys, std::string CatName, std::string Cat ){

TFile *f = new TFile((WorkspaceFile).c_str());
TFile *toys= new TFile((ToysFile).c_str());
TFile *limitfile = new TFile((FitFile).c_str());

RooWorkspace *w = (RooWorkspace *)f->Get("w");
RooRealVar *x = w->var((CatName).c_str());

TTree* limittree = (TTree*) limitfile->Get("limit");
double r=0;
limittree->SetBranchStatus("*",0);
limittree->SetBranchStatus("r",1);
limittree->SetBranchAddress("r", &r);


int varBins_all [88] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058,1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869,5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250, 7500, 7750, 8000};

vector<int> varBins;
TCanvas* c= new TCanvas("c","c", 900,600);


for (int i=0; i<NToys; i++){
    RooDataSet *dataset= (RooDataSet* ) toys->Get(("toys/toy_"+to_string(i+1)).c_str());
    RooDataHist *BinnedHist = dataset->binnedClone("BinnedClone","BinnedClone");
    const RooArgSet *argset = dataset->get();

    RooRealVar* toyvar = (RooRealVar* ) argset->find((CatName).c_str());
    double_t var_min_set = toyvar->getMin();
    double_t var_max_set = toyvar->getMax();


    RooAbsBinning& b = toyvar->getBinning();
    int NvarBins = toyvar->getBins();
    if (i==0){
        for( unsigned int j=0;j<sizeof(varBins_all); j++){
            if (varBins_all[j]<var_min_set) continue;
            else if (varBins_all[j]>var_max_set) continue;
            else varBins.push_back(varBins_all[j]);
        }
        if (var_min_set < varBins[0]) varBins.insert(varBins.begin(),var_min_set);
        if (var_max_set > varBins[-1]) varBins.push_back(var_max_set);

    }
    
    TH1 *hData = BinnedHist->createHistogram(("RooDataHist_"+CatName).c_str(), *toyvar, Binning(b), Cut(("CMS_channel==CMS_channel::"+Cat).c_str()));
    TH1* hData_norm = (TH1*) hData->Clone(("hData_norm_"+CatName).c_str());

    limittree->GetEntry(i);
   

    hData->Draw();
    c->SaveAs("test.png");
    std::cout<<"r= "<<r<<std::endl;
}


}
 
