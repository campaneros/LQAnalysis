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
//#include "RooBernstein.h"
#include "RooAbsPdf.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
//#include "RooPdf.h"
//#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include <time.h>
#include "TBox.h"
#include "TROOT.h"
#include "TRint.h"


using namespace RooFit;


//int AllFits(std::string filename,std::string path, double mass,double ctau,int channel, std::string wd,std::string sel,double lumiRatio, bool binned, bool blinded ){


void Make_x2(std::string ToysFile, std::string WorkspaceFile, std::string FitFile, int NToys, std::string CatName, std::string Cat, std::string Signal, std::string Output, std::string function, std::string expected_signal){

TFile *f = new TFile((WorkspaceFile).c_str());
TFile *toys= new TFile((ToysFile).c_str());
TFile *limitfile = new TFile((FitFile).c_str());

RooWorkspace *w = (RooWorkspace *)f->Get("w");
RooRealVar *x = w->var((CatName).c_str());

//gROOT->SetBatch(true);
//gErrorIgnoreLevel = kFatal;
//RooMsgService::instance().setGlobalKillBelow(WARNING);


//RooMsgService->instance()->setGlobalKillBelow(WARNING);


TTree* limittree = (TTree*) limitfile->Get("limit");
float_t r, p0_postfit, p1_postfit, p2_postfit, p3_postfit;
//limittree->SetBranchStatus("*",0);
//limittree->SetBranchStatus("r",1);
limittree->SetBranchAddress("r", &r);
limittree->SetBranchAddress(("trackedParam_shapeBkg_bkg_"+Cat+"__norm").c_str(), &p0_postfit);
limittree->SetBranchAddress(("trackedParam_p1_std"+Cat).c_str(), &p1_postfit);
limittree->SetBranchAddress(("trackedParam_p2_std"+Cat).c_str(), &p2_postfit);
limittree->SetBranchAddress(("trackedParam_p3_std"+Cat).c_str(), &p3_postfit);


int varBins_all [88] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058,1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250, 7500, 7750, 8000};


vector<double> varBins;
TCanvas* c= new TCanvas(("canvas_"+Cat).c_str(), ("canvas_"+Cat).c_str(), 200, 10, 700, 500);
TPad *fPads1 = new TPad(("pad1_"+Cat).c_str(), ("pad1_"+Cat).c_str(), 0.00, 0.28, 0.99, 0.99);
TPad *fPads2 = new TPad(("pad2_"+Cat).c_str(), ("pad2_"+Cat).c_str(), 0.00, 0.00, 0.99, 0.345);
RooAbsPdf *BkgFit = w->pdf(("ParametricBkgPdf_"+Cat).c_str());
RooRealVar *p1_var    = w->var(("p1_std"+Cat).c_str()); 
RooRealVar *p2_var    = w->var(("p2_std"+Cat).c_str()); 
RooRealVar *p3_var    = w->var(("p3_std"+Cat).c_str());
RooRealVar *norm_var  = w->var(("ParametricBkgPdf_"+Cat+"_norm").c_str());

RooExtendPdf ExtBkgFit = RooExtendPdf(("ExtBkgPdf_"+Cat).c_str(), ("ExtBkgPdf_"+Cat).c_str(), *BkgFit, *norm_var);

RooAbsPdf *signalPdf = w->pdf(("ParametricSignalPdf_"+Signal).c_str());
RooRealVar *nsig = w->var(("ParametricSignalPdf_"+Signal+"_norm").c_str());
RooExtendPdf signalExtendPdf = RooExtendPdf(("ParametricSignalPdf_"+Signal).c_str(), ("ParametricSignalPdf_"+Signal).c_str(), *signalPdf, *nsig);


double Ndof = -3;
double Chi2 = 0;
double GlobalNdof = - 3;
double CombineDof = - 3;
double redChi2 = 0.;
double pull = 0;
double GlobalChi2 = 0;
double CombineGoF = 0;

TH1 *hist_chi2= new TH1D("chi2","chi2",60,0,3);
TH1 *hist_pval= new TH1D("chi2","chi2",50,0,1);

for (int i=0; i<NToys; i++){
    CombineGoF=0;
    RooDataSet *dataset= (RooDataSet* ) toys->Get(("toys/toy_"+to_string(i+1)).c_str());
    RooDataHist *BinnedHist = dataset->binnedClone("BinnedClone","BinnedClone");
    const RooArgSet *argset = dataset->get();

    RooRealVar* toyvar = (RooRealVar* ) argset->find((CatName).c_str());
    double_t var_min_set = toyvar->getMin();
    double_t var_max_set = toyvar->getMax();


    RooAbsBinning& b = toyvar->getBinning();
    int NvarBins = toyvar->getBins();
    if (i==0){
        for( int j : varBins_all){
            //std::cout<<j<<std::endl;
            if (j<var_min_set) continue;
            else if (j>var_max_set) continue;
            else { varBins.push_back(j); 
            //std::cout<<j<<std::endl;
            }
        }
        if (var_min_set < varBins[0]) varBins.insert(varBins.begin(),var_min_set);
        if (var_max_set > varBins[-1]) varBins.push_back(var_max_set);
    }
    double arr_varBins[varBins.size()];
    //std::cout<<varBins.size()<<std::endl;
    std::copy(varBins.begin(),varBins.end(),arr_varBins);
    TH1 *hData = BinnedHist->createHistogram(("RooDataHist_"+CatName).c_str(), *toyvar, Binning(b), Cut(("CMS_channel==CMS_channel::"+Cat).c_str()));
    TH1 *hData_norm = (TH1*) hData->Clone(("hData_norm_"+CatName).c_str());

    TH1 *hSignal = new TH1D(("hSignal_"+Cat).c_str(),"",NvarBins, arr_varBins);
    TH1 *hBkg = new TH1D(("hBkg_"+Cat).c_str(),"",NvarBins, arr_varBins);
    TH1* hist_pull = new TH1D(("pull_"+Cat).c_str(), "", NvarBins, arr_varBins);
    TH1* hist_pull_signal = new TH1D(("pull_signal_"+Cat).c_str(), "", NvarBins, arr_varBins);

    x->setRange("global_range", arr_varBins[0], arr_varBins[NvarBins]);
    RooAbsReal *signal_pdfIntrinsicNorm = signalPdf->createIntegral(RooArgSet(*x),NormSet(RooArgSet(*x)),Range("global_range"));
    RooAbsReal *bkg_pdfIntrinsicNorm    = BkgFit->createIntegral(RooArgSet(*x),NormSet(RooArgSet(*x)),Range("global_range"));
       
    //std::cout<<signal_pdfIntrinsicNorm->getVal()<<std::endl;
    limittree->GetEntry(i);
    //std::cout<<"norm val "<<p0_postfit<<std::endl;
    norm_var->setVal(p0_postfit);
    p1_var->setVal(p1_postfit);
    p2_var->setVal(p2_postfit);
    p3_var->setVal(p3_postfit);


    for (int j=1; j<=NvarBins; j++){
        double bin_low   = hSignal->GetBinLowEdge(j);
        double bin_width = hSignal->GetBinWidth(j);
        double bin_up = bin_low + bin_width;
        x->setRange(("toy_"+to_string(j)).c_str(),bin_low,bin_up);
        RooAbsReal *signal_pdfIntegral = signalPdf->createIntegral(RooArgSet(*x),NormSet(RooArgSet(*x)),Range(("toy_"+to_string(j)).c_str()));
        double signal_pdfIntegral_norm = signal_pdfIntegral->getVal()/signal_pdfIntrinsicNorm->getVal();
        double bin_signalEvents = signal_pdfIntegral_norm*nsig->getVal()*r;
        double bin_signalEvents_norm = signal_pdfIntegral_norm*nsig->getVal()*r/bin_width;
        //std::cout<<bin_signalEvents_norm<<std::endl;

        hSignal->SetBinContent(j,bin_signalEvents_norm);

        RooAbsReal *bkg_pdfIntegral = BkgFit->createIntegral(RooArgSet(*x),NormSet(RooArgSet(*x)),Range(("toy_"+to_string(j)).c_str()));
        double bkg_pdfIntegral_norm = bkg_pdfIntegral->getVal()/bkg_pdfIntrinsicNorm->getVal();
        double bin_bkgEvents = bkg_pdfIntegral_norm*norm_var->getVal();
        double bin_bkgEvents_norm = bkg_pdfIntegral_norm*norm_var->getVal()/bin_width;
        
      
    
            
        hBkg->SetBinContent(j,bin_bkgEvents_norm);
        
        double data = hData->GetBinContent(j);
        hData_norm->SetBinContent(j, data/bin_width);
        hData_norm->SetBinError(j, TMath::Sqrt(data)/bin_width);

        double bkg_err = TMath::Sqrt(bin_bkgEvents_norm*bin_width)/bin_width;

        double bin_data = data/bin_width;


        if (bin_data!=0 || bin_bkgEvents_norm!=0){
            pull = (bin_data-bin_bkgEvents_norm)/bkg_err;
            //std::cout<<"pull "<<pull<<std::endl;
            //std::cout<<"data "<< bin_data <<"bkg"<< bin_bkgEvents_norm<<std::endl;
            hist_pull->SetBinContent(j,pull);
            CombineGoF += 2*(  (bin_bkgEvents_norm+bin_signalEvents_norm *r)-bin_data + bin_data*TMath::Log( bin_data/(bin_bkgEvents_norm+bin_signalEvents_norm*r) )  );
        } else  CombineGoF += 2*(  (bin_bkgEvents_norm+bin_signalEvents_norm *r)-bin_data);
        double pull_signal = 0;
        if (bkg_err != 0 ){
            pull_signal = (bin_signalEvents_norm)/bkg_err;
        }  else{
            pull_signal = 0; 
        }
        hist_pull_signal->SetBinContent(j,pull_signal);

        if ( bin_data>10/bin_width) {
            Chi2 += pull*pull;
            //std::cout<<Chi2<<std::endl;
            Ndof += 1.0;
            GlobalChi2 += pull*pull;
            GlobalNdof += 1.0;
        }
    }
    if (Ndof <  0) Ndof = 0;
    if (Ndof != 0) redChi2 = Chi2/Ndof;

    double pval= TMath::Prob(Chi2, Ndof);


    double frame_Ymax = 0;
    //#frame_Ymin =40
    //double frame_Ymin = hData->GetBinContent(hData->GetNbinsX()-1);
    //std::cout<<redChi2<<std::endl;
    hist_chi2->Fill(redChi2);
    hist_pval->Fill(pval);
    //hData_norm->Draw();
    //hBkg->Draw("SAMEL");
    //c->SaveAs("test.png");
    //std::cout<<"r= "<<r<<std::endl;
   // std::cout<<"norm val "<<p0_postfit<<std::endl;
    //delete[] arr_varBins;
    if (i%100==0){
        std::cout<<i<<std::endl;
    }



}

hist_chi2->Draw();
c->SaveAs((Output+"/"+function+"_"+Signal+"_expect_signal_"+expected_signal+"_toys"+to_string(NToys)+".png").c_str());
hist_pval->Draw();
c->SaveAs((Output+"/pval_"+function+"_"+Signal+"_expect_signal"+expected_signal+"_toys"+to_string(NToys)+".png").c_str());
limitfile->Close();
toys->Close();
f->Close();
}
 
