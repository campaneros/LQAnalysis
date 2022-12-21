#! /usr/bin/env python

from os import *
import os
import sys
import optparse
import datetime
import subprocess
import io

from array import array
from glob import glob
#from ROOT import *
import ROOT as ROOT


sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi
import create_workspaces_and_datacards_utils as cwd_utils



fitFunction_name = "std_4par"
fitparam = []
#ROOT.gROOT.LoadMacro(os.path.dirname(os.path.abspath(__file__))+"/../../src/libCpp/RooDoubleCBFast.cc+")
ROOT.gROOT.SetBatch(True)
gErrorIgnoreLevel = ROOT.kFatal
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

inputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/category"
filenameInput = "test_h1_mmuj_ak4.root"
subDirList = next(os.walk(inputdir))[1]
print(subDirList)

## Output directories  
outputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC"
outputdirdatacards = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards"
weboutputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_plot"
os.system("mkdir -p "+outputdir)
os.system("mkdir -p "+outputdirdatacards)
#os.system(u"rm -f "+outputdirdatacards+"/*")
os.system("mkdir -p "+weboutputdir)
scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
os.system("cp "+scriptsPath+"/index.php "+weboutputdir)
filenameworkspace = "workspace_gen.root"
filenameworkspacePath = outputdir+"/"+filenameworkspace
workspaceName = "w"








## Histogram
#if opt.fitData:
#    generateToy = -1
#    histoname = "h1_mmuj_ak4__DATA"
#else:
generateToy = 1
histoname = "htot_bkg_MC_m_muj_ak4"

statMultiplier = 1
varname = "m_muj_ak4"
vartitle = "m_{\muj}_ak4 [GeV]"

## Signal input
signalInputfilename = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output/signals.txt"


ncategories = len(subDirList)



## Binning
varBins_all = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250, 7500, 7750, 8000]
NvarBins_all = len(varBins_all)-1

## CMS_lumi variables (see CMS_lumi.py)
lumi = 140
CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(lumi)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

##Combine commands: 
combineCardScriptPath = os.path.dirname(os.path.abspath(__file__))+"/../"+"combineCards_modified.py" #allows absolute paths of files

#opt.
fitRanges =0 

### ML Fit with combine ###
#Data no signal
#combine -M MultiDimFit datacards/datacard_gg_M3000_R0p1.txt --rMin -5 --rMax 5 -S 0 --algo singles --cl=0.68 --saveFitResult --trackParameters "rgx{.*}" -t 0 --name _santanas_ --keepFailures (--robustHesse 1)

#Toy no signal
#combine -M MultiDimFit datacards/datacard_gg_M3000_R0p1.txt --rMin -5 --rMax 5 -S 0 --algo singles --cl=0.68 --saveFitResult --saveToys --trackParameters "rgx{.*}" -t 10 --name _santanas_ --keepFailures -s 123456 (--robustHesse 1)

#Toy with signal
#combine -M MultiDimFit datacards/datacard_gg_M3000_R0p1.txt --rMin -5 --rMax 5 -S 0 --algo singles --cl=0.68 --saveFitResult --saveToys --trackParameters "rgx{.*}" -t 10 --name _santanas_ --keepFailures --expectSignal 1 -s 123456 (--robustHesse 1)

#Access fit results (higgsCombineTest.MultiDimFit.mH120.123456.root)
#limit->Scan("iToy:quantileExpected:r:covQual:trackedParam_r:trackedParamErr_r","quantileExpected==-1")

### Old ###
#combine -M FitDiagnostics /afs/cern.ch/work/s/santanas/Releases/CMSSW_9_4_0_trijet/src/CMSJET/RootTreeAnalyzer/output/prod1_26_09_18/fitMjetCategories4/datacards/datacard_gg_M3000_R0p1.txt --verbose 3 -S 0 --rMin -20 --rMax 20 (--robustHesse 1)
#Signal fits

##############################################################################################

# Definition of arrays
rootfile = [None] * ncategories

var = [None] * ncategories

th1_fromFile = [None] * ncategories
th1_original = [None] * ncategories
th1_rebin = [None] * ncategories
th1_rebin_bkg = [None] * ncategories
th1_rebin_pull = [None] * ncategories
rooHist = [None] * ncategories
rooHist_bkg = [None] * ncategories
numberOfEvents = [None] * ncategories

nbkg = [None] * ncategories
p1 = [None] * ncategories
p2 = [None] * ncategories
p3 = [None] * ncategories
bkgPdf = [None] * ncategories
bkgExtPdf = [None] * ncategories
bkgExtPdfTF1 = [None] * ncategories
fitResult = [None] * ncategories
ndof = [None] * ncategories
Chi2 = [None] * ncategories
reducedChi2 = [None] * ncategories
th1_error = [None]*ncategories


ParametricBkgPdf = [None] * ncategories



canvas = [None] * ncategories

## Create Workspace
w = ROOT.RooWorkspace(workspaceName,workspaceName)
    # output tree with chi2
M1_br = array('f', [0])
R_br  = array('f', [0])
icat_br = array('i', [0])
Nevt  = array('i', [0])

ndof = array('i', [0])
chi2 = array('f', [0.])
rchi2 = array('f', [0.])
pval = array('f', [0.])
orig_pval = array('f', [0.])
pval_ok = array('i', [0])

ndof_allbins = array('i', [0])
chi2_allbins = array('f', [0.])
reducedchi2_allbins = array('f', [0.])

Npos   = array('f', [0.])
Nneg   = array('f', [0.])
WWruns = array('f', [0.])
WWpval = array('f', [0.])

fitrange_L = array('f', [0.])
fitrange_R = array('f', [0.])


## Loop over signal categories  
#var_min_limit = M1val*0.85
for icat, cat in enumerate(subDirList):
	#if ("2Muon" in str(cat)):	
	#	continue
        print("\n")
        print("#######################################################")
        print("Cat: "+str(icat))
        print("#######################################################")
    
        counter = 0
        #category_edges = cat.split("_")
        #P3mjet_low = float(category_edges[1])
        #P3mjet_high = float(category_edges[2])
        #Res2mjet_low = float(category_edges[4])
        #Res2mjet_high = float(category_edges[5])

        var_min_set = 453
        var_max_set = 5000
        
        ## Modify variable binning
        i_min_border = 0
        i_max_border = 0
        varBins = []
        for iborder, border in enumerate(varBins_all):
            if border<var_min_set:
                i_min_border += 1
                continue
            if border>var_max_set:
                continue
            else:
                varBins.append(float(border))
                i_max_border = iborder

        if fitRanges:
            print(var_min_set, var_max_set)
            if var_min_set < varBins[0] and abs(var_min_set - varBins[0])>0.001:
                varBins.insert(0, var_min_set)
            else:
                var_min_set = varBins[0]
            if var_max_set > varBins[-1]:
                varBins.append(var_max_set)
        else:
            var_min_set = varBins[0 ] 
            var_max_set = varBins[-1]
        print("var_min_set: ", str(var_min_set), "varBins[0]: ", str(varBins[0]), "varBins[1]: ", str(varBins[1]))

        NvarBins = len(varBins)-1
        canvas[icat] = ROOT.TCanvas("canvas_"+cat, "canvas_"+cat, 200, 10, 700, 500 )


        #print(varBins)
        #prova=array('d',varBins)
        #print(prova)	
        ## Get original TH1 histogram from root file
        rootfilename = inputdir+"/"+cat+"/"+filenameInput
        #print rootfilename    
        rootfile[icat] = ROOT.TFile.Open(rootfilename)
        print("Get "+histoname+" from file "+rootfilename)
        th1_fromFile[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
        th1_fromFile[icat].Draw()
        #canvas[icat].SaveAs("tets.png")
        integral = int(th1_fromFile[icat].Integral())
        print("diocaro", integral)

	#outRoot=ROOT.TFile( outputdir+"/background"+cat+".root","recreate")
	#th1_fromFile[icat].Write()
	#getattr(w,'import')(th1_fromFile[icat])
	#print("dai")
	#print("cazo")
	#clone_histo.Write()

        

      	  
        ## Create ROOT.RooDataHist in fit range from TH1
        #rooHist[icat] = ROOT.RooDataHist("ROOT.RooDataHist_"+cat,"ROOT.RooDataHist_"+cat,ROOT.RooArgList(var[icat]),ROOT.RooFit.Import(th1_rebin[icat]))
        #numberOfEvents[icat] = rooHist[icat].sumEntries()
  
        test=th1_fromFile[icat].ComputeIntegral() 
        print(test)
        ## Generate toy histogram
        gRandom = ROOT.TRandom()
        if(generateToy==1):
            th1_original[icat] = ROOT.TH1D("Toy","Toy", th1_fromFile[icat].GetNbinsX(), th1_fromFile[icat].GetXaxis().GetXmin(), th1_fromFile[icat].GetXaxis().GetXmax())
            #th1_original[icat] = ROOT.TH1D("","", 5000, 0, 5000)
            gRandom.SetSeed(0)
            print("MAREMMA MAIALA")
            th1_original[icat].FillRandom(th1_fromFile[icat],integral)
            print("MAREMMA MAIALA")
            #th1_original[icat].Reset()
            #gRandom.SetSeed(0)
            #th1_original[icat].FillRandom(th1_fromFile[icat],integral)
	    #th1_original[icat].Write()
        else:
            th1_original[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
            # th1_original[icat].Add(th1_injected_2[icat], 1)
            # th1_original[icat].Add(th1_injected_1[icat], -2)

        var[icat] = ROOT.RooRealVar(varname+"_"+cat,vartitle,varBins_all[0],varBins_all[-1])
        var[icat].Print()

        print("------------------------------------------------------------")
        print("                   FINAL FIT                                ")
        print("------------------------------------------------------------")

        var[icat] = ROOT.RooRealVar(varname+"_"+cat,vartitle,var_min_set,var_max_set)
        var[icat].Print()

        ROOT.RooFit.SumW2Error(ROOT.kTRUE)
        ## Create data histogram with coarser binning
        th1_rebin_bkg[icat] = th1_original[icat].Rebin(NvarBins,"th1_rebin_"+cat,array('d',varBins))
        th1_rebin[icat] = th1_fromFile[icat].Rebin(NvarBins,"th1_rebin_"+cat,array('d',varBins))
        ## Create ROOT.RooDataHist in fit range from TH1
        rooHist[icat] = ROOT.RooDataHist("ROOT.RooDataHist_"+cat,"ROOT.RooDataHist_"+cat,ROOT.RooArgList(var[icat]),ROOT.RooFit.Import(th1_rebin[icat]))
        rooHist_bkg[icat] = ROOT.RooDataHist("ROOT.RooDataHist_bkg_"+cat,"ROOT.RooDataHist_bkg_"+cat,ROOT.RooArgList(var[icat]),ROOT.RooFit.Import(th1_rebin_bkg[icat]))
        numberOfEvents[icat] = rooHist_bkg[icat].sumEntries()
        nbkg[icat] = rooHist[icat].sumEntries()
        th1_error[icat] = ROOT.TH1D("Error"+cat,"Error"+cat, NvarBins,array('d',varBins))
        #th1_rebin[icat].SetBinErrorOption(ROOT.TH1.kPoisson)
        #th1_rebin_bkg[icat].SetBinErrorOption(ROOT.TH1.kPoisson)
        for ibin,bin in enumerate(varBins):
            ibin+=1
            bincontent = th1_rebin[icat].GetBinContent(ibin)
            binerror = th1_rebin[icat].GetBinError(ibin)
            error = ROOT.TMath.Sqrt(bincontent)
            bin_low   = th1_rebin[icat].GetBinLowEdge(ibin)
            bin_width = th1_rebin[icat].GetBinWidth(ibin)
            print("Bin low ", bin_low, "  Bin Up: ", bin_low+bin_width, " bin content: ", bincontent)
            print("SumW2: ", binerror, " Sqrt(n):", error)
            if error != 0:
                th1_error[icat].SetBinContent(ibin,binerror/error)
            else:
                th1_error[icat].SetBinContent(ibin,0)
            th1_error[icat].SetBinError(ibin,0)
        canvas[icat].SetLogy(1)
        th1_error[icat].SetMinimum(0.01)
        th1_error[icat].GetYaxis().SetTitle("Mass")
        th1_error[icat].GetYaxis().SetTitle("getbinerror/sqrt(n)")
        th1_error[icat].GetYaxis().SetTitleSize(0.06)
        th1_error[icat].GetYaxis().SetTitleOffset(0.8)
        th1_error[icat].SetMarkerStyle(20)
        th1_error[icat].SetMarkerSize(0.8)
        th1_error[icat].SetMarkerColor(1)
        th1_error[icat].SetLineColor(1)
        th1_error[icat].Draw("P")
        for ext in ['.png','.pdf']:
            canvas[icat].SaveAs("Error_histo"+cat+ext)

            



