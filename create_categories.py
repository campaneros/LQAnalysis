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
#gROOT.LoadMacro(os.path.dirname(os.path.abspath(__file__))+"/../../src/libCpp/RooDoubleCBFast.cc+")
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
os.system(u"rm -f "+outputdirdatacards+"/*")
os.system("mkdir -p "+weboutputdir)
scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
os.system("cp "+scriptsPath+"/index.php "+weboutputdir)
filenameworkspace = "workspace.root"
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

        var_min_set = 450
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

        ## Get original TH1 histogram from root file
        rootfilename = inputdir+"/"+cat+"/"+filenameInput
        #print rootfilename    
        rootfile[icat] = ROOT.TFile.Open(rootfilename)
        print("Get "+histoname+" from file "+rootfilename)
        th1_fromFile[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
        th1_fromFile[icat].Draw()
        canvas[icat].SaveAs("tets.png")
        integral = int(th1_fromFile[icat].Integral())
        print("diocaro", integral)



        test=th1_fromFile[icat].ComputeIntegral()   
        print(test)
        ## Generate toy histogram
        gRandom = ROOT.TRandom()
        if(generateToy==1):
            th1_original[icat] = ROOT.TH1D("Toy","Toy", th1_fromFile[icat].GetNbinsX(), th1_fromFile[icat].GetXaxis().GetXmin(), th1_fromFile[icat].GetXaxis().GetXmax())
            #th1_original[icat] = TH1D("","", 5000, 0, 5000)
            gRandom.SetSeed(0)
            print("MAREMMA MAIALA")
            th1_original[icat].FillRandom(th1_fromFile[icat],integral)
            print("MAREMMA MAIALA")
            #th1_original[icat].Reset()
            #gRandom.SetSeed(0)
            #th1_original[icat].FillRandom(th1_fromFile[icat],integral)
            th1_original[icat].Write()
        else:
            th1_original[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
            # th1_original[icat].Add(th1_injected_2[icat], 1)
            # th1_original[icat].Add(th1_injected_1[icat], -2)

        var[icat] = ROOT.RooRealVar(varname+"_"+cat,vartitle,varBins_all[0],varBins_all[-1])
        var[icat].Print()

        print("------------------------------------------------------------")
        print("                   FINAL FIT                                ")
        print("------------------------------------------------------------")

        ## Main physics observable defined in fit range
        fitrange_L[0] = varBins[ 0]
        fitrange_R[0] = varBins[-1]

        var[icat] = ROOT.RooRealVar(varname+"_"+cat,vartitle,var_min_set,var_max_set)
        var[icat].Print()
        
        ## Create data histogram with coarser binning
        th1_rebin[icat] = th1_original[icat].Rebin(NvarBins,"th1_rebin_"+cat,array('d',varBins))
        print(th1_rebin[icat].GetBin(0))

        ## Create ROOT.RooDataHist in fit range from TH1
        rooHist[icat] = ROOT.RooDataHist("RooDataHist_"+cat,"RooDataHist_"+cat,ROOT.RooArgList(var[icat]),ROOT.RooFit.Import(th1_rebin[icat]))
        numberOfEvents[icat] = rooHist[icat].sumEntries()

        # Set RooPdf and parameters of selected background function, corresponding to fitFunction_name
        cwd_utils.set_bkg_fit_function(icat, cat, var, fitFunction_name, fitparam, nbkg, bkgPdf, ParametricBkgPdf, bkgExtPdf, th1_rebin, numberOfEvents)

        ## Bkg Fit
        ROOT.RooFit.SumW2Error(ROOT.kTRUE)
        fitResult[icat] = bkgExtPdf[icat].fitTo(rooHist[icat], 
                                                ROOT.RooFit.Strategy(1),
                                                ROOT.RooFit.Range(var_min_set,var_max_set),
                                                ROOT.RooFit.Save(ROOT.kTRUE),
                                                ROOT.RooFit.Verbose(ROOT.kFALSE),
                                                ROOT.RooFit.Extended(ROOT.kTRUE)
        )
        fitResult[icat].Print()
        Nevt[0] = int(nbkg[icat].getVal())
                    
        ## Fit quality
        print("Fit status: "+str(fitResult[icat].status()))
        print("Quality of covariance matrix: "+ str(fitResult[icat].covQual()))
        
        ## Create TF1 of background fit function    
        cwd_utils.bkgRooPdf_to_TF1(icat, cat, fitFunction_name, fitparam, bkgExtPdfTF1, nbkg, var_min_set, var_max_set)
        bkgExtPdfTF1[icat].SetLineColor(2)        
        
        ## Create background and pull histogram for final plot
        th1_rebin[icat].SetBinErrorOption(ROOT.TH1.kPoisson)
        th1_rebin_bkg[icat] = ROOT.TH1D("th1_rebin_bkg_"+cat,"th1_rebin_bkg_"+cat,NvarBins,array('d',varBins))
        th1_rebin_pull[icat] = ROOT.TH1D("th1_rebin_pull_"+cat,"th1_rebin_pull_"+cat,NvarBins,array('d',varBins))

        Npos[0] = 0
        Nneg[0] = 0
        WWruns[0]  = 1
        prev_pull_sign = 0
        pull_sign = 0
        for bin in range(NvarBins):
            bin+=1
            
            bin_low = th1_rebin[icat].GetBinLowEdge(bin)
            bin_up = th1_rebin[icat].GetBinLowEdge(bin)+th1_rebin[icat].GetBinWidth(bin)
        
            data = float(th1_rebin[icat].GetBinContent(bin))
            err_data_up = float(th1_rebin[icat].GetBinErrorUp(bin))
            
            exp_ = float(bkgExtPdfTF1[icat].Integral(bin_low,bin_up))

            # bkg histo
            th1_rebin_bkg[icat].SetBinContent(bin, exp_)

            # pull histo
            if data!=0:
                pull = (data - exp_) / err_data_up
                th1_rebin_pull[icat].SetBinContent(bin,pull)
                th1_rebin_pull[icat].SetBinError(bin,1)

                if pull >= 0:
                    Npos[0] += 1
                else:
                    Nneg[0] += 1
                if prev_pull_sign == 0:
                    prev_pull_sign = ROOT.TMath.Sign(1, pull)
                pull_sign = ROOT.TMath.Sign(1, pull)

            if pull_sign != prev_pull_sign and data != 0: 
                WWruns[0] += 1

        ## Chi2 and p-value
        icat_br[0] = icat
        cwd_utils.evaluate_chi2(icat, ndof, chi2, rchi2, ndof_allbins, chi2_allbins, reducedchi2_allbins, th1_rebin, th1_rebin_pull, len(fitparam)+1)
        pval[0] = ROOT.TMath.Prob(chi2[0], ndof[0])

        pval_ok[0] = 0

    ## Loop over signals
        signalInputfile = io.open(signalInputfilename, "r")
        for iline, line in enumerate(signalInputfile):
            line = line.rstrip('\n')
            splitline = line.split(" ")
            modell = splitline[0]
            category = splitline[1]
            M1 = splitline[2]
            L = splitline[3]
            Nsig = splitline[4]
            e_Nsig = splitline[5]
            mu = splitline[6]
            e_mu = splitline[7]
            ssigma = splitline[8]
            e_ssigma = splitline[9]
            a1 = splitline[10]
            e_a1 = splitline[11]
            n1 = splitline[12]
            e_n1 = splitline[13]
            a2 = splitline[14]
            e_a2 = splitline[15]
            n2 = splitline[16]
            e_n2 = splitline[17]
        	
        	#L=L.replace("p",".")
            print(modell)	
            if not category == cat:
                continue
            print("current signal = ", modell, category, M1, L, Nsig, e_Nsig, mu, e_mu, ssigma, e_ssigma, a1, e_a1, n1, e_n1, a2, e_a2, n2, e_n2)
            signalString = modell+"_"+category#+"_"+"M"+str(int(float(M1)))+"_R0p"+str(int(float(L)*10))+"_"+category
            print("signal_string   "+signalString)        	


        	## Signal pdf
            mean = ROOT.RooRealVar("mean_"+signalString,"mean_"+signalString,float(mu)) 
            width = ROOT.RooRealVar("width_"+signalString,"width_"+signalString,float(ssigma)) 
            alpha1 = ROOT.RooRealVar("alpha1_"+signalString,"alpha1_"+signalString,float(a1))
            n1 = ROOT.RooRealVar("n1_"+signalString,"n1_"+signalString,float(n1))
            alpha2 = ROOT.RooRealVar("alpha2_"+signalString,"alpha2_"+signalString,float(a2))
            n2 = ROOT.RooRealVar("n2_"+signalString,"n2_"+signalString,float(n2))
            nsig = ROOT.RooRealVar("ParametricSignalPdf_"+signalString+"_norm","ParametricSignalPdf_"+signalString+"_norm",float(Nsig),0,100000000)

            var_tmp = ROOT.RooRealVar("var_tmp","var_tmp",453, 5000)
            var_tmp.setRange("maximum_range", 453, 5000)  # create range to integrate over
            signalPdf_tmp = ROOT.RooDoubleCBFast("CB_tmp", "CB_tmp", var_tmp, mean, width, alpha1, n1, alpha2, n2)            
            intrinsicNorm = signalPdf_tmp.createIntegral(ROOT.RooArgSet(var_tmp), ROOT.RooFit.NormSet(ROOT.RooArgSet(var_tmp)), ROOT.RooFit.Range("maximum_range")) 
            var_tmp.setRange("fit_range", var_min_set, var_max_set)  # create range to integrate over
            signal_integral = signalPdf_tmp.createIntegral(ROOT.RooArgSet(var_tmp), ROOT.RooFit.NormSet(ROOT.RooArgSet(var_tmp)), ROOT.RooFit.Range("fit_range"))
            nsig.setVal( float(Nsig)*signal_integral.getVal()/intrinsicNorm.getVal() )
            print(str(nsig.getVal()))

        	
            signalPdf = ROOT.RooDoubleCBFast("signalPdf_"+signalString, "signalPdf_"+signalString, var[icat], mean, width, alpha1, n1, alpha2, n2)
        	
        #	meanShape_mu_err_cat = ROOT.RooRealVar("meanShape_mu_err_"+signalString, "meanShape_mu_err_"+signalString, 1, 0.5, 1.5)
        #	meanShape_ssigma_err_cat = ROOT.RooRealVar("meanShape_ssigma_err_"+signalString, "meanShape_ssigma_err_"+signalString, 1, 0.4, 1.6)
        #	mean_sist  = RooFormulaVar("mean_sist_" +signalString, "mean_sist_" +signalString, "mean_" +signalString+"*meanShape_mu_err_"+signalString   , ROOT.RooArgList(mean , meanShape_mu_err_cat   ))
       # 	width_sist = RooFormulaVar("width_sist_"+signalString, "width_sist_"+signalString, "width_"+signalString+"*meanShape_ssigma_err_"+signalString, ROOT.RooArgList(width, meanShape_ssigma_err_cat))
        	
        #	signalPdf = RooDoubleCBFast("signalPdf_"+signalString, "signalPdf_"+signalString, var[icat], mean_sist, width_sist, alpha1, n1, alpha2, n2)            
        	#ParametricSignalPdf_ = RooParametricShapeBinPdf("ParametricSignalPdf_"+signalString, "ParametricSignalPdf_"+signalString, signalPdf, var[icat], ROOT.RooArgList(mean_sist, width_sist, alpha1, n1, alpha2, n2), th1_rebin[icat])
            ParametricSignalPdf_ = ROOT.RooParametricShapeBinPdf("ParametricSignalPdf_"+signalString, "ParametricSignalPdf_"+signalString, signalPdf, var[icat], ROOT.RooArgList(mean, width, alpha1, n1, alpha2, n2), th1_rebin[icat])


            mean.setConstant(ROOT.kTRUE)
            width.setConstant(ROOT.kTRUE)
            alpha1.setConstant(ROOT.kTRUE)
            n1.setConstant(ROOT.kTRUE)
            alpha2.setConstant(ROOT.kTRUE)
            n2.setConstant(ROOT.kTRUE)
            nsig.setConstant(ROOT.kTRUE) #signal normalization should be set to constant (required by combine tool)

            signalPdf.Print()
            nsig.Print()

            ## Fill Workspace for signal
            # Import background pdf and normalization
            getattr(w,'import')(ParametricSignalPdf_)
                #getattr(w,'import')(signalPdf)
            getattr(w,'import')(nsig)

            os.system("mkdir -p "+outputdirdatacards+"/"+modell)    
            os.system("mkdir -p "+outputdirdatacards+"/"+modell+"/categories")    
        	    ## Create datacard for current signal
            datacardfilename = outputdirdatacards+"/"+modell+"/categories/"+"datacard_"+signalString+".txt"
            outputfile = io.open(datacardfilename,'w')

            ## Create datacard for current signal
            #datacardfilename = outputdirdatacards+"/"+"datacard_"+signalString+".txt"
            #outputfile = io.open(datacardfilename,'w')

            outputfile.write( u"## Datacard for "+signalString+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"imax * number of channels"+"\n" )
            outputfile.write( u"jmax * number of backgrounds"+"\n" )
            outputfile.write( u"kmax * number of nuisance parameters (source of systematic uncertainties)"+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"shapes sig "+cat+" "+filenameworkspacePath+" "+workspaceName+":"+ParametricSignalPdf_.GetName()+"\n" )
            outputfile.write( u"shapes bkg "+cat+" "+filenameworkspacePath+" "+workspaceName+":"+ParametricBkgPdf[icat].GetName()+"\n" )
            outputfile.write( u"shapes data_obs "+cat+" "+filenameworkspacePath+" "+workspaceName+":"+rooHist[icat].GetName()+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"bin \t\t"+cat+"\n" )
            outputfile.write( u"observation \t -1"+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"bin \t\t"+cat+"\t"+cat+"\n" )
            outputfile.write( u"process \t"+"sig \t\t"+ "bkg"+"\n" )
            outputfile.write( u"process \t"+ "0" +" \t\t"+"1"+"\n" )
            outputfile.write( u"rate \t\t"+ "1" +" \t\t"+"1"+"\n" )
            #outputfile.write( u"process \t"+"0 \t\t"+ "1"+"\n" )
            #outputfile.write( u"rate \t\t"+str(int(nsig.getValV()))+" \t\t"+ str(int(nbkg[icat].getValV())) +"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
        	#outputfile.write( p1[icat].GetName()+u" flatParam"+"\n" )
        	#outputfile.write( p2[icat].GetName()+u" flatParam"+"\n" )
        	#outputfile.write( p3[icat].GetName()+u" flatParam"+"\n" )
        	#outputfile.write( nbkg[icat].GetName()+u" flatParam"+"\n" )
		    #outputfile.write( u"lumi lnN 1.018 -\n")
            for param in fitparam:
                outputfile.write( param.GetName()+u" flatParam"+"\n" )
                outputfile.write( nbkg[icat].GetName()+u" flatParam"+"\n" )
            outputfile.close()            
        signalInputfile.close()

    ## Fill Workspace for data and background
    # Import data
        getattr(w,'import')(rooHist[icat])
   # Import background pdf and normalization
        getattr(w,'import')(ParametricBkgPdf[icat])
        getattr(w,'import')(nbkg[icat])

   ############### Plots ##################

   ## Create canvas
        canvas[icat] = ROOT.TCanvas("canvas_"+cat, "canvas_"+cat, 200, 10, 700, 500 )
        fPads1 = ROOT.TPad("pad1_"+cat, "pad1_"+cat, 0.00, 0.28, 0.99, 0.99)
        fPads2 = ROOT.TPad("pad2_"+cat, "pad2_"+cat, 0.00, 0.00, 0.99, 0.345)
   #
        fPads1.SetFillColor(0)
        fPads1.SetLineColor(0)
        fPads1.SetTopMargin(0.1)
        fPads1.SetBottomMargin(0.1)
        fPads1.SetRightMargin(0.1)
        fPads1.SetLeftMargin(0.13)
        fPads1.SetLogy()
   #
        fPads2.SetFillColor(0)
        fPads2.SetLineColor(0)
        fPads2.SetGridx()
        fPads2.SetGridy()
        fPads2.SetTopMargin(0.05)
        fPads2.SetBottomMargin(0.3)
        fPads2.SetRightMargin(0.1)
        fPads2.SetLeftMargin(0.13)
   #
        fPads1.Draw()
        fPads2.Draw()

   ## Pad1
        fPads1.cd()
        th1_rebin[icat].SetMinimum(0.03)
        th1_rebin[icat].GetYaxis().SetTitle("Number of events")
        th1_rebin[icat].GetYaxis().SetTitleSize(0.06)
        th1_rebin[icat].GetYaxis().SetTitleOffset(0.8)
        th1_rebin[icat].SetMarkerStyle(20)
        th1_rebin[icat].SetMarkerSize(0.8)
        th1_rebin[icat].SetMarkerColor(1)
        th1_rebin[icat].SetLineColor(1)
        th1_rebin[icat].SetStats(0)
        th1_rebin[icat].Draw("pe")
        th1_rebin_bkg[icat].SetLineColor(2)
        th1_rebin_bkg[icat].SetLineWidth(1)
        th1_rebin_bkg[icat].Draw("histsame")
        
   	#draw the lumi text on the canvas
        CMS_lumi.CMS_lumi(fPads1, iPeriod, iPos) 
        canvas[icat].Modified()
        canvas[icat].Update()
   #
        ## Pad2
        fPads2.cd()
        th1_rebin_pull[icat].GetXaxis().SetTitle(vartitle)
        th1_rebin_pull[icat].GetYaxis().SetTitle("#frac{Data - Fit}{Uncertainty} ")
        th1_rebin_pull[icat].SetTitle("")
        th1_rebin_pull[icat].SetMinimum(-4)    
        th1_rebin_pull[icat].SetMaximum(4)
        th1_rebin_pull[icat].SetLineColor(2)
        th1_rebin_pull[icat].SetFillColor(2)
        th1_rebin_pull[icat].SetMarkerStyle(20)
        th1_rebin_pull[icat].SetMarkerColor(1)
        th1_rebin_pull[icat].SetStats(0)
        th1_rebin_pull[icat].GetYaxis().SetNdivisions(405, ROOT.kTRUE)
        th1_rebin_pull[icat].GetXaxis().SetTitleSize(0.16)
        th1_rebin_pull[icat].GetXaxis().SetLabelSize(0.13)
        th1_rebin_pull[icat].GetXaxis().SetTitleOffset(0.83)
        th1_rebin_pull[icat].GetYaxis().SetTitleSize(0.12)
        th1_rebin_pull[icat].GetYaxis().SetLabelSize(0.11)
        th1_rebin_pull[icat].GetYaxis().SetTitleOffset(0.35)
        th1_rebin_pull[icat].Draw("hist")
   ## Pa
        ## Legend
        fPads1.cd()
        legend = ROOT.TLegend(0.7, 0.7, 0.87, 0.87)
        legend.SetLineColor(0)
        legend.SetLineWidth(0)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        if(generateToy==1):
            legend.AddEntry(th1_rebin[icat], "Toy Data", "p")
        else:
            legend.AddEntry(th1_rebin[icat], "Toy Data", "p")
        legend.AddEntry(th1_rebin_bkg[icat], "Fit", "l")
        legend.Draw()
        
       ## Plot fit results
        pt = ROOT.TPaveText(0.73, 0.58, 0.87, 0.87,"ndc")
        pt.SetFillColor(0)
        
        Chi2Text = "#chi^{2}="+str(round(chi2[0],2))+"   Ndf="+str(ndof[0])
        t1 = pt.AddText(Chi2Text)
        t1.SetTextColor(1)
        t1.SetTextSize( 0.04 )
        Chi2Text = "#chi^{2} / ndf (N_{bin}>10) = "+str(round(rchi2[0],2))+"  pval = "+str(round(pval[0],5))
        t1 = pt.AddText(Chi2Text)
        t1.SetTextColor(1)
        t1.SetTextSize( 0.04 )
        Chi2Text = "#chi^{2} / ndf (N_{bin}>10) = "+str(round(reducedchi2_allbins[0],2))+"  pval = "+str(round(pval[0],5))
        t1 = pt.AddText(Chi2Text)
        t1.SetTextColor(1)
        t1.SetTextSize( 0.04 )
        WWText = "N+ = "+str(int(Npos[0]))+"  N- = "+str(int(Nneg[0]))+"\n Nruns = "+str(int(WWruns[0])) + "  pval = "+str(round(WWpval[0],2)) 
        #t1 = pt.AddText(WWText)
        t1.SetTextColor(1)
        t1.SetTextSize( 0.04 )
        
        t_covqual = pt.AddText("covQual = "+str(round(fitResult[icat].covQual(),1)))    
        t_nbkg = pt.AddText("nbkg = "+str(round(nbkg[icat].getValV(),0))+" #pm "+str(round(nbkg[icat].getError(),0)))    
        
        for i,param in enumerate(fitparam):
            t_p = pt.AddText("p"+str(i+1)+" = "+str(round(param.getValV(),2))+" #pm "+str(round(param.getError(),2)))    
            t_p.SetTextColor(1)
            t_p.SetTextSize( 0.04 )
        t_covqual.SetTextColor(1)
        t_covqual.SetTextSize( 0.04 )
        t_nbkg.SetTextColor(1)
        t_nbkg.SetTextSize( 0.04 )
        
        pt.Draw("same")
        canvas[icat].Update()
	


        #pt = TPaveText(0.54, 0.58, 0.64, 0.87,"ndc")
        #pt.SetFillColor(0)
        
        #Chi2Text = "#chi^{2} / ndf (N_{bin}>10) = "+str(round(reducedChi2[icat],2))
        #t1 = pt.AddText(Chi2Text)
        #t1.SetTextColor(1)
        #t1.SetTextSize( 0.04 )
        
        #t_covqual = pt.AddText("covQual = "+str(round(fitResult[icat].covQual(),1)))    
        #t_nbkg = pt.AddText("nbkg = "+str(round(nbkg[icat].getValV(),0))+" #pm "+str(round(nbkg[icat].getError(),0)))    
        #t_p1 = pt.AddText("p1 = "+str(round(p1[icat].getValV(),2))+" #pm "+str(round(p1[icat].getError(),2)))    
        #t_p2 = pt.AddText("p2 = "+str(round(p2[icat].getValV(),2))+" #pm "+str(round(p2[icat].getError(),2)))    
        #t_p3 = pt.AddText("p3 = "+str(round(p3[icat].getValV(),2))+" #pm "+str(round(p3[icat].getError(),2)))    
        #t_covqual.SetTextColor(1)
        #t_covqual.SetTextSize( 0.04 )
        #t_nbkg.SetTextColor(1)
        #t_nbkg.SetTextSize( 0.04 )
        #t_p1.SetTextColor(1)
        #t_p1.SetTextSize( 0.04 )
        #t_p2.SetTextColor(1)
        #t_p2.SetTextSize( 0.04 )
        #t_p3.SetTextColor(1)
        #t_p3.SetTextSize( 0.04 )
        
        #pt.Draw("same")
        #canvas[icat].Update()
        
        ## Save canvas
        outputfilename_pdf = outputdir+"/"+"c_"+cat+".pdf"
        outputfilename_png = outputdir+"/"+"c_"+cat+".png"
        outputfilename_root = outputdir+"/"+"c_"+cat+".root"
        canvas[icat].SaveAs(outputfilename_pdf)
        canvas[icat].SaveAs(outputfilename_png)
        out=ROOT.TFile(outputfilename_root,"recreate")
        canvas[icat].Write()
        th1_original[icat].Write()
        th1_fromFile[icat].Write()
        os.system("cp "+ outputfilename_pdf + " " + weboutputdir)
        os.system("cp "+ outputfilename_png + " " + weboutputdir)
        
  ##
        rootfile[icat].Close()
        del fitparam[:]

w.writeToFile(filenameworkspacePath)
w.Delete()
print("Workspace at "+filenameworkspacePath)

#### Create final datacards (by combinging all categories for a given signal) ####

## Loop over signals (find the modells)
listOfSignalModels = []

signalInputfile = io.open(signalInputfilename, "r")
for iline, line in enumerate(signalInputfile):
    
    line = line.rstrip('\n')
    splitline = line.split(" ")
    
    modell = splitline[0]
    #category = splitline[1]
    M1 = splitline[2]
    L = splitline[3]

    signalStringModel = modell#+"_"+"M"+str(int(float(M1)))+"_R0p"+str(int(float(R)*10))

    if signalStringModel not in listOfSignalModels:
        listOfSignalModels.append(signalStringModel)

signalInputfile.close()
print(listOfSignalModels)

datacardList = io.open(outputdirdatacards+"/datacardList.txt", 'w')
datacardList_single = io.open(outputdirdatacards+"/datacardList_single_category.txt", 'w')
## Loop over signals and categories (create combined datacard)
for signal in listOfSignalModels:
    currentPath = os.getcwd()
    command = ""
    os.system("mkdir -p "+outputdirdatacards+"/"+signal)
    ## Loop over event categories    
    for icat, cat in enumerate(subDirList):
        signalStringCat = signal+"_"+cat
        datacardfilenameAll = outputdirdatacards+"/"+signal+"/"+"datacard_"+signal+".txt"
        datacardfilename = outputdirdatacards+"/"+signal+"/categories/"+"datacard_"+signalStringCat+".txt"        

        print(signal)
        if icat==0:
            command = "python "+combineCardScriptPath

        if not os.path.isfile(datacardfilename):
            continue        

        ## Select only mjet categories within -2ssigma and +1ssigma of the Res2 mass
        #resolution = 0.1
        #mjetMin = float(cat.split("_")[1])
        #mjetMax = float(cat.split("_")[2])
        #mjetmean = (mjetMax+mjetMin)/2 
        #mass = float(signal.split("_")[1].split("M")[1])
        #R = float((signal.split("_")[2].split("R0p")[1]))*0.1
        #mass2 = mass*R
        #ssigma_mass2 = mass2*resolution
        #print(mjetMin, mjetMax, mjetmean, mass, R, mass2, ssigma_mass2
        #if ( 
        #    ( (mass2-mjetMax)<=2*ssigma_mass2 and (mass2-mjetMax)>=0 ) 
        #    or ( (mjetMin-mass2)<=1*ssigma_mass2 and (mjetMin-mass2)>=0 )
        #    or ( mass2>=mjetMin and mass2<mjetMax )  
        #    ): 
        command += " "+cat+"="+datacardfilename
        print(signalStringCat)

        splitline = signalStringCat.split("_")
        modell = splitline[0]
        M1val = (splitline[1]).strip("M")
        Lval  = (splitline[2]).strip("L")
        Lval =  Lval.replace("p",".")
        category = splitline[3]         
        datacardList_single.write(modell+" "+M1val+" "+Lval+" "+category+" "+datacardfilename+"\    n")

        command += " > "+datacardfilenameAll
        print(command)
        os.system("cd "+outputdirdatacards+" ; "+command+" ; "+ "cd .." )
        print("Created final datacard at "+datacardfilenameAll) 
   
        splitline = signal.split("_")  
        modell = splitline[0]
        M1val = (splitline[1]).strip("M")
        val  = (splitline[2]).strip("L")
        val =  Lval.replace("p",".")

    datacardList.write(modell+" "+M1val+" "+Lval+" "+datacardfilenameAll+"\n")
    print("Created datacard List at ",outputdirdatacards,"/datacardList.txt")

#datacardList.close()


## Make multidimfit with fixed signal strength r=0

