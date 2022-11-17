#! /usr/bin/env python

from os import *
import os
import sys
import optparse
import datetime
import subprocess
import io
from collections import OrderedDict

import numpy as np

from array import array
from glob import glob
from ROOT import *
import ROOT

import create_workspaces_and_datacards_utils as cwd_utils
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi

gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

inputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/category"
filenameInput = "histograms.root"
subDirList = next(os.walk(inputdir))[1]
print(subDirList)

## Output directories  
outputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/output_MC/"
outputdirdatacards = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/datacards"
weboutputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/output_plot"
os.system("mkdir -p "+outputdir)
os.system("mkdir -p "+outputdirdatacards)
os.system(u"rm -f "+outputdirdatacards+"/*")
os.system("mkdir -p "+weboutputdir)
scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
os.system("cp "+scriptsPath+"/index.php "+weboutputdir)
filenameworkspace = "workspace.root"
filenameworkspacePath = outputdir+"/"+filenameworkspace
workspaceName = "w"



##################################################################################################

## Input directories (each folder contains a root file with histograms)
inputdir = opt.inputdir
filenameInput = "test_h1_mmuj_ak4.root"
subDirList = next(os.walk(inputdir))[1]
#print(subDirList

#fit function name
fitFunction_name = opt.fitFunction

## Output directories  
outputdir = opt.outputdir
outputdirplots = outputdir+"/plots_"+fitFunction_name
outputdirdatacards = outputdir+"/datacards_"+fitFunction_name
weboutputdir = opt.weboutputdir+"/plots_"+fitFunction_name+""
os.system("mkdir -p "+outputdir)
os.system("mkdir -p "+outputdirdatacards)
#os.system(u"rm -r "+outputdirdatacards+"/*")
#os.system(u"rm -r "+outputdirdatacards+"/../*")
os.system("mkdir -p "+weboutputdir)
scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
os.system("cp "+scriptsPath+"/index.php "+weboutputdir)
os.system("cp "+scriptsPath+"/index.php "+weboutputdir+"/../")
os.system("cp "+scriptsPath+"/index.php "+weboutputdir+"/../../")

## Histogram
if opt.fitData:
    generateToy = -1
    histoname = "h1_mmuj_ak4__DATA"
else:
    generateToy = 1
    histoname = "htot_bkg_MC_m_muj_ak4"

statMultiplier = 1
varname = "m_muj_ak4"
vartitle = "m_{\muj} [GeV]"

## Signal Input
#signalInputdirname  = opt.signaldir
#signalInputfilename = opt.signalfile
signalInputfilename = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/output/signals.txt"

##############################################################################################
## Parse systematics files
syst_dir   = opt.syst_dir
syst_22cat = OrderedDict()
syst_9cat  = OrderedDict()
syst_1cat  = OrderedDict()

cwd_utils.parse_syst_file( syst_22cat, syst_9cat, syst_1cat, syst_dir )


## Binning
varBins_all = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7150, 7500, 7850, 8250, 8650, 8999, 9500, 9999]
#varBins_all = range(0,8000,150)
NvarBins_all = len(varBins_all)-1

## CMS_lumi variables (see CMS_lumi.py)
#lumi = 4.0
lumi = 137.098
CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(lumi)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

##Combine commands: 
combineCardScriptPath = os.path.dirname(os.path.abspath(__file__))+"/../../"+"combineCards_modified.py" #allows absolute paths of files

##############################################################################################
# Load fit ranges from file 

if opt.fitRanges:
    fitRangesFile = TFile.Open(opt.fitRanges)
    fitRangesTree = fitRangesFile.Get("chi2tree")

##############################################################################################

chi2file = [None] * len(subDirList)
chi2tree = [None] * len(subDirList)

efftreefile = TFile.Open(outputdir+"/final_efficiencies.root", "RECREATE")
efftreefile.cd()
efftree = TTree("efftree", "efftree")

icat_i = array("i", [1 ])
M1_i   = array("d", [0.])
R_i    = array("d", [0.])
nsig_i = array("d", [0.])
nerr_i = array("d", [0.])
Ntot_i = array("d", [0.])
Seff_i = array("d", [0.])
Eeff_i = array("d", [0.])

efftree.Branch("icat", icat_i, "icat/I")
efftree.Branch("M1"  , M1_i  , "M1/D"  )
efftree.Branch("R"   , R_i   , "R/D"   )
efftree.Branch("nsig", nsig_i, "nsig/D")
efftree.Branch("nerr", nerr_i, "nerr/D")
efftree.Branch("Ntot", Ntot_i, "Ntot/D")
efftree.Branch("Seff", Seff_i, "Seff/D")
efftree.Branch("Eeff", Eeff_i, "Eeff/D")

## Loop over signals
for isignal, signalname in enumerate(subDirList):
    ## Get signal model, mass, R
    splitsig = signalname.split("_")
    model    = splitsig[0]
    M1val    = float((splitsig[1]).strip("M"))
    #Lval     = float((splitsig[2]).strip("R0p").strip("\n"))
    #Lval = 
    #print(Lval)



    #print(M1val, Lval)

    # if M1val != 8700:
    #     continue
    # if int(M1val / 100) % 5 != 0:
    #     continue

    ## Create Workspace
    filenameworkspace = "workspace_"+signalname+".root"
    filenameworkspacePath = str(os.environ['CMSSW_BASE'])+"/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/"+outputdirdatacards+"/"+signalname+"/"+filenameworkspace
    workspaceName = "w_"+signalname

    w = RooWorkspace(workspaceName,workspaceName)

    # output sub_directory for each signal
    sub_outputdir = outputdirplots+"/"+signalname
    sub_outputdirdatacards = outputdirdatacards+"/"+signalname
    sub_weboutputdir = weboutputdir+"/"+signalname
    os.system("mkdir -p "+sub_outputdir)
    os.system("mkdir -p "+sub_outputdirdatacards)
    #os.system(u"rm -f "+sub_outputdirdatacards+"/*")
    os.system("mkdir -p "+sub_weboutputdir)
    os.system("cp "+scriptsPath+"/index.php "+sub_weboutputdir)
    
    # Select signal and systematics
    categoriesList = next(os.walk(inputdir+"/"+signalname))[1]
    ncategories = len(categoriesList)
    catType = 1
    signalInputfilename = signalInputdirname+"/signals.txt"
    syst_dict = syst_22cat
    # if M2val <= 600:
    #     catType = 1
    #     signalInputfilename = signalInputdirname+"/interpolated_signals_22cat.txt"
    #     syst_dict = syst_22cat
    # elif M2val <= 1200:
    #     catType = 2
    #     signalInputfilename = signalInputdirname+"/interpolated_signals_9cat.txt"
    #     syst_dict = syst_9cat
    # else:
    #     catType = 4
    #     signalInputfilename = signalInputdirname+"/interpolated_signals_1cat.txt"
    #     syst_dict = syst_1cat

    # Definition of arrays
    rootfile = [None] * ncategories

    var = [None] * ncategories
    
    th1_fromFile = [None] * ncategories
    # th1_injected_1 = [None] * ncategories
    # th1_injected_2 = [None] * ncategories
    th1_original = [None] * ncategories
    th1_rebin = [None] * ncategories
    th1_rebin_bkg = [None] * ncategories
    th1_rebin_pull = [None] * ncategories
    rooHist = [None] * ncategories
    numberOfEvents = [None] * ncategories
    
    nbkg = [None] * ncategories
    fitparam = [] 
    bkgPdf = [None] * ncategories
    bkgExtPdf = [None] * ncategories
    bkgExtPdfTF1 = [None] * ncategories
    ParametricBkgPdf = [None] * ncategories
    bkgExtPdfTF1 = [None] * ncategories
    fitResult = [None] * ncategories

    canvas = [None] * ncategories

    # output tree with chi2
    M1_br = array('f', [M1val])
    R_br  = array('f', [round(Rval,4)])
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

    print("Create Chi2file at "+sub_outputdir+"/chi2_"+signalname+".root")
    chi2file[isignal] = TFile(sub_outputdir+"/chi2_"+signalname+".root","RECREATE")
    chi2file[isignal].cd()
    chi2tree[isignal] = TTree("chi2tree", "chi2tree")

    chi2tree[isignal].Branch("M1 " , M1_br  , "M1/F  ")
    chi2tree[isignal].Branch("rho" , R_br   , "Rho/F ")
    chi2tree[isignal].Branch("iCat", icat_br, "iCat/I")
    chi2tree[isignal].Branch("Nevt", Nevt   , "Nevt/I")
    
    chi2tree[isignal].Branch("Ndof", ndof, "Ndof/I")
    chi2tree[isignal].Branch("Chi2", chi2, "Chi2/F")
    chi2tree[isignal].Branch("ReducedChi2", rchi2, "ReducedChi2/F")
    chi2tree[isignal].Branch("Pval", pval, "Pval/F")
    chi2tree[isignal].Branch("Original_Pval", orig_pval, "Original_Pval/F")
    chi2tree[isignal].Branch("Pval_ok", pval_ok, "Pval_ok/I")

    chi2tree[isignal].Branch("Ndof_allbins", ndof_allbins, "Ndof_allbins/I")
    chi2tree[isignal].Branch("Chi2_allbins", chi2_allbins, "Chi2_allbins/F")
    chi2tree[isignal].Branch("ReducedChi2_allbins", reducedchi2_allbins, "ReducedChi2_allbins/F")

    chi2tree[isignal].Branch("Npos", Npos, "Npos/F")
    chi2tree[isignal].Branch("Nneg", Nneg, "Nneg/F")
    chi2tree[isignal].Branch("WWruns", WWruns, "Wald-Wolfowitz test WWruns/F")
    chi2tree[isignal].Branch("WWpval", WWpval, "Wald-Wolfowitz test WWpval/F")

    chi2tree[isignal].Branch("Fitrange_L", fitrange_L, "Fitrange_L/F")
    chi2tree[isignal].Branch("Fitrange_R", fitrange_R, "Fitrange_R/F")

    ## Loop over signal categories  
    #var_min_limit = M1val*0.85
    for icat, cat in enumerate(categoriesList):

        print("\n")
        print("#######################################################")
        print(signalname, "Cat: "+str(icat))
        print("#######################################################")
    
        counter = 0
        category_edges = cat.split("_")
        P3mjet_low = float(category_edges[1])
        P3mjet_high = float(category_edges[2])
        Res2mjet_low = float(category_edges[4])
        Res2mjet_high = float(category_edges[5])

        if opt.fitRanges:
            if Rval < 0.1:
                selection = "M1==%s && (rho-%f<0.001 & rho-%f>-0.001) && iCat==%s" % (str(int(M1val)), 0.1, 0.1, str(icat))
            else:
                selection = "M1==%s && (rho-%f<0.001 & rho-%f>-0.001) && iCat==%s" % (str(int(M1val)), Rval, Rval, str(icat))
            print(selection
            npoints = fitRangesTree.Draw("Fitrange_L",selection,"goff")
            fitRange_L_fromFile = np.ndarray(npoints, 'd', fitRangesTree.GetV1())
            #var_min_set = max( 1600, min(8500*0.9, fitRange_L_fromFile[0]*0.9) )# nominal fit ranges -10%
            var_min_set = min(8500    , fitRange_L_fromFile[0]    )
            #var_min_set = min(8500*1.1, fitRange_L_fromFile[0]*1.1) # nominal fit ranges +10%
        else:
            var_min_set = max(M1val*0.85, 1800)

        if var_min_set < M1val:
            var_max_set = M1val*1.25
        else:
            var_max_set = var_min_set*1.4
        
        #var_max_set = 8000
        
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

        if opt.fitRanges:
            print(var_min_set, var_max_set
            if var_min_set < varBins[0] and abs(var_min_set - varBins[0])>0.001:
                varBins.insert(0, var_min_set)
            else:
                var_min_set = varBins[0]
            if var_max_set > varBins[-1]:
                varBins.append(var_max_set)
        else:
            var_min_set = varBins[0 ] 
            var_max_set = varBins[-1]
        print("var_min_set: ", str(var_min_set), "varBins[0]: ", str(varBins[0]), "varBins[1]: ", str(varBins[1])

        NvarBins = len(varBins)-1

        ## Get original TH1 histogram from root file
        rootfilename = inputdir+"/"+signalname+"/"+cat+"/"+filenameInput
        rootfile[icat] = TFile.Open(rootfilename)
        rootfile[icat].cd()
        th1_fromFile[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
        # th1_injected_1[icat] = rootfile[icat].Get("h1_mjj__INJECTED_S").Clone("h_injected_1")
        # th1_injected_1[icat] = gROOT.FindObject("h1_mjj__INJECTED_S;1")
        # th1_injected_2[icat] = gROOT.FindObject("h1_mjj__INJECTED_S;1")
        # th1_fromFile[icat].Add(th1_injected_2[icat], -1)
        # th1_fromFile[icat].Add(th1_injected_1[icat], -1)
        integral = int(th1_fromFile[icat].Integral()*statMultiplier)
                
        ## Generate toy histogram
        gRandom = TRandom()
        if(generateToy==1):
            th1_original[icat] = TH1F("","", th1_fromFile[icat].GetNbinsX(), th1_fromFile[icat].GetXaxis().GetXmin(), th1_fromFile[icat].GetXaxis().GetXmax())
            gRandom.SetSeed(0)
            th1_original[icat].FillRandom(th1_fromFile[icat],integral)
            #th1_original[icat].Reset()
            #gRandom.SetSeed(0)
            #th1_original[icat].FillRandom(th1_fromFile[icat],integral)
        else:
            th1_original[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
            # th1_original[icat].Add(th1_injected_2[icat], 1)
            # th1_original[icat].Add(th1_injected_1[icat], -2)

        var[icat] = RooRealVar(varname+"_"+cat,vartitle,varBins_all[0],varBins_all[-1])
        var[icat].Print()

        print("------------------------------------------------------------"
        print("                   FINAL FIT                                "
        print("------------------------------------------------------------"

        ## Main physics observable defined in fit range
        fitrange_L[0] = varBins[ 0]
        fitrange_R[0] = varBins[-1]

        var[icat] = RooRealVar(varname+"_"+cat,vartitle,var_min_set,var_max_set)
        var[icat].Print()
        
        ## Create data histogram with coarser binning
        th1_rebin[icat] = th1_original[icat].Rebin(NvarBins,"th1_rebin_"+cat,array('d',varBins))

        ## Create RooDataHist in fit range from TH1
        rooHist[icat] = RooDataHist("RooDataHist_"+cat,"RooDataHist_"+cat,RooArgList(var[icat]),RooFit.Import(th1_rebin[icat]))
        numberOfEvents[icat] = rooHist[icat].sumEntries()

        # Set RooPdf and parameters of selected background function, corresponding to fitFunction_name
        cwd_utils.set_bkg_fit_function(icat, cat, var, fitFunction_name, fitparam, nbkg, bkgPdf, ParametricBkgPdf, bkgExtPdf, th1_rebin, numberOfEvents)

        ## Bkg Fit
        RooFit.SumW2Error(kTRUE)
        fitResult[icat] = bkgExtPdf[icat].fitTo(rooHist[icat], 
                                                RooFit.Strategy(1),
                                                RooFit.Range(var_min_set,var_max_set),
                                                RooFit.Save(kTRUE),
                                                RooFit.Verbose(kFALSE),
                                                RooFit.Extended(kTRUE)
        )
        fitResult[icat].Print()
        Nevt[0] = int(nbkg[icat].getVal())
                    
        ## Fit quality
        print("Fit status: "+str(fitResult[icat].status()) 
        print("Quality of covariance matrix: "+ str(fitResult[icat].covQual())
        
        ## Create TF1 of background fit function    
        cwd_utils.bkgRooPdf_to_TF1(icat, cat, fitFunction_name, fitparam, bkgExtPdfTF1, nbkg, var_min_set, var_max_set)
        bkgExtPdfTF1[icat].SetLineColor(2)        
        
        ## Create background and pull histogram for final plot
        th1_rebin[icat].SetBinErrorOption(TH1.kPoisson)
        th1_rebin_bkg[icat] = TH1F("th1_rebin_bkg_"+cat,"th1_rebin_bkg_"+cat,NvarBins,array('d',varBins))
        th1_rebin_pull[icat] = TH1F("th1_rebin_pull_"+cat,"th1_rebin_pull_"+cat,NvarBins,array('d',varBins))

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
                    prev_pull_sign = TMath.Sign(1, pull)
                pull_sign = TMath.Sign(1, pull)

            if pull_sign != prev_pull_sign and data != 0: 
                WWruns[0] += 1

        ## Chi2 and p-value
        icat_br[0] = icat
        cwd_utils.evaluate_chi2(icat, ndof, chi2, rchi2, ndof_allbins, chi2_allbins, reducedchi2_allbins, th1_rebin, th1_rebin_pull, len(fitparam)+1)
        pval[0] = TMath.Prob(chi2[0], ndof[0])

        pval_ok[0] = 0
        if var_min_set < 0.9*M1val:
            pval_ok[0] = 1

        chi2tree[isignal].Fill()

        # Wald-Wolfowitz runs test
        # Npulls   = float( Npos[0] + Nneg[0] )
        
        # WWmean  = float( 2*Npos[0]*Nneg[0] / Npulls +1 )
        # WWsigma = float( TMath.Sqrt( (WWmean-1)*(WWmean-2)/(Npulls-1) ) )
        # WWnorm  = 1 #/ ( WWsigma * TMath.Sqrt(2*TMath.Pi()) )
            
        # gaussian = TF1("WWpdf", "gaus", WWmean-10*WWsigma, WWmean+10*WWsigma)
        # gaussian.SetParameters(WWnorm, WWmean, WWsigma)
            
        # WWpval[0] = gaussian.Integral(WWruns[0], WWmean+10*WWsigma)


        ##############################################################
        ## Loop over signals                                        ##
        ##############################################################
        print(signalname
        print(signalInputfilename
        signalInputfile = io.open(signalInputfilename, "r")

        for iline, line in enumerate(signalInputfile):
            if not pval_ok[0]:
                break

            line = line.rstrip('\n')
            splitline = line.split(" ")
            
            modillo = splitline[0]
            category = splitline[1]
            M1 = splitline[2]
            R = splitline[3]
            Nsig = splitline[4]
            e_Nsig = splitline[5]
            muillo = splitline[6]
            e_muillo = splitline[7]
            sigma_ = splitline[8]
            e_sigma_ = splitline[9]
            a1 = splitline[10]
            e_a1 = splitline[11]
            n1 = splitline[12]
            e_n1 = splitline[13]
            a2 = splitline[14]
            e_a2 = splitline[15]
            n2 = splitline[16]
            e_n2 = splitline[17]

            # use R=0.1 signals for R<0.1
            if float(R) < 0.1:
                if not modillo+"_M"+str(int(float(M1)))+"_R0p"+str(int(round(float(0.1)+0.00001, 3)*1000)) == signalname:
                    counter += 1
                    continue
                if counter > 0:
                    continue
            else:
                if not modillo+"_M"+str(int(float(M1)))+"_R0p"+str(int(round(float(R)+0.00001, 3)*1000)) == signalname:
                    continue
                if not category == cat:
                    continue

            print(signalname
            print(modillo+"_M"+str(int(float(M1)))+"_R0p"+str(int(round(float(R)+0.00001, 3)*1000))

            signalString = modillo+"_"+"M"+str(int(float(M1)))+"_R0p"+str(int(round(float(R)+0.00001, 3)*1000))+"_"+category

            M2value = float(M1)*float(R)

            ## Signal pdf
            mean = RooRealVar("mean_"+signalString,"mean_"+signalString,float(muillo)) 
            width = RooRealVar("width_"+signalString,"width_"+signalString,float(sigma_), float(muillo)*0.004, float(muillo)*0.4) 
            alpha1 = RooRealVar("alpha1_"+signalString,"alpha1_"+signalString,float(a1))
            n1 = RooRealVar("n1_"+signalString,"n1_"+signalString,float(n1))
            alpha2 = RooRealVar("alpha2_"+signalString,"alpha2_"+signalString,float(a2))
            n2 = RooRealVar("n2_"+signalString,"n2_"+signalString,float(n2))
            nsig = RooRealVar("ParametricSignalPdf_"+signalString+"_norm","ParametricSignalPdf_"+signalString+"_norm",float(Nsig),0,100000000)

            print("ParametricSignalPdf_"+signalString+"_norm", Nsig, nsig.getVal()

            # Evaluate nsig after fit ranges cut
            var_tmp = RooRealVar("var_tmp","var_tmp",1607, 10000)
            var_tmp.setRange("maximum_range", 1607, 10000)  # create range to integrate over
            signalPdf_tmp = RooDoubleCBFast("CB_tmp", "CB_tmp", var_tmp, mean, width, alpha1, n1, alpha2, n2)            
            intrinsicNorm = signalPdf_tmp.createIntegral(ROOT.RooArgSet(var_tmp), RooFit.NormSet(ROOT.RooArgSet(var_tmp)), RooFit.Range("maximum_range")) 
            var_tmp.setRange("fit_range", var_min_set, var_max_set)  # create range to integrate over
            signal_integral = signalPdf_tmp.createIntegral(ROOT.RooArgSet(var_tmp), RooFit.NormSet(ROOT.RooArgSet(var_tmp)), RooFit.Range("fit_range"))
            nsig.setVal( float(Nsig)*signal_integral.getVal()/intrinsicNorm.getVal() )
            
            # print("###########################################################################"
            # print(Nsig, signal_integral.getVal(), intrinsicNorm.getVal(), nsig.getVal()
            # print("###########################################################################"

            icat_i[0] = icat
            M1_i[0]   = float(M1)
            R_i[0]    = float(R)
            nsig_i[0] = nsig.getVal()

            Ntot_i[0] = lumi * 1000
            Seff_i[0] = nsig_i[0]/Ntot_i[0]
            Eeff_i[0] = Seff_i[0]*0.2

            efftree.Fill()

            ## Signal pdf sistematics uncertainties
            #JES_uncertainty = RooRealVar("JES_uncertainty", "JES_uncertainty", 1, 0.916, 1.08)
            #JER_uncertainty = RooRealVar("JER_uncertainty", "JER_uncertainty", 1, 0.5, 1.5)
            meanShape_mu_err_cat = RooRealVar("meanShape_mu_err_"+signalString, "meanShape_mu_err_"+signalString, 1, 0.5, 1.5)
            meanShape_sigma_err_cat = RooRealVar("meanShape_sigma_err_"+signalString, "meanShape_sigma_err_"+signalString, 1, 0.4, 1.6)

            ## Parameters with ALL sistematics uncertainties
            #mean_sist = RooFormulaVar("mean_sist_"+signalString,"mean_sist_"+signalString,"mean_"+signalString+"*JES_uncertainty*meanShape_mu_err_"+signalString, RooArgList(mean,JES_uncertainty,meanShape_mu_err_cat))
            #width_sist = RooFormulaVar("width_sist_"+signalString,"width_sist_"+signalString,"width_"+signalString+"*JER_uncertainty*meanShape_sigma_err_"+signalString, RooArgList(width,JER_uncertainty,meanShape_sigma_err_cat))

            ## no mu syst (independent)
            #mean_sist = RooFormulaVar("mean_sist_"+signalString,"mean_sist_"+signalString,"mean_"+signalString+"*JES_uncertainty", RooArgList(mean,JES_uncertainty))
            #width_sist = RooFormulaVar("width_sist_"+signalString,"width_sist_"+signalString,"width_"+signalString+"*meanShape_sigma_err_"+signalString, RooArgList(width,meanShape_sigma_err_cat))

            ## no sigma syst (independent)
            #mean_sist = RooFormulaVar("mean_sist_"+signalString,"mean_sist_"+signalString,"mean_"+signalString+"*JES_uncertainty*meanShape_mu_err_"+signalString, RooArgList(mean,JES_uncertainty,meanShape_mu_err_cat))
            #width_sist = RooFormulaVar("width_sist_"+signalString,"width_sist_"+signalString,"width_"+signalString+"*JER_uncertainty", RooArgList(width,JER_uncertainty))
            
            ## Parameters with only Jet Corrections uncertainties
            # mean_sist = RooFormulaVar("mean_sist_"+signalString,"mean_sist_"+signalString,"mean_"+signalString+"*JES_uncertainty", RooArgList(mean,JES_uncertainty))
            # width_sist = RooFormulaVar("width_sist_"+signalString,"width_sist_"+signalString,"width_"+signalString+"*JER_uncertainty", RooArgList(width,JER_uncertainty))

            ## Parameters with only independent uncertainties on mu and sigma
            mean_sist  = RooFormulaVar("mean_sist_" +signalString, "mean_sist_" +signalString, "mean_" +signalString+"*meanShape_mu_err_"+signalString   , RooArgList(mean , meanShape_mu_err_cat   ))
            width_sist = RooFormulaVar("width_sist_"+signalString, "width_sist_"+signalString, "width_"+signalString+"*meanShape_sigma_err_"+signalString, RooArgList(width, meanShape_sigma_err_cat))

            signalPdf = RooDoubleCBFast("signalPdf_"+signalString, "signalPdf_"+signalString, var[icat], mean_sist, width_sist, alpha1, n1, alpha2, n2)            
            ParametricSignalPdf_ = RooParametricShapeBinPdf("ParametricSignalPdf_"+signalString, "ParametricSignalPdf_"+signalString, signalPdf, var[icat], RooArgList(mean_sist, width_sist, alpha1, n1, alpha2, n2), th1_rebin[icat])

            
            # signalPdf = RooDoubleCBFast("signalPdf_"+signalString, "signalPdf_"+signalString, var[icat], mean_sist, width_sist, alpha1, n1, alpha2, n2)            
            # ParametricSignalPdf = RooParametricShapeBinPdf("ParametricSignalPdf_"+signalString, "ParametricSignalPdf_"+signalString, signalPdf, var[icat], RooArgList(mean_sist, width_sist, alpha1, n1, alpha2, n2), th1_rebin[icat])

            ## signal pdf without systematic uncertainties
            #signalPdf = RooDoubleCBFast("signalPdf_"+signalString, "signalPdf_"+signalString, var[icat], mean, width, alpha1, n1, alpha2, n2)            
            #ParametricSignalPdf = RooParametricShapeBinPdf("ParametricSignalPdf_"+signalString, "ParametricSignalPdf_"+signalString, signalPdf, var[icat], RooArgList(mean, width, alpha1, n1, alpha2, n2), th1_rebin[icat])

            mean.setConstant(kTRUE)
            width.setConstant(kTRUE)
            alpha1.setConstant(kTRUE)
            n1.setConstant(kTRUE)
            alpha2.setConstant(kTRUE)
            n2.setConstant(kTRUE)
            nsig.setConstant(kTRUE) #signal normalization should be set to constant (required by combine tool)
            
            #signalPdf.Print()
            #nsig.Print()

            ## Fill Workspace for signal
            # Import background pdf and normalization
            getattr(w,'import')(ParametricSignalPdf_)
            #getattr(w,'import')(signalPdf)
            getattr(w,'import')(nsig)
            
            ## Create datacard for current signal
            datacardfilename = sub_outputdirdatacards+"/"+"datacard_"+signalString+".txt"
            outputfile = io.open(datacardfilename,'w')

            outputfile.write( u"## Datacard for "+signalString+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"imax * number of channels"+"\n" )
            outputfile.write( u"jmax * number of backgrounds"+"\n" )
            outputfile.write( u"kmax * number of nuisance parameters (source of systematic uncertainties)"+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"shapes sig "+cat+" "+filenameworkspace+" "+workspaceName+":"+ParametricSignalPdf_.GetName()+"\n" )
            outputfile.write( u"shapes bkg "+cat+" "+filenameworkspace+" "+workspaceName+":"+ParametricBkgPdf[icat].GetName()+"\n" )
            outputfile.write( u"shapes data_obs "+cat+" "+filenameworkspace+" "+workspaceName+":"+rooHist[icat].GetName()+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"bin \t\t"+cat+"\n" )
            outputfile.write( u"observation \t -1"+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"bin \t\t"+cat+"\t"+cat+"\n" )
            outputfile.write( u"process \t"+"sig \t\t"+ "bkg"+"\n" )
            outputfile.write( u"process \t"+"0 \t\t"+ "1"+"\n" )
            outputfile.write( u"rate \t\t"+"1 \t\t"+ "1"+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )

            ## Sistematic uncertainties for mean pulse shape and JEC
            ## Nsig (normalization)
            ## lumi
            outputfile.write( u"lumi lnN 1.018 -\n")

            ## JMS
            JMS_Nsig_down = syst_dict["JMS_down_Nsig"][icat]
            JMS_Nsig_up   = syst_dict["JMS_up_Nsig"  ][icat]
            outputfile.write( u"JMS_uncertainty lnN "+JMS_Nsig_down+"/"+JMS_Nsig_up+" -\n")

            ## JMR
            JMR = syst_dict["JMR_Nsig"][icat]
            outputfile.write( u"JMR_uncertainty lnN 1/"+JMR+" -\n")

            ## Tau21
            T21_shift = syst_dict["T21_shift_Nsig"][icat]
            outputfile.write( u"T21_shift lnN 1/"+T21_shift+" -\n")

            ## Uncertainty on Mean
            # JES + meanShape_mean_err = sqrt( (0.02)^2 + (0.2)^2 ) =~ 0.028
            outputfile.write( u"meanShape_mu_err_"+signalString+" param 1 0.028\n")

            ## Uncertainty on Sigma
            ## JER + meanShape_sigma_err
            #JER_syst_cat = TMath.Sqrt( pow( float(syst_dict["JER_width"][icat])-1 , 2 ) + pow( 0.3, 2 ) ) # JER from file
            JER_syst_cat = TMath.Sqrt( pow( 0.15 , 2 ) + pow( 0.3, 2 ) ) # JER fixed
            outputfile.write( u"meanShape_sigma_err_"+signalString+" param 1 "+str(JER_syst_cat)+"\n")
            
            for param in fitparam:
                outputfile.write( param.GetName()+u" flatParam"+"\n" )

            outputfile.write( nbkg[icat].GetName()+u" flatParam"+"\n" )
            outputfile.close()
            print("Saved datacard "+datacardfilename
        signalInputfile.close()
            
        ## Fill Workspace for data and background
        # Import data
        getattr(w,'import')(rooHist[icat])
        # Import background pdf and normalization
        getattr(w,'import')(ParametricBkgPdf[icat])
        getattr(w,'import')(nbkg[icat])
        
        ############### Plots ##################
        
        ## Create canvas
        canvas[icat] = TCanvas("canvas_"+cat, "canvas_"+cat, 200, 10, 700, 500 )
        fPads1 = TPad("pad1_"+cat, "pad1_"+cat, 0.00, 0.28, 0.99, 0.99)
        fPads2 = TPad("pad2_"+cat, "pad2_"+cat, 0.00, 0.00, 0.99, 0.345)
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
        #th1_rebin[icat].SetMinimum(0.03)
        ymin  = th1_rebin[icat].GetBinContent(th1_rebin[icat].FindBin(var_max_set))*0.1
        if ymin >0:
            th1_rebin[icat].GetYaxis().SetRangeUser( th1_rebin[icat].GetBinContent(th1_rebin[icat].FindBin(var_max_set))*0.1,  th1_rebin[icat].GetMaximum()*5 )
        else:
            th1_rebin[icat].GetYaxis().SetRangeUser( 0.03,  th1_rebin[icat].GetMaximum()*5 )

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
        th1_rebin_bkg[icat].GetYaxis().SetRangeUser( 1000, 100000)
        #th1_rebin_bkg[icat].GetYaxis().SetRangeUser( th1_rebin[icat].GetBinContent(th1_rebin[icat].FindBin(var_max_set)),  th1_rebin[icat].GetMaximum()*5 )
        th1_rebin_bkg[icat].Draw("histsame")

        #draw the lumi text on the canvas
        CMS_lumi.CMS_lumi(fPads1, iPeriod, iPos)    #canvas.Update()
        canvas[icat].Modified()
        canvas[icat].Update()
        
        ## Pad2
        fPads2.cd()
        th1_rebin_pull[icat].GetXaxis().SetTitle(vartitle)
        th1_rebin_pull[icat].GetYaxis().SetTitle("#frac{Bin - Fit}{Uncertainty} ")
        th1_rebin_pull[icat].SetTitle("")
        th1_rebin_pull[icat].SetMinimum(-4)    
        th1_rebin_pull[icat].SetMaximum(4)
        th1_rebin_pull[icat].SetLineColor(2)
        th1_rebin_pull[icat].SetFillColor(2)
        th1_rebin_pull[icat].SetMarkerStyle(20)
        th1_rebin_pull[icat].SetMarkerColor(1)
        th1_rebin_pull[icat].SetStats(0)
        th1_rebin_pull[icat].GetYaxis().SetNdivisions(405, kTRUE)
        th1_rebin_pull[icat].GetXaxis().SetTitleSize(0.16)
        th1_rebin_pull[icat].GetXaxis().SetLabelSize(0.13)
        th1_rebin_pull[icat].GetXaxis().SetTitleOffset(0.83)
        th1_rebin_pull[icat].GetYaxis().SetTitleSize(0.12)
        th1_rebin_pull[icat].GetYaxis().SetLabelSize(0.11)
        th1_rebin_pull[icat].GetYaxis().SetTitleOffset(0.35)
        th1_rebin_pull[icat].Draw("hist")
        
        ## Legend
        fPads1.cd()
        legend = TLegend(0.45, 0.7, 0.7, 0.87)
        legend.SetLineColor(0)
        legend.SetLineWidth(0)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        if(opt.fitData):
            legend.AddEntry(th1_rebin[icat], "Data", "p")
        else:
            legend.AddEntry(th1_rebin[icat], "Toy MC", "p")
        legend.AddEntry(th1_rebin_bkg[icat], "Fit", "l")
        legend.Draw()

        ## Plot fit results
        pt = TPaveText(0.73, 0.58, 0.87, 0.87,"ndc")
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
        t1 = pt.AddText(WWText)
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

        if not pval_ok[0]:
            pt2 = TPaveText(0.15, 0.3, 0.55, 0.15,"ndc")
            pt2.SetFillColor(0)
        
            text = "NOT INCLUDED IN FINAL FIT"
            t1 = pt2.AddText(text)
            t1.SetTextColor(1)
            t1.SetTextSize( 0.06 )
            pt2.Draw("same")
            canvas[icat].Update()
        
        ## Save canvas
        outputfilename_pdf = sub_outputdir+"/"+"c_"+cat+".pdf"
        outputfilename_png = sub_outputdir+"/"+"c_"+cat+".png"
        outputfilename_root = sub_outputdir+"/"+"c_"+cat+".root"
        canvas[icat].SaveAs(outputfilename_pdf)
        canvas[icat].SaveAs(outputfilename_png)
        canvas[icat].SaveAs(outputfilename_root)
        ##os.system("cp "+ outputfilename_pdf + " " + sub_weboutputdir)
        os.system("cp "+ outputfilename_png + " " + sub_weboutputdir)
    
        ##
        rootfile[icat].Close()

        del fitparam[:]

    chi2file[isignal].cd()
    chi2tree[isignal].Write()
    chi2file[isignal].Save()
    chi2file[isignal].Close()

    w.writeToFile(filenameworkspacePath)
    w.Delete()
    print("Workspace at "+filenameworkspacePath

efftree.Write()
efftreefile.Write()
#### Create final datacards (by combinging all categories for a given signal) ####

## Loop over inputdir (find the models)
listOfSignalModels = []
for signalstring in subDirList:
    listOfSignalModels.append(unicode(signalstring, "utf-8"))
#print(listOfSignalModels

## Loop over signals and categories (create combined datacard AND datacard list)
datacardList = io.open(outputdirdatacards+"/datacardList.txt", 'w')

for signal in listOfSignalModels:
    currentPath = os.getcwd() # get current working directory (combineCard only work for absolute path)
    datacardfilenameAll = currentPath+"/"+outputdirdatacards+"/datacard_"+signal+".txt"
    command = ""

    splitline = signal.split("_")
    model = splitline[0]
    M1val = (splitline[1]).strip("M")
    Rval  = "0."+(splitline[2]).strip("R0p")

    # if float(M1val) != 8700:
    #     continue
    # if int(float(M1val) / 100) % 5 != 0:
    #     continue

    ## Loop over event categories    
    categoriesList = next(os.walk(inputdir+"/"+signal))[1]
    ncategories = len(categoriesList)

    for icat, cat in enumerate(categoriesList):
        signalStringCat = signal+"_"+cat
        datacardfilename = currentPath+"/"+outputdirdatacards+"/"+signal+"/"+"datacard_"+signalStringCat+".txt"        

        if icat==0:
            command = "python "+combineCardScriptPath

        if not os.path.isfile(datacardfilename):
            continue        
	
        command += " "+cat+"="+datacardfilename

    command += " > "+datacardfilenameAll

    #print(command
    os.system("cd "+outputdirdatacards+"/ ; "+command+" ; "+ "cd .." )
    print("Created final datacard at "+datacardfilenameAll

    datacardList.write(model+" "+M1val+" "+Rval+" "+datacardfilenameAll+"\n")
    print("Created datacard List at "+outputdirdatacards+"/datacardList.txt"
datacardList.close()


## Make multidimfit with fixed signal strength r=0
print("\n\n#######################################################################################\n"
print("                Starting MultiDimFit in bkg-only hipothesis \n"
print("#######################################################################################\n"

os.system("mkdir -p "+outputdir+"/multidimfit_"+fitFunction_name)
for signal in listOfSignalModels:
    splitline = signal.split("_")
    model = splitline[0]
    M1val = (splitline[1]).strip("M")
    Rval  = "0."+(splitline[2]).strip("R0p")
    # if float(M1val) != 8700:
    #     continue
    # if int(float(M1val) / 100) % 5 != 0:
    #     continue

    datacardPathandName = currentPath+"/"+outputdirdatacards+"/datacard_"+signal+".txt"
    datacardName = "datacard_"+signal
    multidimfit_outputdir = outputdir+"/multidimfit_"+fitFunction_name+"/"+datacardName
    os.system("mkdir -p "+multidimfit_outputdir)
    outputlabel = datacardName+"_toy0_robustHesse0"#_rfixed0"


    rMin = 0
    rMax = 0
    if float(M1val) < 3000:
        rMin = -1
        rMax = 1
    elif float(M1val) < 5000:
        rMin = -0.1
        rMax = 0.1
    elif float(M1val) < 6400:
        rMin = -0.05
        rMax = 0.05
    elif float(M1val) < 7000:
        rMin = 0
        rMax = 0.05
    elif float(M1val) < 9000:
        rMin = 0
        rMax = 0.005

    command = ("combine -M MultiDimFit "+datacardPathandName
               +" --rMin "+str(rMin) 
               +" --rMax "+str(rMax)
               +" -S 1"
               +" --algo singles"
               +" --cl 0.68"
               +" -t 0"
               +" --expectSignal 0"
               +" --name _"+outputlabel
               +" --saveFitResult"
               +" --saveWorkspace"
               +" --trackParameters \"rgx{.*}\""
               +" --robustHesse 0"
               #+" --verbose 9"
               #+" --setParameters r=0"
               #+" --freezeParameters \"rgx{meanShape_mu_err.*}\""
               #+" --freezeParameters r"
               #+" --verbose 1"
           )
    print(command
    os.system(command)
    
    moveFit  = "mv multidimfit_"+outputlabel+".root "+outputdir+"/multidimfit_"+fitFunction_name+"/"+datacardName
    moveTree = "mv higgsCombine_"+outputlabel+".MultiDimFit.mH120.root "+outputdir+"/multidimfit_"+fitFunction_name+"/"+datacardName
    #print(moveFit
    os.system(moveFit)
    #print(moveTree
    os.system(moveTree)

    filenameworkspacePath = str(os.environ['CMSSW_BASE'])+"/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/"+outputdirdatacards+"/"+signal+"/workspace_"+signal+".root"
    multidimfit_plotsdir = outputdir+"/plotSimFit_"+fitFunction_name+"/"+signal
    multidimfit_plotswebdir = opt.weboutputdir+"/plotSimFit_"+fitFunction_name+"/"+signal
    ## Plot multidimfit and compute chi2
    command = ("python plotSimultaneousFit_ParametricShape.py"
               +" -t "+filenameworkspacePath
               +" -f "+multidimfit_outputdir+"/higgsCombine_"+outputlabel+".MultiDimFit.mH120.root"
               +" -n -1"
               +" -c "+outputdirdatacards+"/"+signal
               +" -o "+multidimfit_plotsdir
               +" -w "+multidimfit_plotswebdir
               +" -F "+fitFunction_name
               #+" --draw_limit_exp /afs/cern.ch/work/c/cquarant/Trijet_Analysis/CMSSW_8_1_0_Trijet_test/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_FULLDATA_lowMjjExtended/limits/"
               #+" --draw_limit_obs /afs/cern.ch/work/c/cquarant/Trijet_Analysis/CMSSW_8_1_0_Trijet_test/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data_ARCreview_noBCFilter/limits_std_3par/"
               )
    if opt.fitData:
        command += " --fit_to_data"

    print("--------------------------------------"
    print("      Simultaneous fit plot           "
    print("--------------------------------------"
    print(command
    os.system(command)
    print(""
    print("multidimflit plots at "+multidimfit_plotsdir
    print("multidimflit plots at "+multidimfit_plotswebdir
    print(""


