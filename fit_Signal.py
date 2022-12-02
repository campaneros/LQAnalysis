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
from ROOT import *

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi

#gROOT.LoadMacro(os.path.dirname(os.path.abspath(__file__))+"/../../src/libCpp/RooDoubleCBFast.cc+")
gROOT.SetBatch(True)

## Input directories (each folder contains a root file with histograms)
inputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/category"
subDirList = next(os.walk(inputdir))[1]
print(subDirList)

## Output directories  
outputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output"
weboutputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output"
os.system("mkdir -p "+outputdir)
os.system("mkdir -p "+weboutputdir)
scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
os.system("cp "+scriptsPath+"/index.php "+weboutputdir)

## Output txt
signals_filename = "signals.txt"

## Histogram
varname = "m_muj_ak"
vartitle = "m_muj_ak [GeV]"

## Signal Histogram
histoname_signal = ["h1_mmuj_ak4__LQumu_M3000_L1p0", "h1_mmuj_ak4__LQumu_M3000_L0p1", "h1_mmuj_ak4__LQumu_M2000_L1p0", "h1_mmuj_ak4__LQumu_M2000_L0p1", "h1_mmuj_ak4__LQumu_M1000_L1p0", "h1_mmuj_ak4__LQumu_M1000_L0p1"]

## Fit ranges (!!should match the categories!!)
ncategories = len(subDirList)
var_min = []
var_max = [] 

for cat in subDirList:
    var_min.append(453)
    var_max.append(4000)

## Binning
varBins_all = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250, 7500, 7750, 8000]
NvarBins_all = len(varBins_all)-1

## CMS_lumi variables (see CMS_lumi.py)
#lumi = 40.926
lumi = 140
CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(lumi)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

##############################################################################################

# Definition of arrays
rootfile = [None] * ncategories
var = [None] * ncategories

## Output txt file with signals

outputfilename = outputdir+"/"+signals_filename
outputfile = io.open(outputfilename,'w')

## Loop over event categories  
for icat, cat in enumerate(subDirList):

    #if icat > 1:
    #    continue

    print("\n")
    print("######################################")
    print(icat)
    print("######################################")

    ## Input Rootfile
    rootfilename = inputdir+"/"+cat+"/"+"test_h1_mmuj_ak4.root"
    print(rootfilename)
    #print rootfilename    
    rootfile[icat] = TFile.Open(rootfilename)

    ## Main physics observable defined in fit range
    var[icat] = RooRealVar(varname+"_"+cat,vartitle,var_min[icat],var_max[icat])
    var[icat].Print()

    ## Modify variable binning
    varBins = []
    for border in varBins_all:
        #print border
        if border<var_min[icat]:
            continue
        else:
            varBins.append(float(border))
    NvarBins = len(varBins)-1            
    #print NvarBins
    #print varBins

    ## Loop over signals
    for isignal, sighistname in enumerate(histoname_signal):
        #os.system("mkdir -p "+weboutputdir+"/"+sighistname+"/")
        #print isignal, sighistname
        sighistname_split = sighistname.split("__")
        model = sighistname_split[1]
        M1value = (sighistname_split[1].split("_")[1]).split("M")[1]
        Lvaluep = (sighistname_split[1].split("_")[2]).split("L")[1]
        Lvalue = Lvaluep.replace("p",".")
        print("Signal model "+model+" with MRes1="+str(M1value)+" and L="+str(Lvalue))

        ## Get original TH1 histogram from root file
        print("Get "+sighistname+" from file "+rootfilename)
        th1_fromFile_signal = rootfile[icat].Get(sighistname) # 1 GeV bin histogram
        outRoot=TFile(outputdir+"/"+sighistname+"_"+cat+".root","recreate")
        th1_fromFile_signal.Write()	
        outRoot.Close()

        ## Create RooDataHist for signal
        signalString = model+"_M"+str(M1value)+"_"+str(Lvaluep)+"_"+cat
        rooHist_signal = RooDataHist("RooDataHist"+"_"+signalString,
                                     "RooDataHist"+"_"+signalString,
                                     RooArgList(var[icat]),
                                     RooFit.Import(th1_fromFile_signal)
                                     )
        rooHist_signal.Print()

        ## Signal pdf
        mean = RooRealVar("mean_"+signalString,"mean_"+signalString,float(M1value),float(M1value)-300,float(M1value)+300) 
        width = RooRealVar("width_"+signalString,"width_"+signalString,float(M1value)*0.1,float(M1value)*0.1*0.1,float(M1value)*0.1*3) 
        alpha1 = RooRealVar("alpha1_"+signalString,"alpha1_"+signalString,1,0,5)
        n1 = RooRealVar("n1_"+signalString,"n1_"+signalString,1,0,50)
        alpha2 = RooRealVar("alpha2_"+signalString,"alpha2_"+signalString,1,0,5)
        n2 = RooRealVar("n2_"+signalString,"n2_"+signalString,5,0,15)
        nsig = RooRealVar("signalPdf_"+signalString+"_norm","signalPdf_"+signalString+"_norm",1000,0,1000000)
        signalPdf = RooDoubleCBFast("signalPdf_"+signalString, "signalPdf_"+signalString, var[icat], mean, width, alpha1, n1, alpha2, n2)
        
        # Signal extended pdf
        signalExtPdf = RooExtendPdf("signalExtPdf_"+signalString,"signalExtPdf_"+signalString,signalPdf,nsig) 
        
        ## Fit to signal
        if( float(M1value) > var_min[icat] and float(M1value) < var_max[icat]):
            fitResult_signal = signalExtPdf.fitTo(rooHist_signal, 
                                                  RooFit.Strategy(1),
                                                  RooFit.Range(var_min[icat],var_max[icat]),
                                                  RooFit.SumW2Error(kTRUE),
                                                  RooFit.Save(kTRUE),
                                                  RooFit.Verbose(kFALSE),
                                                  RooFit.Extended(kTRUE)
                                                  )
            
            fitResult_signal.Print()

        ## Draw plot
        canvas = TCanvas("canvasSig_"+signalString, "canvasSig_"+signalString, 200, 10, 700, 500 )
        frame = var[icat].frame()
        #rooHist_signal.plotOn(frame,RooFit.Binning(100))
        #if( float(M1value) > var_min[icat] and float(M1value) < var_max[icat]):
        #    signalPdf.plotOn(frame,RooFit.Binning(100))
        #frame.Draw()
        frame = var[icat].frame()
        var[icat].setRange("signal",  max(float(M1value)*0.65, var_min[icat]), min(float(M1value)*1.45, 4000) )
        var[icat].setRange("full", var_min[icat], var_max[icat] )
        rooHist_signal.plotOn(frame,RooFit.Binning(100), RooFit.Range("full"))
        if( float(M1value) > var_min[icat] and float(M1value) < var_max[icat]):
            signalExtPdf.plotOn(frame,RooFit.Binning(100), RooFit.Range("signal"))#,RooFit.VisualizeError(fitResult_signal))
        Chi2 = frame.chiSquare()
        frame.Draw()
        frame.SetAxisRange(0.001,20000,"Y")
        canvas.SetLogy(1)

        pt = TPaveText(0.4, 0.1, 0.6, 0.3,"ndc")
        pt.SetFillColor(0)
        
        Chi2Text = "#chi^{2}/ndof="+str(round(Chi2,2))
        t1 = pt.AddText("mean "+str(round(mean.getValV(),2))+"#pm"+str(round(mean.getError(),2)))
        t1.SetTextColor(1)
        t1.SetTextSize( 0.02 )
        t2 = pt.AddText("width "+str(round(width.getValV(),2))+"#pm"+str(round(width.getError(),2)))
        t2.SetTextColor(1)
        t2.SetTextSize( 0.02 )
        t3 = pt.AddText("alpha1 "+str(round(alpha1.getValV(),2))+"#pm"+str(round(alpha1.getError(),2)))
        t3.SetTextColor(1)
        t3.SetTextSize( 0.02 )
        t4 = pt.AddText("n1 "+str(round(n1.getValV(),2))+"#pm"+str(round(n1.getError(),2)))
        t4.SetTextColor(1)
        t4.SetTextSize( 0.02 )
        t5 = pt.AddText("alpha2 "+str(round(alpha2.getValV(),2))+"#pm"+str(round(alpha2.getError(),2)))
        t5.SetTextColor(1)
        t5.SetTextSize( 0.02 )
        t6 = pt.AddText("n2 "+str(round(n2.getValV(),2))+"#pm"+str(round(n2.getError(),2)))
        t6.SetTextColor(1)
        t6.SetTextSize( 0.02 )

        pt.Draw("SAME")
        
        # integral inside fit range
        var[icat].setRange("my_range", var_min[icat], var_max[icat])  # create range to integrate over
        intrinsicNorm = signalPdf.createIntegral(RooArgSet(var[icat]), RooArgSet(var[icat]), "my_range")
        #print(var_min[icat], var_min_eff[icat])
        var[icat].setRange("CBfit_range", max(float(M1value)*0.5, var_min[icat]), min(float(M1value)*1.4, 4000) )  # create range to integrate over
        CBfitNorm = signalPdf.createIntegral(RooArgSet(var[icat]), RooArgSet(var[icat]), "CBfit_range")
        #var[icat].setRange("my_range_2", var_min_eff[icat], var_max_eff[icat])  # create range to integrate over
        #integral = signalPdf.createIntegral(RooArgSet(var[icat]), RooArgSet(var[icat]),"my_range_2")

        print(nsig.getVal())
        print(rooHist_signal.sumEntries(varname+"_"+cat+" > 0", "my_range"))
        print(rooHist_signal.sumEntries(varname+"_"+cat+" > 0", "CBfit_range"))
        print(rooHist_signal.sumEntries(varname+"_"+cat+" > 0", "my_range_2"))
        NSIG = rooHist_signal.sumEntries(varname+"_"+cat+" > 0", "my_range_2")




        ## Save canvas
        outputfilename_pdf = outputdir+"/"+"cSig_"+signalString+".pdf"
        outputfilename_png = outputdir+"/"+"cSig_"+signalString+".png"
        outputfilename_root = outputdir+"/"+"cSig_"+signalString+".root"
        canvas.SaveAs(outputfilename_pdf)
        canvas.SaveAs(outputfilename_png)
        canvas.SaveAs(outputfilename_root)
        #os.system("cp "+ outputfilename_pdf + " " + weboutputdir+"/"+sighistname)
        #os.system("cp "+ outputfilename_png + " " + weboutputdir+"/"+sighistname)
        #os.system("cp "+scriptsPath+"/index.php "+weboutputdir+"/"+sighistname)

        print(model)
        ## Write signal parameters on txt file
        if( float(M1value) > var_min[icat] and float(M1value) < var_max[icat]):
            outputfile.write(model+" "
                             + cat+" "
                             + str(float(M1value))+" "
                             + str(float(Lvalue))+" "
                             + str(nsig.getValV())+" "+ str(nsig.getError())+" "
                             + str(mean.getValV())+" "+ str(mean.getError())+" "
                             + str(width.getValV())+" "+ str(width.getError())+" "
                             + str(alpha1.getValV())+" "+ str(alpha1.getError())+" "
                             + str(n1.getValV())+" "+ str(n1.getError())+" "
                             + str(alpha2.getValV())+" "+ str(alpha2.getError())+" "
                             + str(n2.getValV())+" "+ str(n2.getError())
                             +u"\n")
    rootfile[icat].Close()


outputfile.close()
print("All signal parameters are in: "+outputfilename)



