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
import create_workspaces_and_datacards_utils as cwd_utils

usage ="test"

parser = optparse.OptionParser(usage)


parser.add_option("-t", "--toysfile", dest="toys", default=1,
                  help="number of toys genarted")

(opt, args) = parser.parse_args()


fitFunction_name = "std_4par"
fitparam = []
#gROOT.LoadMacro(os.path.dirname(os.path.abspath(__file__))+"/../../src/libCpp/RooDoubleCBFast.cc+")
gROOT.SetBatch(True)
gErrorIgnoreLevel = kFatal
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

inputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/category"
filenameInput = "test_h1_mmuj_ak4.root"
subDirList = next(os.walk(inputdir))[1]
print(subDirList)

## Output directories  
outputdir = "/data/mcampana/CMS/CMSSW_10_2_13/src/Fit_Signal/output_MC"
outputdirdatacards = "/data/mcampana/CMS/CMSSW_10_2_13/src/Fit_Signal/datacards"
weboutputdir = "/data/mcampana/CMS/CMSSW_10_2_13/src/Fit_Signal/output_plot"
os.system("mkdir -p "+outputdir)
os.system("mkdir -p "+outputdirdatacards)
filenameworkspace = "workspace.root"
filenameworkspacePath = outputdir+"/"+filenameworkspace
workspaceName = "w"


signalInputfilename = "/data/mcampana/CMS/CMSSW_10_2_13/src/Fit_Signal/output/signals.txt"

listOfSignalModels = []

signalInputfile = io.open(signalInputfilename, "r")
for iline, line in enumerate(signalInputfile):
    
    line = line.rstrip('\n')
    splitline = line.split(" ")
    
    modell = splitline[0]
    category = splitline[1]
    M1 = splitline[2]
    L = splitline[3]

    signalStringModel = modell#+"_"+"M"+str(int(float(M1)))+"_R0p"+str(int(float(R)*10))

    if signalStringModel not in listOfSignalModels:
        listOfSignalModels.append(signalStringModel)

signalInputfile.close()
print(listOfSignalModels)
#datacardfilename = outputdirdatacards+"/"+modell+"/categories/"+"datacard_"+signalString+".txt"

print "\n\n#######################################################################################\n"
print "                Starting MultiDimFit \n"
print "#######################################################################################\n"

os.system("mkdir -p "+outputdir+"/multidimfit_"+fitFunction_name)
for signal in listOfSignalModels:
    splitline = signal.split("_")
    model = splitline[0]
    M1val = (splitline[1]).strip("M")
    Lval  = (splitline[2]).strip("L")
    Lval =  Lval.replace("p",".")
    # if float(M1val) != 8700:
    #     continue
    # if int(float(M1val) / 100) % 5 != 0:
    #     continue

    datacardPathandName = outputdirdatacards+"/"+signal+"/datacard_"+signal+".txt"
    datacardName = "datacard_"+signal
    multidimfit_outputdir = outputdir+"/multidimfit_"+fitFunction_name+"/"+datacardName
    os.system("mkdir -p "+multidimfit_outputdir)
    outputlabel = datacardName+"_toy0_robustHesse0"#_rfixed0"
    print(outputlabel)


    rMin = 0
    rMax = 0
    if float(M1val) < 2000:
        rMin = 0
        rMax = 0.01
    elif float(M1val) < 4000:
        rMin = -0.1
        rMax = 0.01
    #elif float(M1val) < 6400:
    #    rMin = -0.05
    #    rMax = 0.05
    #elif float(M1val) < 7000:
    #    rMin = 0
    #    rMax = 0.05
    #elif float(M1val) < 9000:
    #    rMin = 0
    #    rMax = 0.005

    t = int(opt.toys)
    print(t)
    if t != 0:
	fitTreename = "higgsCombine_"+outputlabel+".MultiDimFit.mH120.123456.root"
	expectSignal= 0.0001
	name = "expected_signal"
    	command = ("combine -M MultiDimFit "+datacardPathandName
               +" --rMin "+str(rMin) 
               +" --rMax "+str(rMax)
               +" -S 0"
               #+" --algo singles"
               #+" --cl 0.68"
               +" -t "+str(t)
               +" --expectSignal "+ str(expectSignal)
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
    else:
	fitTreename = "higgsCombine_"+outputlabel+".MultiDimFit.mH120.root"
	expectSignal= 0
	name = "bkg_only"
    	command = ("combine -M MultiDimFit "+datacardPathandName
               +" --rMin "+str(0) 
               +" --rMax "+str(0.01)
               +" -S 0"
               #+" --algo singles"
               #+" --cl 0.68"
               #+" -t "+str(t)
               #+" --expectSignal 0"
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


    print(command)
    os.system(command)
    moveFit  = "mv multidimfit_"+outputlabel+".root "+outputdir+"/multidimfit_"+fitFunction_name+"/"+datacardName
    moveTree = "mv "+fitTreename+" "+outputdir+"/multidimfit_"+fitFunction_name+"/"+datacardName
    #print(moveFit
    os.system(moveFit)
    #print(moveTree
    os.system(moveTree)

    #filenameworkspacePath = outputdir+"/workspace.root"
    multidimfit_plotsdir = outputdir+"/plotSimFit_"+fitFunction_name+"/"+signal
    multidimfit_plotswebdir = weboutputdir+"/plotSimFit_"+fitFunction_name+"/"+signal
    ## Plot multidimfit and compute chi2
    command = ("python plotSimultaneousFit_ParametricShape.py"
               +" -t "+filenameworkspacePath
               +" -f "+multidimfit_outputdir+"/"+fitTreename
               +" -n -1"
               +" -c "+outputdirdatacards+"/"+signal+"/categories"
               +" -o "+multidimfit_plotsdir
               +" -w "+multidimfit_plotswebdir
               +" -F "+fitFunction_name
	       +" -s" +signal
	       +" -b" + name
               #+" --draw_limit_exp /afs/cern.ch/work/c/cquarant/Trijet_Analysis/CMSSW_8_1_0_Trijet_test/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_FULLDATA_lowMjjExtended/limits/"
               #+" --draw_limit_obs /afs/cern.ch/work/c/cquarant/Trijet_Analysis/CMSSW_8_1_0_Trijet_test/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data_ARCreview_noBCFilter/limits_std_3par/"
               )
    #if opt.fitData:
    #    command += " --fit_to_data"

    print("--------------------------------------")
    print("      Simultaneous fit plot           ")
    print("--------------------------------------")
    print(command)
    os.system(command)
    print("")
    print("multidimflit plots at "+multidimfit_plotsdir)
    print("multidimflit plots at "+multidimfit_plotswebdir)
