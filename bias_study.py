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

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

usage = "usage: To be run from trijetana: python bias_study.py -g /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/datacards_gen/datacard_Res1ToRes2ToGluGlu_M7000_R0p1.txt -d /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/datacards_fit/datacard_Res1ToRes2ToGluGlu_M7000_R0p1.txt -o /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/multiDimFit -t 100 --expectSignal 0 -l 0.001 -S 1"

parser = optparse.OptionParser(usage)

parser.add_option("-g", "--gendatacard", dest="gendatacard", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M1000_L0p1/datacard_gen_LQumu_M1000_L0p1.txt",
                  help= "input combine datacard you want to generate from")

parser.add_option("-d", "--datacard", dest="datacard", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M1000_L0p1/datacard_std_4par_LQumu_M1000_L0p1.txt",
                  help="input combine datacard you want to use for fit")

parser.add_option("-o", "--output", dest="outputdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/", 
                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

parser.add_option("-S", dest="syst", default=0,
                  help="fit systematics")

parser.add_option("-F", dest="fitFunction", default="std_4par",
                  help="fit function name")

parser.add_option("-t", dest="toys", default=1, 
                  help="0=data, -1=asimov, >0=run toys")

parser.add_option("--gen", dest="gentoys", default=0, 
                  help="0= not generate toys, !=0 generate toys")

parser.add_option("--rfreeze", dest="freezer", default=0, 
                  help="0= r not freeze, !=0 rfreezed to 0")


parser.add_option("-s", "--seed", dest="seed", default=123456,
                  help="seed for toy generation")

parser.add_option("--expectSignal", dest="expectSignal", default=0.001,
                  help="expected signal strenght (r)")

parser.add_option("-l", dest="limit", default=0.01,
                  help="expected signal strenght (r)")

(opt, args) = parser.parse_args()

if not opt.datacard:   
    parser.error('datacard not provided')

if not opt.outputdir:   
    parser.error('output dir not provided')

if not opt.toys or int(opt.toys)<=0:
    parser.error('specify a number of toys >0')

##################################################################################################

datacardname = (opt.datacard.split("/")[-1]).split(".")[0]
print(datacardname)

mass = float( (datacardname.split("_")[4]).replace("M","") )

genoutputlabel = datacardname+"_genToys_"+str(opt.toys)+"_syst"+str(opt.syst)+"_seed"+str(opt.seed)
outputlabel = datacardname+"_t_"+str(opt.toys)+"_syst"+str(opt.syst)+"_seed"+str(opt.seed)
outputDir = opt.outputdir+"/"+outputlabel
os.system("mkdir -p "+outputDir)
#os.system("rm -f "+outputDir+"/*")

#if opt.gendatacard.startswith("/afs"):
gendatacard = opt.gendatacard
#else:
#    gendatacard = pwd+"/"+opt.gendatacard

#if opt.datacard.startswith("/afs"):
datacard = opt.datacard
#index = datacard.find("_LQu")
#datacard = datacard[:index]+"_"+opt.fitFunction+datacard[index:]
print(datacard)
#else:
#    datacard = pwd+"/"+opt.datacard

pwd = os.environ['PWD']
os.chdir(outputDir)
os.environ['PWD'] = outputDir

#############################################
# Generate toys with bkg pdf in gendatacard #
#############################################
if opt.gentoys:
    command = ("combine -M GenerateOnly -d "+str(gendatacard)
           +" -t "+str(opt.toys)               
           +" -S "+str(opt.syst)
           +" --saveToys"
           +" -n _toys%s_expectSignal%s" % ( str(opt.toys), str(opt.expectSignal) )
           #+" -s "+str(opt.seed)
           +" --expectSignal "+str(opt.expectSignal)
           #+" --expectSignalMass 0"
	        #+" -v 4"
           #+" --toysFrequentist "
           #+" --freezeParameters \"rgx{meanShape_sigma_err_*}\",JES_uncertainty"
           #+" --setParameters pdf_index=0"
           #+" --freezeParameters pdf_index"
           )

    print(command)
    os.system(command)

limit = float(opt.limit)
if float(opt.expectSignal)==0:
    rMin = -limit*10
    #rMin = 0
    rMax = limit*10

    #if mass >= 2100:
    #    rMin = 0
else:
    rMin = -float(opt.expectSignal)*10
    rMax = float(opt.expectSignal)*10
    #if mass >= 7000:
    #    rMin = 0
#rMin=0 # require rMin>0 ALWAYS, only to rerun failed jobs
#	rMin= -0.01
#	rMax=0.01
print("RMIN ", rMin, "      RMAX",rMax)
#command = ("combine -M FitDiagnostics "+str(datacard)
print("higgsCombine_toys"+str(opt.toys)+"_expectSignal"+str(opt.expectSignal)+".GenerateOnly.mH120."+str(opt.seed)+".root")
if opt.freezer:
    command = ("combine -M MultiDimFit "+str(datacard)
           +" -t "+str(opt.toys)
           +" --toysFile "+outputDir +"/higgsCombine_toys"+str(opt.toys)+"_expectSignal"+str(opt.expectSignal)+".GenerateOnly.mH120."+str(opt.seed)+".root"
           +" --rMin "+str(rMin)
           +" --rMax "+str(rMax)
           +" -S "+str(opt.syst)
           +" --trackParameters \"rgx{.*}\""
           +" --robustHesse 1"
           +" --saveFitResult"
          # +" --expectSignal "+str(opt.expectSignal)
           +" -n _toys%s_expectSignal%s_%s_rfreezed" % ( str(opt.toys), str(opt.expectSignal), str(opt.fitFunction) )
            #+ " -s "+ str(opt.seed) 
           #+" --cminDefaultMinimizerStrategy=0"
           #+" -v 4"
           #+" --freezeParameters \"rgx{meanShape_sigma_err_*}\"""
#           +" --noErrors"
           +" --setParameters r=0"
           +" --freezeParameters r"
          )
else:
    command = ("combine -M MultiDimFit "+str(datacard)
           +" -t "+str(opt.toys)
           +" --toysFile "+outputDir +"/higgsCombine_toys"+str(opt.toys)+"_expectSignal"+str(opt.expectSignal)+".GenerateOnly.mH120."+str(opt.seed)+".root"
           +" --rMin "+str(rMin)
           +" --rMax "+str(rMax)
           +" -S "+str(opt.syst)
           +" --trackParameters \"rgx{.*}\""
           +" --robustHesse 1"
           +" --saveFitResult"
          # +" --expectSignal "+str(opt.expectSignal)
           +" -n _toys%s_expectSignal%s_%s" % ( str(opt.toys), str(opt.expectSignal), str(opt.fitFunction) )
            #+ " -s "+ str(opt.seed) 
           #+" --cminDefaultMinimizerStrategy=0"
           #+" -v 4"
           #+" --freezeParameters \"rgx{meanShape_sigma_err_*}\"""
#           +" --noErrors"
           #+" --setParameters r=0"
           #+" --freezeParameters r"
          )
print(command)  
os.system(command)

#moveTree = "mv "+"higgsCombineTest.FitDiagnostics.mH120."+str(opt.seed)+".root"+" "+outputDir
#print moveTree    
#os.system(moveTree)
#moveFitRes = "mv fitDiagnostics.root "+outputDir
#print command
#os.system(moveFitRes)

