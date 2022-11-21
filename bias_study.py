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
import ROOT as ROOT
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi

gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

usage = "usage: To be run from trijetana: python bias_study.py -g /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/datacards_gen/datacard_Res1ToRes2ToGluGlu_M7000_R0p1.txt -d /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/datacards_fit/datacard_Res1ToRes2ToGluGlu_M7000_R0p1.txt -o /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/multiDimFit -t 100 --expectSignal 0 -l 0.001 -S 1"

parser = optparse.OptionParser(usage)

parser.add_option("-g", "--gendatacard", dest="gendatacard",
                  help= "input combine datacard you want to generate from")

parser.add_option("-d", "--datacard", dest="datacard",
                  help="input combine datacard you want to use for fit")

parser.add_option("-o", "--output", dest="outputdir",
                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

parser.add_option("-S", dest="syst", default=0,
                  help="fit systematics")

parser.add_option("-t", dest="toys", default=0, 
                  help="0=data, -1=asimov, >0=run toys")

parser.add_option("-s", "--seed", dest="seed", default=1123456,
                  help="seed for toy generation")

parser.add_option("--expectSignalFile", dest="expectSignal", default=0,
                  help="expected signal strenght (r)")

parser.add_option("-l", dest="limit", default=0.1,
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
#print datacardname

mass = float( (datacardname.split("_")[2]).replace("M","") )

genoutputlabel = datacardname+"_genToys_"+str(opt.toys)+"_syst"+str(opt.syst)+"_seed"+str(opt.seed)
outputlabel = datacardname+"_t_"+str(opt.toys)+"_syst"+str(opt.syst)+"_seed"+str(opt.seed)
outputDir = opt.outputdir+"/"+outputlabel
os.system("mkdir -p "+outputDir)
#os.system("rm -f "+outputDir+"/*")

pwd = os.environ['PWD']
#if opt.gendatacard.startswith("/afs"):
gendatacard = opt.gendatacard
#else:
#    gendatacard = pwd+"/"+opt.gendatacard

#if opt.datacard.startswith("/afs"):
datacard = opt.datacard
#else:
#    datacard = pwd+"/"+opt.datacard

os.chdir(outputDir)
os.environ['PWD'] = outputDir

#############################################
# Generate toys with bkg pdf in gendatacard #
#############################################

command = ("combine -M GenerateOnly -d "+str(gendatacard)
           +" -t "+str(opt.toys)               
           +" -S "+str(opt.syst)
           +" --saveToys"
           +" -n _toys%s_expectSignal%s" % ( str(opt.toys), str(opt.expectSignal) )
           +" -s "+str(opt.seed)
           +" --expectSignal "+str(opt.expectSignal)
           #+" --toysFrequentist "
           #+" --freezeParameters \"rgx{meanShape_sigma_err_*}\",JES_uncertainty"
           #+" --setParameters pdf_index=0"
           #+" --freezeParameters pdf_index"
           )

print command
os.system(command)

limit = float(opt.limit)
if float(opt.expectSignal)==0:
    rMin = -limit*10
    rMax = limit*10

    if mass >= 7000:
        rMin = 0
else:
    rMin = -float(opt.expectSignal)*10
    rMax = float(opt.expectSignal)*10
    #if mass >= 7000:
    #    rMin = 0
#rMin=0 # require rMin>0 ALWAYS, only to rerun failed jobs
#	rMin= -0.01
#	rMax=0.01
#command = ("combine -M FitDiagnostics "+str(datacard)
print "higgsCombine_toys"+str(opt.toys)+"_expectSignal"+str(opt.expectSignal)+".GenerateOnly.mH120."+str(opt.seed)+".root"
command = ("combine -M MultiDimFit "+str(datacard)
           +" -t "+str(opt.toys)
           +" --toysFile higgsCombine_toys"+str(opt.toys)+"_expectSignal"+str(opt.expectSignal)+".GenerateOnly.mH120."+str(opt.seed)+".root"
           +" --rMin "+str(rMin)
           +" --rMax "+str(rMax)
           +" -S "+str(opt.syst)
           +" --trackParameters \"rgx{.*}\""
           +" --robustHesse 1"
           +" --saveFitResult"
           +" --expectSignal "+str(opt.expectSignal)
           +" -n _toys%s_expectSignal%s" % ( str(opt.toys), str(opt.expectSignal) )
	   + " -s "+ str(opt.seed) 
           #+" --cminDefaultMinimizerStrategy=0"
           #+" -v 4"
           #+" --freezeParameters \"rgx{meanShape_sigma_err_*}\""
#           +" --noErrors"
#           +" --setParameters pdf_index=1"
#           +" --freezeParameters r"
          )

print command    
os.system(command)

#moveTree = "mv "+"higgsCombineTest.FitDiagnostics.mH120."+str(opt.seed)+".root"+" "+outputDir
#print moveTree    
#os.system(moveTree)
#moveFitRes = "mv fitDiagnostics.root "+outputDir
#print command
#os.system(moveFitRes)

