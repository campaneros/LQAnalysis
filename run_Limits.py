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

gROOT.SetBatch(True)

usage = "usage: python runLimits.py -d /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/datacards/datacardList_Res1ToRes2ToGluGlu.txt -o /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/limit -l Res1ToRes2ToGluGlu"

parser = optparse.OptionParser(usage)

parser.add_option("-d", "--datacard", dest="datacardFileName",
                  help="input list of combine datacards to be processes (ful path is required for each line in the text file)")

parser.add_option("-o", "--output", dest="outputdir",
                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

parser.add_option("-l", "--label", dest="label",
                  help="the name of the output subdirectory")

parser.add_option("-v", "--verbose", dest="verbose", default=1,
                  help="verbosity level")

parser.add_option("--rMin", dest="rMin", default=-0.1,
                  help="rMin value")

parser.add_option("--rMax", dest="rMax", default=0.1,
                  help="rMax value")

parser.add_option("-S", dest="syst", default=0,
                  help="fit systematics")

parser.add_option("-s", dest="single", default=False, action = 'store_true',
                  help="single category limit")


(opt, args) = parser.parse_args()

if not opt.datacardFileName:   
    parser.error('input datacards not provided')

if not opt.outputdir:   
    parser.error('output dir not provided')

if not opt.label:   
    parser.error('label not provided')

##################################################################################################

outputDir = opt.outputdir+"/"+opt.label
#os.system("mkdir -p "+outputDir)

if opt.single:
	lenline = 5
else:
	lenline = 4

datacardFile = io.open(opt.datacardFileName, "r")
for line in datacardFile:

    line = line.rstrip('\n')
    splitline = line.split(" ")
    
    
    if len(splitline) != lenline:
        print("ERROR: wrong file with datacard list")
        sys.exit()
        
    #print splitline
    model = splitline[0]
    mass = splitline[1]
    R = splitline[2].replace(".","p")
    if opt.single:
        category = splitline[3]
        datacard = splitline[4]
    else:
        category = "all"
        datacard = splitline[3]
	
    label = "_"+str(model)+"_"+str(mass)+"_"+str(R)+"_"+category

    ## First run MultiDimFit (to get an hint of the limit)
    #currentScriptDir = os.path.dirname(os.path.abspath(__file__))
    #commandFit = "python "+currentScriptDir+"/"+"runMultiDimFit.py -d "+datacard+" -o "+outputDir+" "+" -t 0 --robustHesse 0 -S 0"
    #print commandFit
    #os.system(commandFit)

    #datacardname = (datacard.split("/")[-1]).split(".")[0]
    #rootFile = outputDir+"/"+datacardname+"_toy0_syst0_robustHesse0"+"/"+"higgsCombine_"+datacardname+"_toy0_syst0_robustHesse0.MultiDimFit.mH120.root"
    ##print rootFile
    #tfileinput = TFile.Open(rootFile)
    #tree = tfileinput.Get("limit")
    #rExpMin = 0
    #rExpMax = 0
    #for event in tree:
    #    r = event.trackedParam_r
    #    er = event.trackedParamErr_r
    #    quantile = event.quantileExpected
    #    if quantile == -1:            
    #        rExpMin = r - 10*er
    #        rExpMax = r + 10*er
    #print rExpMin, rExpMax
    #tfileinput.Close()

    ## Then run AsymptoticLimits    
    command = ("combine -M AsymptoticLimits "
               +datacard
               +" --verbose "+str(opt.verbose)
               +" --name "+label
               +" --rMin "+str(opt.rMin)
               +" --rMax "+str(opt.rMax)               
               +" --rAbsAcc 0.0000001"
               +" -S "+str(opt.syst)
              # +" --expectSignal 1"
               )
    print(command)
    os.system(command)

    os.system("mkdir -p "+outputDir+"/"+R)
    os.system("mkdir -p "+outputDir+"/"+R+"/"+category)

    outputrootfile = "higgsCombine"+label+".AsymptoticLimits.mH120.root"
    commandOutput = "mv "+outputrootfile+" "+outputDir+"/"+R+"/"+category
    os.system(commandOutput) 

print("Output in "+outputDir)

