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

usage = "usage: To be run from trijetana: python bias_study.py "

parser = optparse.OptionParser(usage)

parser.add_option("-d", "--datacard", dest="datacard", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/new_umu_datacards_finalcat_STD_450/umuLQumu_M1000_L1p0/categories/datacard_std_4par_umuLQumu_M1000_L1p0_category1Muon_BDT_tight_btag.txt",
                  help="input combine datacard you want to use for fit")

parser.add_option("-t", dest="toys", default=1,
                  help = "number of toys to run")

parser.add_option("--expectSignal", dest="expectSignal", default=0.001,
                    help="expected signal strenght (r)")

parser.add_option("--rfreeze", dest="freezer", default=0,
                    help="0= r not freeze, !=0 rfreezed to 0")

parser.add_option("-s", "--seed", dest="seed", default=123456,
                    help="seed for toy generation")

parser.add_option("-F", dest="fitFunction", default="std_3par",
                  help="fit function name")

parser.add_option("-o", "--output", dest="outputdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_test/new_umu_datacards_finalcat_STD_450/umuLQumu_M1000_L1p0/categories/",
                  help="the web output directory. Sub-directories are created automatically.")

parser.add_option("-f", "--file", dest="outputfile", default="goodness_of_fit_study.root",
                    help="name of the output file")

parser.add_option("-a", "--algo", dest="algo", default="saturated",
                    help="algo for goodness of fit test")

parser.add_option("--syst", dest="syst", default=0,
                  help="0= no syst, !=0 syst")

(opt, args) = parser.parse_args()


if not opt.datacard:   
    parser.error('datacard not provided')

if not opt.outputdir:   
    parser.error('output dir not provided')

if not opt.toys or int(opt.toys)<=0:
    parser.error('specify a number of toys >0')

#################

datacardname = (opt.datacard.split("/")[-1]).split(".")[0]
print(datacardname)


#mass = float( (datacardname.split("_")[4]).replace("M","") )
genoutputlabel = datacardname+"_genToys_"+str(opt.toys)+"_syst"+str(opt.syst)+"_seed"+str(opt.seed)
outputlabel = datacardname+"_t_"+str(opt.toys)+"_syst"+str(opt.syst)+"_seed"+str(opt.seed)

outputDir = opt.outputdir+"/"+outputlabel
os.system("mkdir -p "+outputDir)


datacard = opt.datacard

pwd = os.environ['PWD']
os.chdir(outputDir)
os.environ['PWD'] = outputDir

if opt.freezer:
    command = ("combine -M GoodnessOfFiti " + str(datacard)
               +" -t "+str(opt.toys)
               +" --algo="+str(opt.algo)
               +" --fixedSignalStrength="+str(opt.expectSignal)
               +" --seed="+str(opt.seed)
               +" -S "+str(opt.syst)
               #+" --toysFreq"
               +" -n _toys%s_expectSignal%s_%s_rfreezed" % ( str(opt.toys), str(opt.expectSignal), str(opt.fitFunction) )
    )
else:
    command = ("combine -M GoodnessOfFit " + str(datacard)
                +" -t "+str(opt.toys)
                +" --algo="+str(opt.algo)
                +" --seed="+str(opt.seed)
                +" -S "+str(opt.syst)
                #+" --toysFreq"
                +" -n _toys%s_expectSignal%s_%s" % ( str(opt.toys), str(opt.expectSignal), str(opt.fitFunction) )
    )
print(command)  
os.system(command)


