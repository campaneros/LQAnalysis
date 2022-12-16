from os import *
import os
import sys
import optparse
import datetime
import subprocess
import io
import ROOT as ROOT
import re

from array import array
from glob import glob
#from ROOT import *
#from ROOT import AddressOf
#from ROOT import ROOT.RooFit

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

usage = "usage: python plotSimFit.py -t toysfile/workspacefile -f fitFile -n 1 -c catdir -o outputdir -w weboutputdir -F fitFunction --draw_limit_(exp/obs) limit_dir"

parser = optparse.OptionParser(usage)




parser.add_option("-t", "--toysfile", dest="toysfile", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_LQumu_M3000_L1p0_t_1000_syst0_seed123456/higgsCombine_toys1000_expectSignal0.GenerateOnly.mH120.123456.root",
                  help="input file with fitted toys")

parser.add_option("-f", "--fitFile", dest="fitFile", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_LQumu_M3000_L1p0_t_1000_syst0_seed123456/higgsCombine_toys1000_expectSignal0.0_gen.MultiDimFit.mH120.123456.root",
                  help="input file with tree of post-fit parameters.")

parser.add_option("-c", "--catdir", dest="catdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M3000_L1p0/categories/",
                  help="name of directory containing categories dirs")
                  
parser.add_option("-o", "--outputdir", dest="outputdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/plotSimFit_std_4par/LQumu_M1000_L0p1_bis_test",
                  help="name of the output directory")

parser.add_option("-b", "--outputFilename", dest="outputFile", default="expect_signal",
                  help="name of the output file")

parser.add_option("-w", "--weboutputdir", dest="weboutputdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_plot/plotSimFit_std_4par/LQumu_M1000_L0p1_bis",
                  help="name of the web output directory")

parser.add_option("-s", "--signaname", dest="signaname", default="LQumu_M3000_L1p0",
                  help="signale name model")

parser.add_option("-n", dest="nToy", default=1,
                  help="number of the toy to plot")

parser.add_option("-F", dest="fitFunction", default="std_3par",
                  help="fit function name")

parser.add_option("--fit_to_data", action="store_true", dest="fitData",
                  help="Fit data. If not specified, fit will be performed on toy thrown from MC")

parser.add_option("--draw_limit_exp", dest="draw_limit_exp", default=False,
                  help="Draw signal with expected limit on cross section")

parser.add_option("--draw_limit_obs", dest="draw_limit_obs", default=False,
                  help="Draw signal with observed limit on cross section")
parser.add_option("--save", dest="save", default=False,
                  help="Draw signal with observed limit on cross section")
parser.add_option("--nc", dest="ncat", default="all",
                  help="fcategories to be run")
parser.add_option("--debug", action="store_true", dest="debug", default=False,
                  help="fcategories to be run")


(opt, args) = parser.parse_args()

if not opt.toysfile:   
    parser.error('input toy file not provided')

if not opt.fitFile:   
    parser.error('input fit file not provided')

if not opt.outputdir:
    parser.error('output dir name not provided')

if not opt.weboutputdir:
    parser.error('web output dir name not provided')

######################################################################

## CMS_lumi variables (see CMS_lumi.py)
#lumi = 4.0
lumi = 140.000
CMS_lumi.lumi_13TeV = "%.0f fb^{-1}"%(lumi)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4


os.system("mkdir -p "+opt.outputdir)

#Define binning
#varBins_all = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058,1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869,5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7150, 7500, 7850, 8250, 8650, 8999, 9500, 9999]
varBins_all = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058,1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869,5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250, 7500, 7750, 8000]

print(len(varBins_all))

fitFunction = opt.fitFunction

#Load toy data histogram
nToy = int(opt.nToy)
toysfilenamePath = opt.toysfile
toysfilename     = toysfilenamePath.split("/")[-1]
toysfile = ROOT.TFile.Open(toysfilenamePath)
#print(toysfile)



#signalname  = (toysfilename.replace(".root", "")).replace("workspace_", "")
signalname = opt.signaname
signalname  = re.sub( "_toy[0-9999999]", "", signalname )
splitsignal = signalname.split("_")
model       = splitsignal[0]
M1          = float(splitsignal[1].replace("M",""))
L           = splitsignal[2].strip("L")
L           = float(L.replace("p","."))

#workspacename = toysfilename.replace("workspace_","w_").replace(".root","")

workspacePath = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/workspace.root"
#workspacePath = workspacePath.split("/")[-1]
workspacefile = ROOT.TFile.Open(workspacePath)

workspacename = "w"
workspace = ROOT.RooWorkspace()

workspace=workspacefile.Get(workspacename)

#toy = toysfile.Get("toys/toy_1")
#print(toy)
#else:
#    toy = toysfile.Get("toys/toy_"+str(nToy) )

#Create categories list
categoriesList = next(os.walk(opt.catdir))[2]
for icat, cat in enumerate(categoriesList):
    categoriesList[icat] = cat.replace("datacard_Res1ToRes2ToGluGlu_"+splitsignal[1]+"_"+splitsignal[2]+"_", "").replace(".txt", "")
ncategories = len(categoriesList)

#catType = 0

#if M2val <= 600:
#    catType = 1
#elif M2val <= 1200:
#    catType = 2
#else:
#    catType = 4

#Create output rooROOT.TFile with chi2
indexcat = array('i', [0 ])
Chi2     = array('d', [0.])
Ndof     = array('d', [0.])
Mass  = array('d', [0.])
L    = array('d', [0.])

GlobalChi2 = array('d', [0. ])
GlobalNdof = array('d', [-1.]) #subtract r parameter
CombineGoF = array('d', [0. ])
CombineDof = array('d', [-1.]) #subtract r parameter

chi2file = ROOT.TFile.Open(opt.outputdir+"/test_statistics.root", "RECREATE")

chi2tree = ROOT.TTree("chi2tree","chi2tree")
icatbranch = chi2tree.Branch("icat", indexcat, "icat/I")
chi2branch = chi2tree.Branch("chi2", Chi2, "chi2/D")
ndofbranch = chi2tree.Branch("ndof", Ndof, "ndof/D")
chi2branch = chi2tree.Branch("Mass", Mass, "Mass/D")
ndofbranch = chi2tree.Branch("L", L, "L/D")

globchi2tree = ROOT.TTree("globchi2tree","globchi2tree")
globchi2branch = globchi2tree.Branch("globchi2"  , GlobalChi2, "globchi2/D"  )
globndofbranch = globchi2tree.Branch("globndof"  , GlobalNdof, "globndof/D"  )
combiGoFbranch = globchi2tree.Branch("CombineGoF", CombineGoF, "CombineGoF/D")
combiDoFbranch = globchi2tree.Branch("CombineDof", CombineDof, "CombineDof/D")

outputrooTFile = [None] * len(categoriesList)

if opt.debug:
    ROOT.gROOT.LoadMacro("./Loop_toy.C++")
else:
    ROOT.gROOT.LoadMacro("./Loop_toy.C")


for icat,cat in enumerate(categoriesList):
        if "workspace_"  in cat:
            continue
        elif "gen_" in cat:
                continue
        outputrooTFile[icat] = ROOT.TFile.Open(opt.outputdir+"/"+cat+".root","RECREATE")
        outputrooTFile[icat].cd()
        indexcat[0]+=1

        #Define mjj category range
        splitname = cat.split("_")
        Mass[0] = float(splitname[2].strip("M"))
        L[0] = float(splitname[3].strip("L").replace("p","."))

        app=(cat.split("_"))[4]
        print(app)
        sign_=(cat.strip("datacard_"))

        varname = "m_muj_ak4_"+app
        print(app, sign_)


        ROOT.Make_x2(toysfilenamePath, workspacePath, opt.fitFile, nToy, varname, app, sign_)