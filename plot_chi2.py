
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

usage = "usage: python plotSimFit.py -t toysfile/workspacefile -f fitFile -n 1 -c catdir -o outputdir -w weboutputdir -F fitFunction --draw_limit_(exp/obs) limit_dir"

parser = optparse.OptionParser(usage)

parser.add_option("-F", dest="fitFunction", default="std_4par",
                  help="fit function name")

parser.add_option("-s", dest="start_range", default="450",
                    help="start range for graph")

(opt, args) = parser.parse_args()


#fitFunction_name = "std_4par"
fitparam = []
#gROOT.LoadMacro(os.path.dirname(os.path.abspath(__file__))+"/../../src/libCpp/RooDoubleCBFast.cc+")
ROOT.gROOT.SetBatch(True)
gErrorIgnoreLevel = ROOT.kFatal
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

inputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna_new/src/RootTreeAnalyzer/all_years/category_BDT_data_all/"
filenameInput = "test_h1_mmuj_ak4.root"
subDirList = next(os.walk(inputdir))[1]
print(subDirList)

pval_dict= {}
pvalall_dict= {}
chi2_dict= {}
chi2all_dict= {}
chi2_reduced_dict = {}
chi2all_reduced_dict = {}
ndof_dict = {}
ndofall_dict = {}

canvas_dict = {}


## Output directories
range_list        = [ 300, 320, 350, 385, 400, 450, 480, 525,]
fit_function_list = ["std_2par","std_3par","std_4par","std_5par"]
histos={}
for ilists,lists in enumerate(fit_function_list):
    canvas_dict["%s"%lists] = ROOT.TCanvas("canvas_%s"%lists,"canvas_%s"%lists, 200,10,900,500)
    fitFunction_name = lists
    outputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/output_MC_"+str(lists)
    outputdirdatacards = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/datacards_"+str(lists)
    weboutputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/output_plot_"+str(lists)
    os.system("mkdir -p "+outputdir)
    os.system("mkdir -p "+outputdirdatacards)
    #os.system(u"rm -f "+outputdirdatacards+"/*")
    os.system("mkdir -p "+weboutputdir)
    scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
    
    
    
    
    
    
    
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
    
    
    # Definition of arrays
    

    chi2_dict ["%s"%lists] = ROOT.TH2D("%s"%lists,"%s"%lists, len(subDirList), 0, len(subDirList), len(range_list), 0,len(range_list))
    chi2_reduced_dict ["%s"%lists] = ROOT.TH2D("%s_reduced"%lists,"%s_reduced"%lists, len(subDirList), 0, len(subDirList), len(range_list), 0,len(range_list))
    for index,rang in enumerate(range_list):
        outroot=ROOT.TFile.Open("/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/output_MC_"+str(rang)+"std/outFile_fisher.root")
        histos["%s"%lists]=outroot.Get("chi2_%s"%lists) 
        histos["%s_reduced"%lists]=outroot.Get("chi2red_%s"%lists) 
        for icat, cat in enumerate(subDirList):
            chi2_dict ["%s"%lists].SetBinContent(icat+1,index+1, histos["%s"%lists].GetBinContent(icat+1))
            chi2_dict ["%s"%lists].GetXaxis().SetBinLabel(icat+1,(cat.strip("category")).replace("_BDT_","_"))
            chi2_dict ["%s"%lists].GetYaxis().SetBinLabel(index+1,str(rang))
            chi2_reduced_dict ["%s"%lists].SetBinContent(icat+1,index+1, histos["%s_reduced"%lists].GetBinContent(icat+1))
            chi2_reduced_dict ["%s"%lists].GetXaxis().SetBinLabel(icat+1,(cat.strip("category")).replace("_BDT_","_"))
            chi2_reduced_dict ["%s"%lists].GetYaxis().SetBinLabel(index+1,str(rang))
        outroot.Close()




for ilists,lists in enumerate(fit_function_list):
    outputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/output_MC_"+str(lists)
    outputdirdatacards = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/datacards_"+str(lists)
    weboutputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/output_plot_"+str(lists)
    canvas_dict["%s"%lists] = ROOT.TCanvas("canvas_%s"%lists,"canvas_%s"%lists, 200,10,900,500)
    chi2_dict ["%s"%lists].Draw("COLZ TEXT")
    for ext in [".png",".pdf"]:
        canvas_dict["%s"%lists].SaveAs(outputdir+"/"+str(lists)+ext)
    chi2_reduced_dict ["%s"%lists].Draw("COLZ TEXT")
    for ext in [".png",".pdf"]:
        canvas_dict["%s"%lists].SaveAs(outputdir+"/chi2red"+str(lists)+ext)
