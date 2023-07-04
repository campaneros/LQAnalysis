import os
import argparse
import itertools
import sys
import optparse
import datetime
import subprocess
import io
from multiprocessing import Pool
import create_workspaces_and_datacards_utils as cwd_utils


from array import array
from glob import glob
from ROOT import *

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi

dict_lambda = {
    "0p1": 0,
    "1p0": 1,
    "1p5": 2,
    "2p0": 3,
}

dict_mass = {
    "700": 0,
    "1000": 1,
    "2000": 2,
    "3000": 3,
    "4000": 4,
    "5000": 5,
}


gROOT.SetBatch(True)
usage = "usage: python bias_study.py"

parser = optparse.OptionParser(usage)

parser.add_option("-o", "--output", dest="outputdir",
                  help="the web output directory. Sub-directories are created automatically.")

parser.add_option("-f", "--file", dest="outputfile",
                  help="name of the output file")

parser.add_option("-t", "--toys", dest="toys",
                  help="toys of the bias study")

parser.add_option("-a", dest="all", default=0,
                  help="all category toghether (1) or not (0)")

(opt, args) = parser.parse_args()


if not opt.outputdir:   
    parser.error('weboutput dir not provided')

if not opt.outputfile:
    parser.error('output file name not provided')







def main():
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    c1 = TCanvas("c1","c1")
    lumi = 138.000
    os.system("mkdir -p "+opt.outputdir)
    CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(lumi)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    iPos = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.15
    iPeriod = 4
    Mass = [700, 1000, 2000, 3000, 4000, 5000]
    L = ['0p1','1p0','1p5','2p0']
    categories = ['category2Muon_BDT_loose_btag','category2Muon_BDT_loose_nobtag','category2Muon_BDT_tight_nobtag','category2Muon_BDT_tight_btag','category1Muon_BDT_loose_btag','category1Muon_BDT_loose_nobtag','category1Muon_BDT_tight_nobtag','category1Muon_BDT_tight_btag']
    #pool = Pool()
    tasks = list(itertools.product(Mass, L))
    #bias = run_bias_study(tasks)
    print(tasks)
    histos = {}
    outputfile = TFile.Open(opt.outputfile+".root", "RECREATE")
    fitting_function = ['STD','UA2']
    for function_gen, function_fit in itertools.product(fitting_function, fitting_function):
        t=opt.toys
        for cat in categories:
            fit_fucntion_name_gen= cwd_utils.nested_dict["function_%s"%function_gen]["%s"%cat]
            fit_fucntion_name_fit= cwd_utils.nested_dict["function_%s"%function_fit]["%s"%cat]
            exp = 0
            histos["%s_%s_%s"%(function_gen, function_fit,cat)] = TH2D("bias_signal_%s_%s_%s"%(function_gen, function_fit,cat), "bias_signal_%s_%s_%s"%(function_gen, function_fit,cat), len(Mass), 0, len(Mass), len(L), 0, len(L))
            histos["hr_%s_%s_%s"%(function_gen, function_fit,cat)] = TH2D("hr_bias_signal_%s_%s_%s"%(function_gen, function_fit,cat), "hr_bias_signal_%s_%s_%s"%(function_gen, function_fit,cat), len(Mass), 0, len(Mass), len(L), 0, len(L))
            for m,l in tasks:
                print(m,l,t,cat)
                if not "%s_%s_%s"%(function_gen, function_fit,cat) in histos:
                    print("cazzo")
                    print(histos["%s_%s_%s"%(function_gen, function_fit,cat)])
                inputfile ="../Fit_Signal/output_MC/umu_gen_"+function_gen+"_fit_"+function_fit+"/datacard_"+fit_fucntion_name_fit+"_umuLQumu_M"+str(m)+"_L"+str(l)+"_"+cat+"_t_"+str(t)+"_syst0_seed123456/higgsCombine_toys"+str(t)+"_expectSignal"+str(exp)+"_std_4par.MultiDimFit.mH120.123456.root"
                cat_name = os.path.dirname(inputfile)
                print(cat_name)
                print(histos["%s_%s_%s"%(function_gen, function_fit,cat)])
                cat_parts = cat_name.split("_")
                if not opt.all:
                    if "2Muon" in cat_name:
                        index = cat_parts.index('category2Muon')
                    else:
                        index = cat_parts.index('category1Muon') 
               #cates = cat_parts[index]+"_"+cat_parts[index+1]+"_"+cat_parts[index+2]+"_"+cat_parts[index+3]
                #print(cat)
                    mmass = float(cat_parts[index-2].replace("M",""))
                else:
                #cates = "all"
                    mmass=   0
                print(inputfile)
                inputfile = TFile.Open(inputfile)
                t_mu=99.
                tree = inputfile.Get("limit")
                if tree:
                    if tree.GetEntries() > 0:
                        h_pull = TH1D("h_pull", "", 32, -4, 4)
                        tree.Draw("(r-"+str(exp)+")/trackedParamErr_r>>h_pull", "trackedParamErr_r>0")
                        if ("2Muon" in cat_name and mmass<2000) or ("1Muon" in cat_name and mmass<3100) or opt.all:
                            h_pull.Fit("gaus")
                            gaus = h_pull.GetListOfFunctions().FindObject("gaus")
                            print(gaus)
                        if ("2Muon" in cat_name and mmass<2000) or ("1Muon" in cat_name and mmass<3100) or opt.all:
                            t_mu     = gaus.GetParameter(1)
                        else:
                            t_mu   = h_pull.GetMean() 
                        hr = TH1D("hr", "", 150000, -10, 10)
                        tree.Draw("r>>hr", "trackedParamErr_r>0")
                        hr.Rebin( int(round( (hr.FindLastBinAbove(0)-hr.FindFirstBinAbove(0))/32 )) )
                        hr.GetXaxis().SetRangeUser( hr.GetMean()-5*hr.GetRMS(), hr.GetMean()+5*hr.GetRMS())
                        appp = hr.GetMean()
                #t_mu     = pt.AddText("#mu    = %s #pm %s" % ( str(round(gaus.GetParameter(1),2)) , str(round(gaus.GetParError(1),2))      )    )
                #t_sigma  = pt.AddText("#sigma = %s #pm %s" % ( str(round(gaus.GetParameter(2),2)) , str(round(gaus.GetParError(2),2))      )    )   
            #print(histos["%s_%s_%s"%(function_gen, function_fit,cat)]) 
            #t_mean   = pt.AddText("mean   = %s " % ( str(round(h_pull.GetMean(), 3))    ) )    
            #t_StdDev = pt.AddText("StdDev = %s " % ( str(round(h_pull.GetStdDev(),2))   ) )    
            #t_r      = pt.AddText("signal inj r = %s"  % ( str( format(float(opt.rExpected), ".2E") ) ) )
                #else:
                #    t_mu = -99
                print(histos["%s_%s_%s"%(function_gen, function_fit, cat)])
                histos["%s_%s_%s"%(function_gen, function_fit, cat)].SetBinContent(dict_mass["%s"%m]+1, dict_lambda["%s"%l]+1, float(t_mu))
                histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].SetBinContent(dict_mass["%s"%m]+1, dict_lambda["%s"%l]+1, float(appp))
                histos["%s_%s_%s"%(function_gen, function_fit, cat)].GetXaxis().SetBinLabel(dict_mass["%s"%m]+1, str(m))
                histos["%s_%s_%s"%(function_gen, function_fit, cat)].GetYaxis().SetBinLabel(dict_lambda["%s"%l]+1, str(l))
                histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].GetXaxis().SetBinLabel(dict_mass["%s"%m]+1, str(m))
                histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].GetYaxis().SetBinLabel(dict_lambda["%s"%l]+1, str(l))
            histos["%s_%s_%s"%(function_gen, function_fit, cat)].GetXaxis().SetTitle("Mass")
            histos["%s_%s_%s"%(function_gen, function_fit, cat)].GetYaxis().SetTitle("Coupling")
            histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].GetXaxis().SetTitle("Mass")
            histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].GetYaxis().SetTitle("Coupling")
            #histos["%s_%s_%s"%(function_gen, function_fit,cat)].SetMaximum(3)
            gStyle.SetPalette(kRainBow) #// set a color palette
            #histos["%s_%s_%s"%(function_gen, function_fit,cat)].SetContour(999) #set a large number of contour levels

            histos["%s_%s_%s"%(function_gen, function_fit,cat)].GetZaxis().SetNdivisions(4)
            histos["hr_%s_%s_%s"%(function_gen, function_fit,cat)].GetZaxis().SetNdivisions(4)
            levels = array('d') #  // define an array for the contour levels
            levels.append(-3) # // first contour level at 0
            #fill the array
            for i in range(0,8):
                levels.append(-1 + (i/8.)*2) # // map [0,1000] to [-5,5]
            levels.append(1) # // last contour level at 1
            levels.append(3) # // last contour level at 10
            histos["%s_%s_%s"%(function_gen, function_fit,cat)].SetContour(len(levels), levels);            
            histos["%s_%s_%s"%(function_gen, function_fit, cat)].GetZaxis().SetRangeUser(-1, 1)

            histos["%s_%s_%s"%(function_gen, function_fit, cat)].SetTitle("Bias %s_%s_%s"%(function_gen, function_fit, cat))
            histos["%s_%s_%s"%(function_gen, function_fit, cat)].Draw("COLZ TEXT")
            for ext in ['png', 'pdf','C']:
                c1.SaveAs("%s/umu_bias_signal_%s_%s_%s.%s"%(opt.outputdir, function_gen, function_fit, cat, ext))
            histos["hr_%s_%s_%s"%(function_gen, function_fit,cat)].SetContour(len(levels), levels);            
            histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].GetZaxis().SetRangeUser(-1, 1)

            histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].SetTitle("hr_Bias %s_%s_%s"%(function_gen, function_fit, cat))
            histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].Draw("COLZ TEXT")
            for ext in ['png', 'pdf','C']:
                c1.SaveAs("%s/hr_umu_bias_signal_%s_%s_%s.%s"%(opt.outputdir, function_gen, function_fit, cat, ext))
            outputfile.cd()
            histos["%s_%s_%s"%(function_gen, function_fit, cat)].Write()
            histos["hr_%s_%s_%s"%(function_gen, function_fit, cat)].Write()

if __name__ == '__main__':
    main()