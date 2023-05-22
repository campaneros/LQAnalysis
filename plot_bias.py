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

usage = "usage: python plot_bias.py -i /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/bias_f1_vs_f2/fitDiagnostics.root -o /afs/cern.ch/user/c/cquarant/trijet_boosted/Res1ToRes2ToGluGlu/bias_study/pull -f gen_f1_fit_f2 -r 0"

parser = optparse.OptionParser(usage)

parser.add_option("-i", "--inputfile", dest="inputfile",
                  help="input file with fitDisgnostics results")

parser.add_option("-o", "--weboutput", dest="weboutputdir",
                  help="the web output directory. Sub-directories are created automatically.")

parser.add_option("-f", "--file", dest="outputfile",
                  help="name of the output file")
                  
parser.add_option("-r", dest="rExpected", default=0,
                  help="expected signal strenght (r)")

parser.add_option("-a", dest="all", default=0,
                  help="all category toghether (1) or not (0)")

(opt, args) = parser.parse_args()

if not opt.inputfile:   
    parser.error('input file not provided')

if not opt.weboutputdir:   
    parser.error('weboutput dir not provided')

if not opt.outputfile:
    parser.error('output file name not provided')

if not opt.rExpected:   
    parser.error('expected value of r (toys), not provided')

## CMS_lumi variables (see CMS_lumi.py)
lumi = 138.000
CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(lumi)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.15
iPeriod = 4

######################################################################
cat_name = os.path.dirname(opt.inputfile)
print(cat_name)
cat_parts = cat_name.split("_")
if not opt.all:
    if "2Muon" in cat_name:
        index = cat_parts.index('category2Muon')
    else:
        index = cat_parts.index('category1Muon') 

    cat = cat_parts[index]+"_"+cat_parts[index+1]+"_"+cat_parts[index+2]+"_"+cat_parts[index+3]
    print(cat)
    Mass = float(cat_parts[index-2].replace("M",""))
else:
    cat = "all"
    Mass = 0

index = cat_name.index('umuLQumu')


inputfile = TFile.Open(opt.inputfile)
tree = inputfile.Get("limit")

h_pull = TH1D("h_pull", "", 32, -4, 4)
tree.Draw("(r-"+str(opt.rExpected)+")/trackedParamErr_r>>h_pull", "trackedParamErr_r>0")
if ("2Muon" in cat_name and Mass<2000) or ("1Muon" in cat_name and Mass<3100) or opt.all:
        h_pull.Fit("gaus")
        gaus = h_pull.GetListOfFunctions().FindObject("gaus")



gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
c1 = TCanvas("c1","c1")
c1.SetBottomMargin(0.12)
c1.SetLeftMargin(0.12)

h_pull.GetYaxis().SetTitle("toys")
h_pull.GetYaxis().SetTitleSize(0.06)
h_pull.GetYaxis().SetTitleOffset(0.7)
h_pull.GetXaxis().SetRangeUser( h_pull.GetMean()-5*h_pull.GetRMS(), h_pull.GetMean()+5*h_pull.GetRMS())

h_pull.GetXaxis().SetTitle("(#sigma-#sigma_{inj})/#sigma_{Err}")
h_pull.GetXaxis().SetTitleSize(0.06)
h_pull.GetXaxis().SetTitleOffset(0.7)
h_pull.Draw()

pt = TPaveText(0.65, 0.55, 0.89, 0.87,"ndc")
pt.SetFillColor(0)
pt.SetLineColor(1)
pt.SetLineStyle(1)
        
t1 = pt.AddText("toys number: "+str(h_pull.GetEntries()))
t1.SetTextColor(1)
t1.SetTextSize( 0.04 )

if ("2Muon" in cat_name and Mass<2000) or ("1Muon" in cat_name and Mass<3100) or opt.all:
    t_mu     = pt.AddText("#mu    = %s #pm %s" % ( str(round(gaus.GetParameter(1),2)) , str(round(gaus.GetParError(1),2))      )    )
    t_sigma  = pt.AddText("#sigma = %s #pm %s" % ( str(round(gaus.GetParameter(2),2)) , str(round(gaus.GetParError(2),2))      )    )
    t_mu.SetTextColor(1)
    t_mu.SetTextSize( 0.04 )
    t_sigma.SetTextColor(1)
    t_sigma.SetTextSize( 0.04 )
t_mean   = pt.AddText("mean   = %s " % ( str(round(h_pull.GetMean(), 3))    ) )    
t_StdDev = pt.AddText("StdDev = %s " % ( str(round(h_pull.GetStdDev(),2))   ) )    
t_r      = pt.AddText("signal inj r = %s"  % ( str( format(float(opt.rExpected), ".2E") ) ) )



t_mean.SetTextColor(1)
t_mean.SetTextSize( 0.04 )
t_StdDev.SetTextColor(1)
t_StdDev.SetTextSize( 0.04 )
t_r.SetTextColor(1)
t_r.SetTextSize( 0.04 )

pt.Draw("SAME")

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c1, iPeriod, iPos)    #canvas.Update()
c1.Modified()
c1.Update()
c1.SaveAs(opt.weboutputdir+"/"+opt.outputfile+".png")

#draw the r histogram
hr = TH1D("hr", "", 150000, -10, 10)
tree.Draw("r>>hr", "trackedParamErr_r>0")
hr.Rebin( int(round( (hr.FindLastBinAbove(0)-hr.FindFirstBinAbove(0))/32 )) )
hr.GetXaxis().SetRangeUser( hr.GetMean()-5*hr.GetRMS(), hr.GetMean()+5*hr.GetRMS())
#hr.Fit("gaus")
#gaus = hr.GetListOfFunctions().FindObject("gaus")

hr.GetYaxis().SetTitle("toys")
hr.GetYaxis().SetTitleSize(0.06)
hr.GetYaxis().SetTitleOffset(0.7)

hr.GetXaxis().SetTitle("#sigma")
hr.GetXaxis().SetTitleSize(0.06)
hr.GetXaxis().SetTitleOffset(0.7)
hr.Draw()

pt = TPaveText(0.65, 0.55, 0.89, 0.87,"ndc")
pt.SetFillColor(0)
pt.SetLineColor(1)
pt.SetLineStyle(1)
        
t1 = pt.AddText("toys number: "+str(hr.GetEntries()))
t1.SetTextColor(1)
t1.SetTextSize( 0.04 )
        
#t_mu     = pt.AddText("#mu    = %s #pm %s" % ( str( format(gaus.GetParameter(1), ".2E") ), str( format(gaus.GetParError(1), ".2E"))      )    )
#t_sigma  = pt.AddText("#sigma = %s #pm %s" % ( str( format(gaus.GetParameter(2), ".2E") ), str( format(gaus.GetParError(2), ".2E"))      )    )
t_mean   = pt.AddText("mean   = %s "       % ( str( format(hr.GetMean()        , ".2E") ) ) )    
t_StdDev = pt.AddText("StdDev = %s "       % ( str( format(hr.GetStdDev()      , ".2E") ) ) )    
t_r      = pt.AddText("signal inj r = %s"  % ( str( format(float(opt.rExpected), ".2E") ) ) )

#t_mu.SetTextColor(1)
#t_mu.SetTextSize( 0.04 )
#t_sigma.SetTextColor(1)
#t_sigma.SetTextSize( 0.04 )
t_mean.SetTextColor(1)
t_mean.SetTextSize( 0.04 )
t_StdDev.SetTextColor(1)
t_StdDev.SetTextSize( 0.04 )
t_r.SetTextColor(1)
t_r.SetTextSize( 0.04 )

pt.Draw("SAME")

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c1, iPeriod, iPos)    #canvas.Update()
c1.Modified()
c1.Update()
c1.SaveAs(opt.weboutputdir+"/hr_"+opt.outputfile+".png")

# correct bias with experimental sigma
h_pull_sigmaexp = TH1D("h_pull_sigmaexp", "", 400, -4, 4)

sigma = hr.GetRMS()
print(sigma)
inputfile.Close()
#tree.Draw("(r-"+str(opt.rExpected)+")/"+str(sigma)+">>h_pull_sigmaexp", "trackedParamErr_r>0")
#h_pull_sigmaexp.GetXaxis().SetRangeUser(-1,1)
#h_pull_sigmaexp.Fit("gaus")
#gaus = h_pull_sigmaexp.GetListOfFunctions().FindObject("gaus")

#h_pull_sigmaexp.GetYaxis().SetTitle("toys")
#h_pull_sigmaexp.GetYaxis().SetTitleSize(0.06)
#h_pull_sigmaexp.GetYaxis().SetTitleOffset(0.7)
#
#h_pull_sigmaexp.GetXaxis().SetTitle("(#sigma-#sigma_{inj})/#sigma_{err}")
#h_pull_sigmaexp.GetXaxis().SetTitleSize(0.06)
#h_pull_sigmaexp.GetXaxis().SetTitleOffset(0.7)
#h_pull_sigmaexp.GetXaxis().SetRangeUser(-4,4)
#h_pull_sigmaexp.Draw()
#
#pt = TPaveText(0.65, 0.55, 0.89, 0.87,"ndc")
#pt.SetFillColor(0)
#pt.SetLineColor(1)
#pt.SetLineStyle(1)
#        
#t1 = pt.AddText("toys number: "+str(h_pull_sigmaexp.GetEntries()))
#t1.SetTextColor(1)
#t1.SetTextSize( 0.04 )
#
##t_mu     = pt.AddText("#mu    = %s #pm %s" % ( str(round(gaus.GetParameter(1),2)) , str(round(gaus.GetParError(1),2))      )    )
##t_sigma  = pt.AddText("#sigma = %s #pm %s" % ( str(round(gaus.GetParameter(2),2)) , str(round(gaus.GetParError(2),2))      )    )
#t_mean   = pt.AddText("mean   = %s " % ( str(round(h_pull_sigmaexp.GetMean(), 3))    ) )    
#t_StdDev = pt.AddText("StdDev = %s " % ( str(round(h_pull_sigmaexp.GetStdDev(),2))   ) )    
#t_r      = pt.AddText("signal inj r = %s"  % ( str( format(float(opt.rExpected), ".2E") ) ) )
#
##t_mu.SetTextColor(1)
##t_mu.SetTextSize(0.04)
##t_sigma.SetTextColor(1)
##t_sigma.SetTextSize( 0.04 )
#t_mean.SetTextColor(1)
#t_mean.SetTextSize( 0.04 )
#t_StdDev.SetTextColor(1)
#t_StdDev.SetTextSize( 0.04 )
#t_r.SetTextColor(1)
#t_r.SetTextSize( 0.04 )
#
#pt.Draw("SAME")
#
##draw the lumi text on the canvas
#CMS_lumi.CMS_lumi(c1, iPeriod, iPos)    #canvas.Update()
#c1.Modified()
#c1.Update()
#
#c1.SaveAs(opt.weboutputdir+"/"+opt.outputfile+".pdf")
#c1.SaveAs(opt.weboutputdir+"/"+opt.outputfile+"_sigmaexp.png")

