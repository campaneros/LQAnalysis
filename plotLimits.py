#! /usr/bin/env python

from os import *
import os
import sys
import optparse
import datetime
import subprocess
import io
import array

from array import array
from glob import glob
import ROOT as R
#from ROOT import R.TGraph, R.TFile, R.gROOT, R.TCanvas, R.gPad, R.TLegend, TColor

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi

R.gROOT.SetBatch(True)

usage = "usage: python plotLimits.py -i /afs/cern.ch/work/s/santanas/Releases/CMSSW_8_1_0_Trijet/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna/output/fit_data/limit/Res1ToRes2ToGluGlu -l 40926"

parser = optparse.OptionParser(usage)

parser.add_option("-i", "--inputdir", dest="inputdir",
                  help="input directory with limits")

parser.add_option("-l", "--lumi", dest="lumi",
                  help="input directory with limits")

parser.add_option("-s", dest="single", default=False, action = 'store_true',
                  help="single category limit")

(opt, args) = parser.parse_args()

if not opt.inputdir:   
    parser.error('input directory not provided')

if not opt.lumi:   
    parser.error('lumi not provided')

##################################################################################################

CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(float(opt.lumi)/1000.)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

var = {}
limit_2sigmaUp = {}
limit_1sigmaUp = {}
limit_central = {}
limit_1sigmaDown = {}
limit_2sigmaDown = {}
limit_observed = {}
sigma = {}

counter = 0
#file_list=next(os.walk(opt.inputdir))[2]

#file_list.sort()
cross = array('d')
file_list=[]
for path, subdirs, files in os.walk(opt.inputdir):
	for name in files:
		#print(os.path.join(path, name))
		file_list.append(os.path.join(path, name))	

file_list.sort()

print(file_list)

app=""

for fl in file_list:
        #print fl

        if "AsymptoticLimits" not in fl:
            continue
        
        #firstword_fl=fl.split("_")[0]        
        #if firstword_fl != "higgsCombine":
        #    continue
        filename=os.path.basename(fl)

        model = (filename.split(".AsymptoticLimits")[0]).split("_") [1]
        mass = ( filename.split(".AsymptoticLimits")[0]).split("_") [2]
        L = ( filename.split(".AsymptoticLimits")[0]).split("_") [3]
        cat = (filename.split(".AsymptoticLimits")[0]).split("_") [4]
        print(model, mass, L, cat)


        #fl_fullpath = opt.inputdir+"/"+fl
        print(fl)

        inputfile = R.TFile.Open(str(fl))
        #inputfile.ls()
        tree = inputfile.Get("limit")
 
        if app != cat:
            app=cat
            var[cat] = []
            limit_central[cat] = []
            limit_1sigmaUp[cat] = []
            limit_1sigmaDown[cat] = []
            limit_2sigmaUp[cat] = []
            limit_2sigmaDown[cat] = []
            limit_observed[cat] = []
            sigma[cat] = []
        if "0p1" in L:
            sigma[cat].append(9.51E-04)
            sigma[cat].append(5.16E-05)
            sigma[cat].append(5.73E-06)
        elif "1p0"in L:
            sigma[cat].append(9.53E-02)
            sigma[cat].append(5.27E-03)
            sigma[cat].append(6.09E-04)
        var[cat].append( float(mass) )
        for event in tree:
            limit = event.limit
            quantile = event.quantileExpected
            print(limit, quantile)

            if quantile>0.95 and quantile<1:
                limit_2sigmaUp[cat].append(limit)
            if quantile>0.80 and quantile<0.85:
                limit_1sigmaUp[cat].append(limit)
            if quantile==0.5:
                limit_central[cat].append(limit)
            if quantile>0.15 and quantile<0.17:
                limit_1sigmaDown[cat].append(limit)
            if quantile>0.02 and quantile<0.03:
                limit_2sigmaDown[cat].append(limit)
            if quantile == -1:
                limit_observed[cat].append(limit)

        inputfile.Close()
        counter+=1

print(var)
yellow={}
green={}
median={}
observed={}
cross = {}

for key in var:
    
    N = len(var[key])
    #print N
    yellow[key] = R.TGraph(2*N)    # yellow band
    green[key] = R.TGraph(2*N)     # green band
    median[key] = R.TGraph(N)      # median line
    observed[key] = R.TGraph(N)
    cross[key] =  R.TGraph(N)   # median line

    for i in range(N):
	#if opt.single and "all" in cat:
        #print var[key][i]
            yellow[key].SetPoint(    i,    float(var[key][i]), float(limit_2sigmaUp[key][i]) ) # + 2 sigma
            green[key].SetPoint(     i,    float(var[key][i]), float(limit_1sigmaUp[key][i]) ) # + 1 sigma
            median[key].SetPoint(    i,    float(var[key][i]), float(limit_central[key][i]) ) # median[key]
            green[key].SetPoint(  2*N-1-i, float(var[key][i]), float(limit_1sigmaDown[key][i]) ) # - 1 sigma
            yellow[key].SetPoint( 2*N-1-i, float(var[key][i]), float(limit_2sigmaDown[key][i]) ) # - 2 sigma
            observed[key].SetPoint(  i,    float(var[key][i]), float(limit_observed[key][i]) ) # observed
            cross[key].SetPoint(i, float(var[key][i]), float(sigma[key][i]))
	#elif opt.single and (not "all" in cat):
        #	median[key].SetPoint(    i,    float(var[key][i]), float(limit_central[key][i]) ) # median[key]
	#else:
        #	yellow[key].SetPoint(    i,    float(var[key][i]), float(limit_2sigmaUp[key][i]) ) # + 2 sigma
        #	green[key].SetPoint(     i,    float(var[key][i]), float(limit_1sigmaUp[key][i]) ) # + 1 sigma
        #	median[key].SetPoint(    i,    float(var[key][i]), float(limit_central[key][i]) ) # median[key]
        #	green[key].SetPoint(  2*N-1-i, float(var[key][i]), float(limit_1sigmaDown[key][i]) ) # - 1 sigma
        #	yellow[key].SetPoint( 2*N-1-i, float(var[key][i]), float(limit_2sigmaDown[key][i]) ) # - 2 sigma
        #	observed[key].SetPoint(  i,    float(var[key][i]), float(limit_observed[key][i]) ) # observed
		

    W = 800
    H = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    L = 0.04*W
    c = R.TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( L/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()

    frame = c.DrawFrame(1.4,0.001, 4.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% upper limit on #sigma / #sigma_{exp}")
    frame.GetXaxis().SetTitle("Resonance mass [GeV]")
    frame.SetMinimum(0)
    frame.SetMaximum(max(limit_2sigmaUp[key])*1.05)
    frame.GetXaxis().SetLimits(float(min(var[key])),float(max(var[key])))

    yellow[key].SetFillColor(R.kOrange)
    yellow[key].SetLineColor(R.kOrange)
    yellow[key].SetFillStyle(1001)
    yellow[key].Draw('F')
 
    green[key].SetFillColor(R.kGreen+1)
    green[key].SetLineColor(R.kGreen+1)
    green[key].SetFillStyle(1001)
    green[key].Draw('Fsame')
 
    median[key].SetLineColor(1)
    median[key].SetLineWidth(2)
    median[key].SetLineStyle(2)
    median[key].Draw('Lsame')

    observed[key].SetLineColor(1)
    observed[key].SetLineWidth(2)
    observed[key].SetLineStyle(1)
    observed[key].Draw('Lsame')
 
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)   
    R.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')
 
    x1 = 0.45
    x2 = x1 + 0.24
    y2 = 0.76
    y1 = 0.60
    legend = R.TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(observed[key], "Asymptotic CL_{s} observed",'L')
    legend.AddEntry(median[key], "Asymptotic CL_{s} expected",'L')
    legend.AddEntry(green[key], "#pm 1 std. deviation",'f')
    legend.AddEntry(yellow[key],"#pm 2 std. deviation",'f')
    legend.Draw()
    print(" ")
    canvasfilename_png = "UpperLimits_R"+str(key)+".png"
    canvasfilename_pdf = "UpperLimits_R"+str(key)+".pdf"
    c.SaveAs(canvasfilename_png)
    c.SaveAs(canvasfilename_pdf)
    c.Close()

    os.system("mv "+canvasfilename_png+" "+opt.inputdir)
    os.system("mv "+canvasfilename_pdf+" "+opt.inputdir)



W = 800
H = 600
T = 0.08*H
B = 0.12*H
L = 0.12*W
L = 0.04*W
c = R.TCanvas("c","c",100,100,W,H)
c.SetFillColor(0)
c.SetBorderMode(0)
c.SetFrameFillStyle(0)
c.SetFrameBorderMode(0)
c.SetLeftMargin( L/W )
c.SetRightMargin( L/W )
c.SetTopMargin( T/H )
c.SetBottomMargin( B/H )
c.SetTickx(0)
c.SetTicky(0)
c.SetGrid()
c.cd()


for count,key in enumerate(sorted(var)):
    print(key)
    if "all" in key:
        frame = c.DrawFrame(1.4,0.001, 4.1, 10)
        frame.GetYaxis().CenterTitle()
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleOffset(0.9)
        frame.GetXaxis().SetNdivisions(508)
        frame.GetYaxis().CenterTitle(True)
        frame.GetYaxis().SetTitle("95% upper limit on #sigma / #sigma_{exp}")
        frame.GetXaxis().SetTitle("Resonance mass [GeV]")
        frame.SetMinimum(min(min(limit_2sigmaUp[key]),min(sigma[key]))*0.10)
        frame.SetMaximum(max(max(limit_2sigmaUp[key]),max(sigma[key]))*10)
        frame.GetXaxis().SetLimits(float(min(var[key])),float(max(var[key])))
        c.SetLogy(1)        
        yellow[key].SetFillColor(R.kOrange)
        yellow[key].SetLineColor(R.kOrange)
        yellow[key].SetFillStyle(1001)
        yellow[key].Draw('F')

        green[key].SetFillColor(R.kGreen+1)
        green[key].SetLineColor(R.kGreen+1)
        green[key].SetFillStyle(1001)
        green[key].Draw('Fsame')

        median[key].SetLineColor(1)
        median[key].SetLineWidth(2)
        median[key].SetLineStyle(2)
        median[key].Draw('Lsame')      
        observed[key].SetLineColor(1)
        observed[key].SetLineWidth(2)
        observed[key].SetLineStyle(1)
        observed[key].Draw('Lsame')

        cross[key].SetLineColor(R.kBlue)
        cross[key].SetLineWidth(2)
        cross[key].SetLineStyle(1)
        cross[key].Draw('Lsame')


        CMS_lumi.CMS_lumi(c, iPeriod, iPos)   
        R.gPad.SetTicks(1,1)
        frame.Draw('sameaxis')
        x1 = 0.40
        x2 = x1 + 0.24
        y2 = 0.90
        y1 = 0.70
        legend = R.TLegend(x1,y1,x2,y2)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.041)
        legend.SetTextFont(42)
        legend.AddEntry(observed[key], "Asymptotic CL_{s} observed",'L')
        legend.AddEntry(median[key], "Asymptotic CL_{s} expected %s"%key,'L')
        legend.AddEntry(green[key], "#pm 1 std. deviation",'f')
        legend.AddEntry(yellow[key],"#pm 2 std. deviation",'f')
        legend.AddEntry(cross[key],"theoretical cross section",'L')
    else:
        median[key].SetLineColor(count*2)
        median[key].SetLineWidth(2)
        median[key].SetLineStyle(2)
        median[key].Draw('Lsame')
        legend.AddEntry(median[key], " Median asymptotic expected %s"%(key.strip("category")),'L')

    legend.Draw()
for ext in ['png','pdf']:
	canvasfilename = "UpperLimits_all_cat."+ext
#canvasfilename_pdf = "UpperLimits_all_cat.pdf"
	c.SaveAs(canvasfilename)
	os.system("mv "+canvasfilename+" "+opt.inputdir)
#c.SaveAs(canvasfilename_pdf)
c.Close()

#os.system("mv "+canvasfilename_pdf+" "+opt.inputdir)

	
print("Output in: "+opt.inputdir)

