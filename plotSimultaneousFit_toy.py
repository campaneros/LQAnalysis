#! /usr/bin/env python

from os import *
import os
import sys
import optparse
import datetime
import subprocess
import io
import ROOT
import re

from array import array
from glob import glob
from ROOT import *
from ROOT import AddressOf
from ROOT import RooFit

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../")
import tdrstyle
import CMS_lumi

gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

usage = "usage: python plotSimFit.py -t toysfile/workspacefile -f fitfile -n 1 -c catdir -o outputdir -w weboutputdir -F fitFunction --draw_limit_(exp/obs) limit_dir"

parser = optparse.OptionParser(usage)




parser.add_option("-t", "--toysfile", dest="toysfile", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_LQumu_M3000_L1p0_t_100_syst0_seed1123456/higgsCombine_toys100_expectSignal0.0.GenerateOnly.mH120.123456.root",
                  help="input file with fitted toys")

parser.add_option("-f", "--fitfile", dest="fitfile", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_LQumu_M3000_L1p0_t_100_syst0_seed1123456/higgsCombine_toys100_expectSignal0.0_gen.MultiDimFit.mH120.123456.root",
                  help="input file with tree of post-fit parameters.")

parser.add_option("-c", "--catdir", dest="catdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M1000_L0p1/categories/",
                  help="name of directory containing categories dirs")
                  
parser.add_option("-o", "--outputdir", dest="outputdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/plotSimFit_std_4par/LQumu_M1000_L0p1_bis",
                  help="name of the output directory")

parser.add_option("-b", "--outputfilename", dest="outputfile", default="expect_signal",
                  help="name of the output file")

parser.add_option("-w", "--weboutputdir", dest="weboutputdir", default="/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_plot/plotSimFit_std_4par/LQumu_M1000_L0p1_bis",
                  help="name of the web output directory")

parser.add_option("-s", "--signaname", dest="signaname", default="LQumu_M1000_L0p1",
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

parser.add_option("--nc", dest="ncat", default="all",
                  help="fcategories to be run")


(opt, args) = parser.parse_args()

if not opt.toysfile:   
    parser.error('input toy file not provided')

if not opt.fitfile:   
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
varBins_all = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058,1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869,5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7150, 7500, 7850, 8250, 8650, 8999, 9500, 9999]
NvarBins_all = len(varBins_all)-1

fitFunction = opt.fitFunction

#Load toy data histogram
nToy = int(opt.nToy)
toysfilenamePath = opt.toysfile
toysfilename     = toysfilenamePath.split("/")[-1]
toysfile = TFile.Open(toysfilenamePath)
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
workspacefile = TFile.Open(workspacePath)

workspacename = "w"
workspace = RooWorkspace()

workspace=workspacefile.Get(workspacename)

toy = toysfile.Get("toys/toy_1")
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

#Create output rootfile with chi2
indexcat = array('i', [0 ])
Chi2     = array('d', [0.])
Ndof     = array('d', [0.])
Mass  = array('d', [0.])
L    = array('d', [0.])

GlobalChi2 = array('d', [0. ])
GlobalNdof = array('d', [-1.]) #subtract r parameter
CombineGoF = array('d', [0. ])
CombineDof = array('d', [-1.]) #subtract r parameter

chi2file = TFile.Open(opt.outputdir+"/test_statistics.root", "RECREATE")

chi2tree = TTree("chi2tree","chi2tree")
icatbranch = chi2tree.Branch("icat", indexcat, "icat/I")
chi2branch = chi2tree.Branch("chi2", Chi2, "chi2/D")
ndofbranch = chi2tree.Branch("ndof", Ndof, "ndof/D")
chi2branch = chi2tree.Branch("Mass", Mass, "Mass/D")
ndofbranch = chi2tree.Branch("L", L, "L/D")

globchi2tree = TTree("globchi2tree","globchi2tree")
globchi2branch = globchi2tree.Branch("globchi2"  , GlobalChi2, "globchi2/D"  )
globndofbranch = globchi2tree.Branch("globndof"  , GlobalNdof, "globndof/D"  )
combiGoFbranch = globchi2tree.Branch("CombineGoF", CombineGoF, "CombineGoF/D")
combiDoFbranch = globchi2tree.Branch("CombineDof", CombineDof, "CombineDof/D")

outputrootfile = [None] * len(categoriesList)
print(categoriesList)

for number in range(nToy):
    toy = toysfile.Get("toys/toy_"+str(number+1))
    #Make plot of data + fit for each category
    #var_min_limit = M1*0.85

    for icat,cat in enumerate(categoriesList):

        if "workspace_"  in cat:
            continue
        elif "gen_" in cat:
    #           print("maremma boia")
    #           print(cat)
                continue
        elif opt.ncat != "all" and opt.ncat not in cat:
            continue
        print(cat)
        outputrootfile[icat] = TFile.Open(opt.outputdir+"/"+cat+".root","RECREATE")
        outputrootfile[icat].cd()
        indexcat[0]+=1

        #Define mjj category range
        #print cat
        splitname = cat.split("_")
        Mass[0] = float(splitname[2].strip("M"))
        L[0] = float(splitname[3].strip("L").replace("p","."))
        #L[0] = float(L.replace("p","."))

        app=(cat.split("_"))[4]
        sign_=(cat.strip("datacard_"))
        print(app)
        print(sign_)
        #x = RooRealVar("m_muj_ak4_category1Muon","m_muj_ak4_category1Muon",41,453,3854)
        #x = workspace.var("m_muj_ak4_"+app)

        t= toy.get()
        toyvar= t.find("m_muj_ak4_"+app)
        print(toyvar)
        test = toy.binnedClone("test","test")
        #print(test.binVolume())
        #print(test.Print())
        var_min_set = toyvar.getMin()
        var_max_set = toyvar.getMax()

        print(var_min_set)
        print(var_max_set)

        #var_max_set=4000    
        #var_min_set=500   

        varBins = []
        for border in varBins_all:
            if border<var_min_set:
                continue
            if border>var_max_set:
                continue
            else:
                varBins.append(float(border))
        if var_min_set < varBins[0]:
            varBins.insert(0,var_min_set)
        if var_max_set > varBins[-1]:
            varBins.append(var_max_set)

        a = array('d', varBins)
        #b = 1
        #NvarBins = 41
        b = toyvar.getBinning()
        NvarBins = toyvar.getBins()
        print(NvarBins)


        x = workspace.var("m_muj_ak4_"+app)
        
        #Get toy data RooDataHist from file
        #roohist_data = workspace.data("RooDataHist_"+app)
        hData = test.createHistogram("hData_"+app, toyvar, RooFit.Binning(b))
        hData_norm = hData.Clone("hData_norm_"+app)

        #Get ParametricShapeBin pdf and parameters from workspace and extend with norm
        BkgFit    = workspace.pdf("ParametricBkgPdf_"+app)

        if "std" in fitFunction:
            p1_var    = workspace.var("p1_std"+app) 
            p2_var    = workspace.var("p2_std"+app) 
            p3_var    = workspace.var("p3_std"+app) 
        elif "UA2" in fitFunction:
            p1_var    = workspace.var("p1_UA2"+app) 
            p2_var    = workspace.var("p2_UA2"+app) 
        elif "modExp" in fitFunction:
            p1_var    = workspace.var("p1_modExp"+app) 
            p2_var    = workspace.var("p2_modExp"+app) 

        norm_var  = workspace.var("ParametricBkgPdf_"+app+"_norm")
        ExtBkgFit = RooExtendPdf("ExtBkgPdf_"+app, "ExtBkgPdf_"+app, BkgFit, norm_var)

        ## Get signal shape and floating parameters (only sys floating in the fit=
        #if signalname == "":
        #    signalname = opt.appdir.split("/")[-2]
        #signalstring = signalname+"_"+app
        signalPdf = workspace.pdf("ParametricSignalPdf_"+sign_)

        #JES_sys = workspace.var("JES_uncertainty")
        #JER_sys = workspace.var("JER_uncertainty")
        #meanShape_mu_err = workspace.var("meanShape_mu_err_"+signalstring)
        #meanShape_sigma_err = workspace.var("meanShape_sigma_err_"+signalstring)

        nsig = workspace.var("ParametricSignalPdf_"+sign_+"_norm")
        print("maremma maiala",sign_)
        signalExtendPdf = RooExtendPdf("ParametricSignalExtPdf_"+sign_, "ParametricSignalExtPdf_"+sign_, signalPdf, nsig)



        #Load post-simultaneous fit parameters
        fitfile = TFile.Open(opt.fitfile)
        fitTree = fitfile.Get("limit")

        # Bkg parameters
        p0_postfit = array('f', [0.])
        p1_postfit = array('f', [0.])
        p2_postfit = array('f', [0.])
        p3_postfit = array('f', [0.])

        fitTree.SetBranchAddress("trackedParam_shapeBkg_bkg_"+app+"__norm", p0_postfit)

        if "std" in fitFunction:
            fitTree.SetBranchAddress("trackedParam_p1_std"+app, p1_postfit)
            fitTree.SetBranchAddress("trackedParam_p2_std"+app, p2_postfit)
            fitTree.SetBranchAddress("trackedParam_p3_std"+app, p3_postfit)
        if "UA2" in fitFunction:
            fitTree.SetBranchAddress("trackedParam_p1_UA2"+app, p1_postfit)
            fitTree.SetBranchAddress("trackedParam_p2_UA2"+app, p2_postfit)
        if "modExp" in fitFunction:
            fitTree.SetBranchAddress("trackedParam_p1_modExp"+app, p1_postfit)
            fitTree.SetBranchAddress("trackedParam_p2_modExp"+app, p2_postfit)

        fitTree.GetEntry(number)

        norm_var.setVal(p0_postfit[0])
        p1_var.setVal(p1_postfit[0])
        p2_var.setVal(p2_postfit[0])
        p3_var.setVal(p3_postfit[0])


        print("P1 post fit ",p1_var.getVal())
        print("P0 post fit ",norm_var.getVal())
        print("P2 post fit ",p2_var.getVal())
        print("P3 post fit ",p3_var.getVal())
        # Signal parameters
        #JES_sys_postfit = array('f', [0.])
        #JER_sys_postfit = array('f', [0.])
        #meanShape_mu_err_postfit = array('f', [0.])
        #meanShape_sigma_err_postfit = array('f', [0.])

        #fitTree.SetBranchAddress("trackedParam_JES_uncertainty", JES_sys_postfit)
        #fitTree.SetBranchAddress("trackedParam_JER_uncertainty", JER_sys_postfit)
        # if fitTree.GetListOfBranches().FindObject("trackedParam_meanShape_mu_err_"+signalstring): 
        #     fitTree.SetBranchAddress("trackedParam_meanShape_mu_err_"+signalstring, meanShape_mu_err_postfit)
        # else:
        #     meanShape_mu_err_postfit[0] = 1
        #fitTree.SetBranchAddress("trackedParam_meanShape_sigma_err_"+signalstring, meanShape_sigma_err_postfit)

        fitTree.GetEntry(0)

        #JES_sys.setVal(JES_sys_postfit[0])
        #JER_sys.setVal(JER_sys_postfit[0])
        #meanShape_mu_err.setVal(meanShape_mu_err_postfit[0])
        #meanShape_sigma_err.setVal(meanShape_sigma_err_postfit[0])

        # Signal r = fit cross section
        r = fitTree.r
        #r=0.0001586 
        #r=0.0026406
        err_r = fitTree.trackedParamErr_r
        #err_r=0 

        # Signal r = expected limit on cross section
        if opt.draw_limit_exp != False:
            limitdir  = opt.draw_limit_exp
            limitfile = TFile.Open(limitdir+"/higgsCombine_Res1ToRes2ToGluGlu_"+str(M1)+"_"+str(R)+".AsymptoticLimits.mH120.root")
            limittree = limitfile.Get("limit")
            limittree.GetEntry(2)
            r = limittree.limit

        # Signal r = observed limit on cross section
        if opt.draw_limit_obs != False:
            limitdir  = opt.draw_limit_obs
            limitfile = TFile.Open(limitdir+"/higgsCombine_Res1ToRes2ToGluGlu_"+str(M1)+"_"+str(R)+".AsymptoticLimits.mH120.root")
            limittree = limitfile.Get("limit")
            limittree.GetEntry(5)
            r = limittree.limit

        fitfile.Close()

        #Make histogram from signal and background shape
        hSignal = TH1D("hSignal_"+app,"",NvarBins,a)
        hBkg    = TH1D("hBkg_"+app,"",NvarBins,a)

        x.setRange("global_range", a[0], a[NvarBins])
        print(a[0], a[NvarBins])
        signal_pdfIntrinsicNorm = signalPdf.createIntegral(ROOT.RooArgSet(x),RooFit.NormSet(ROOT.RooArgSet(x)),RooFit.Range("global_range"))
        bkg_pdfIntrinsicNorm    = BkgFit.createIntegral(ROOT.RooArgSet(x),RooFit.NormSet(ROOT.RooArgSet(x)),RooFit.Range("global_range"))
        print(signal_pdfIntrinsicNorm.getVal(), nsig.getVal())


        Sum=0

        for ibin, bincenter in enumerate(a):
            ibin+=1
            bin_low   = hSignal.GetBinLowEdge(ibin)
            bin_width = hSignal.GetBinWidth(ibin)
            bin_up    = bin_low + bin_width
            
        #x.setMin(bin_low)
        #x.setMax(bin_up)
            
            x.setRange("toy_"+str(ibin),bin_low,bin_up)

            #print(x.setRange("binrange",bin_low,bin_up))
            
    #	print("bin range ",bin_low,"  ", bin_up)
    #	print("X range ",x.getMin(),"  ", x.getMax())
            signal_pdfIntegral = signalPdf.createIntegral(ROOT.RooArgSet(x),RooFit.NormSet(ROOT.RooArgSet(x)),RooFit.Range("toy_"+str(ibin)))
            signal_pdfIntegral_norm = signal_pdfIntegral.getVal()/signal_pdfIntrinsicNorm.getVal()
            bin_signalEvents = signal_pdfIntegral_norm*nsig.getVal()*r
            bin_signalEvents_norm = signal_pdfIntegral_norm*nsig.getVal()*r/bin_width
    #	print("##########################################################")
    #	print("Fraction of signal pdf of  bin area", signal_pdfIntegral_norm)
    #	print("Event expected", nsig.getVal())
    #	#print("Signal strenght", r)
    #	#print("Events per bin normalised per bin width ",bin_signalEvents_norm)
    #	#print("Events per bin ",bin_signalEvents)
    #	#print("Bin width ", bin_width)
    #	#print("norm", signal_pdfIntegral_norm)
    #	print("integral", signal_pdfIntegral.getVal())
    #	#print("signal_pdfIntrinsicNorm.getVal()", signal_pdfIntrinsicNorm.getVal())
    #	##print("signal_pdfIntegral.getVal()", signal_pdfIntegral.getVal())
    #	print("##########################################################")
            #
            Sum+=bin_signalEvents

        # 
            hSignal.SetBinContent(ibin,bin_signalEvents_norm)

            bkg_pdfIntegral = BkgFit.createIntegral(ROOT.RooArgSet(x),RooFit.NormSet(ROOT.RooArgSet(x)),RooFit.Range("toy_"+str(ibin)))
            bkg_pdfIntegral_norm = bkg_pdfIntegral.getVal()/bkg_pdfIntrinsicNorm.getVal()
            bin_bkgEvents = bkg_pdfIntegral_norm*norm_var.getVal()
            bin_bkgEvents_norm = bkg_pdfIntegral_norm*norm_var.getVal()/bin_width
        
        #print("background ", bkg_pdfIntegral, " ", bin_bkgEvents_norm)       
    
            #hBkg.SetBinContent(ibin,bin_bkgEvents)
            hBkg.SetBinContent(ibin,bin_bkgEvents_norm)
            
            bin_data = hData.GetBinContent(ibin)
            #print(bin_data)
            hData_norm.SetBinContent(ibin, bin_data/bin_width)
            hData_norm.SetBinError(ibin, TMath.Sqrt(bin_data)/bin_width)
        


        print("#################################################################")
        print("#################################################################")
        print("######          #######   ###     ##   ####              ########")
        print("######          ##        ## ##   ##   ##  ##            ########")
        print("######          #######   ##  ##  ##   ##   ##           ########")
        print("######          ##        ##   ## ##   ##  ##            ########")
        print("######          #######   ##     ###   ####              ########")
        print("#################################################################")
        print("#################################################################")
        print("Sum of all events", Sum)
        print("Expected Events * Signal strenght ", nsig.getVal()*r )
        print("Signal strenght", r)

        #Pull histo2 and chi square computation
        Ndof[0] = -3
        Chi2[0] = 0
        GlobalNdof[0] -= 3
        CombineDof[0] -= 3
        redChi2 = 0.
        hist_pull = TH1D("pull_"+app, "", NvarBins, a)
        hist_pull_signal = TH1D("pull_signal_"+app, "", NvarBins, a)

        frame_Ymax = 0
        #frame_Ymin =40
        frame_Ymin = hData.GetBinContent(hData.GetNbinsX()-1)
        #print(hData.GetBinContent(hData.GetNbinsX()-1))
        for ibin in range(hData.GetNbinsX()):
            ibin  +=1
            bin_width = float(hBkg.GetBinWidth(ibin))
            #data  = float(hData.GetBinContent(ibin))
            data  = float(hData_norm.GetBinContent(ibin))
            bkg   = float(hBkg.GetBinContent(ibin) )
            bkg_err = TMath.Sqrt(bkg*bin_width)/bin_width
            sign  = float(hSignal.GetBinContent(ibin))

        ##print(data)
            if data > frame_Ymax:
                frame_Ymax = data
            #print("Eureka  ",frame_Ymax)
            if data < frame_Ymin:
                frame_Ymin = data

            CombineDof[0] += 1
            if data!=0 and bkg!=0:
                pull = (data-bkg)/bkg_err
                hist_pull.SetBinContent(ibin,pull)
                CombineGoF[0] += 2*(  (bkg+sign*r)-data + data*TMath.Log( data/(bkg+sign*r) )  )
            #print data, bkg, sign*r, 2*(  (bkg+sign*r)-data + data*TMath.Log( data/(bkg+sign*r) )  ), GlobalNdof[0], CombineGoF[0]
            else:
                CombineGoF[0] += 2*( (bkg+sign*r)-data )
            #print data, bkg, sign*r, 2*( (bkg+sign*r)-data ), GlobalNdof[0], CombineGoF[0]

            if bkg_err !=0:
                pull_signal = (sign)/bkg_err
            else:
                pull_signal = 0
            hist_pull_signal.SetBinContent(ibin,pull_signal)

            if data>10/bin_width:
                Chi2[0] += pull*pull
                Ndof[0] += 1.0
                GlobalChi2[0] += pull*pull
                GlobalNdof[0] += 1.0

        #Fill chi2 tree
        chi2tree.Fill()
        
        if Ndof[0] <  0:
            Ndof[0] = 0
        if Ndof[0] != 0:
            redChi2 = Chi2[0]/Ndof[0]

        ## Create canvas
        canvas = TCanvas("canvas_"+app, "canvas_"+app, 200, 10, 700, 500 )

        fPads1 = TPad("pad1_"+app, "pad1_"+app, 0.00, 0.28, 0.99, 0.99)
        fPads2 = TPad("pad2_"+app, "pad2_"+app, 0.00, 0.00, 0.99, 0.345)
        #
        fPads1.SetFillColor(0)
        fPads1.SetLineColor(0)
        fPads1.SetTopMargin(0.1)
        fPads1.SetBottomMargin(0.1)
        fPads1.SetRightMargin(0.1)
        fPads1.SetLeftMargin(0.13)
        fPads1.SetLogy()
        #
        fPads2.SetFillColor(0)
        fPads2.SetLineColor(0)
        fPads2.SetGridx()
        fPads2.SetGridy()
        fPads2.SetTopMargin(0.05)
        fPads2.SetBottomMargin(0.3)
        fPads2.SetRightMargin(0.1)
        fPads2.SetLeftMargin(0.13)
        #
        fPads1.Draw()
        fPads2.Draw()

        fPads1.cd()
        ## Pad1
        frame = x.frame()

        #frame_Ymax = hData.GetMaximum(1)
        print(frame_Ymax)
        print(frame_Ymin)
        if frame_Ymin == 0:
            frame.GetYaxis().SetRangeUser(0.0005, frame_Ymax*6)
        else:
            frame.GetYaxis().SetRangeUser(frame_Ymin*0.5, frame_Ymax*6)
            #frame.GetYaxis().SetRangeUser(0.0005, frame_Ymax*6)


        #frame.GetYaxis().SetRangeUser(0.0005, 500000)
            
        frame.GetYaxis().SetTitle("#events / bin width (GeV)")
        frame.GetYaxis().SetTitleSize(0.06)
        frame.GetYaxis().SetTitleOffset(0.8)
        frame.GetYaxis().SetLabelSize(0.04)

        frame.SetTitle(" ")
            
        frame.Draw()
        #test.plotOn(frame)

        hBkg.SetLineColor(kRed)
        hBkg.SetLineWidth(2)
        hBkg.Draw("same")

        hSignal.SetLineColor(kBlue)
        hSignal.SetLineStyle(kDashed)
        hSignal.SetLineWidth(2)
        hSignal.Draw("same")

        hData.SetLineColor(kBlack)
        hData.SetMarkerStyle(8)
        #hData.Draw("same")
        gStyle.SetErrorX(0.0001)
        hData_norm.SetLineColor(kBlack)
        hData_norm.SetMarkerStyle(8)
        hData_norm.Draw("E1 SAME")

        #line = TLine(var_min_set, 0.05, var_min_set, frame_Ymax*6)
        #line.SetLineColor(kRed)
        #line.SetLineWidth(2)
        #line.Draw("same")

        canvas.Modified()
        canvas.Update()

        #draw the lumi text on the canvas
        CMS_lumi.CMS_lumi(fPads1, iPeriod, iPos)    
        canvas.Modified()
        canvas.Update()    

        # Pad2
        fPads2.cd()
        hist_pull.GetXaxis().SetTitle("m_{\muj} [GeV]")
        hist_pull.GetYaxis().SetTitle("#frac{Bin - Fit}{Uncertainty} ")
        hist_pull.SetTitle("")
        hist_pull.SetMinimum(-4)    
        hist_pull.SetMaximum(4)
        hist_pull.SetLineColor(2)
        hist_pull.SetFillColor(2)
        hist_pull.SetMarkerStyle(20)
        hist_pull.SetMarkerColor(1)
        hist_pull.SetStats(0)
        hist_pull.GetYaxis().SetNdivisions(405, kTRUE)
        hist_pull.GetXaxis().SetTitleSize(0.16)
        hist_pull.GetXaxis().SetLabelSize(0.13)
        hist_pull.GetXaxis().SetTitleOffset(0.83)
        hist_pull.GetYaxis().SetTitleSize(0.12)
        hist_pull.GetYaxis().SetLabelSize(0.11)
        hist_pull.GetYaxis().SetTitleOffset(0.35)
        hist_pull.Draw("hist")
        hist_pull_signal.SetLineColor(kBlue)
        hist_pull_signal.SetLineStyle(kDashed)
        hist_pull_signal.SetLineWidth(2)
        hist_pull_signal.Draw("same")


        ## Legend
        fPads1.cd()
        legend = TLegend(0.64, 0.65, 0.87, 0.87)
        legend.SetLineColor(0)
        legend.SetLineWidth(0)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)

        #RooPlot objects
        if(opt.fitData):
            legend.AddEntry(hData_norm, "Data", "p")
        else:
            legend.AddEntry(hData_norm, "Toy MC", "p")
        legend.AddEntry(hBkg, "Fit", "l")
        #if r < 0.001:
        #    legend.AddEntry(hSignal, "signal (fb) r="+str(round(r*1000,10)), "l")
        #else:
        legend.AddEntry(hSignal, "signal (pb) r="+str(round(r,6)) +"#pm  "+str(round(err_r,5)), "l")
        legend.Draw()
        
        ## Plot fit results
        pt = TPaveText(0.44, 0.58, 0.54, 0.87,"ndc")
        pt.SetFillColor(0)

        if nToy>0:
            t1 = pt.AddText("toy "+str(nToy))
            t1.SetTextColor(1)
            t1.SetTextSize( 0.04 )

        Chi2Text = "#chi^{2} / ndf (N_{bin}>10) = "+str(round(redChi2,2))
        tchi2 = pt.AddText(Chi2Text)
        tchi2.SetTextColor(1)
        tchi2.SetTextSize( 0.04 )
        tnorm = pt.AddText("nbkg "+str( round(norm_var.getVal()) )+" +/- "+str( round(norm_var.getError()) ))
        tnorm.SetTextColor(1)
        tnorm.SetTextSize( 0.04 )
        tp1 = pt.AddText("p1 "+str( round(p1_var.getVal(),2) )+" +/- "+str( round(p1_var.getError(),2) ))
        tp1.SetTextColor(1)
        tp1.SetTextSize( 0.04 )
        tp2 = pt.AddText("p2 "+str( round(p2_var.getVal(),2) )+" +/- "+str( round(p2_var.getError(),2) ))
        tp2.SetTextColor(1)
        tp2.SetTextSize( 0.04 )
        tp3 = pt.AddText("p3 "+str( round(p3_var.getVal(),2) )+" +/- "+str( round(p3_var.getError(),2) ))
        tp3.SetTextColor(1)
        tp3.SetTextSize( 0.04 )
        
        pt.Draw("same")
        canvas.Update()

        canvas.Update()

        #canvas.SaveAs(opt.outputdir+"/canvas_toy_"+str(nToy)+"_"+app+".pdf")
        canvas.SaveAs(opt.outputdir+"/canvas_toy_"+str(nToy)+"_"+app+"_"+opt.outputfile+"toyesnumber"+str(number)+".png")

        #Write to file
        outputrootfile[icat].cd()
        hBkg.Write()
        hSignal.Write()
        #hData.Write()
        #hData_norm.Write()
        #hist_pull_signal.Write()
        outputrootfile[icat].Close()

    #Fill Global chi2 tree
    if GlobalNdof[0]>0:
        GlobalReducedChi2 = GlobalChi2[0]/GlobalNdof[0]
    else:
        GlobalNdof[0] = 0
    if CombineDof[0] < 0:
        CombineDof[0] = 0
    globchi2tree.Fill()

    ## Save Chi2 in rootfile
    
  
    #Save on web dir
    #xscriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
    #os.system("mkdir -p "+ opt.weboutputdir)
    #os.system("cp "+opt.outputdir+"/*.p* "+opt.weboutputdir)
    #os.system("cp "+scriptsPath+"/index.php "+opt.weboutputdir)
    #os.system("cp "+scriptsPath+"/index.php "+opt.weboutputdir+"/../")
    #os.system("cp "+scriptsPath+"/index.php "+opt.weboutputdir+"/../../")

globchi2tree.Write()

chi2file.Write()
chi2file.Close()