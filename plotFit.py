
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
import create_workspaces_and_datacards_utils as cwd_utils





gROOT.SetBatch(True)
gErrorIgnoreLevel = kFatal
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)




signalInputfilename = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output/signals.txt"



c=TCanvas("t","t",900,600)

signalInputfile = io.open(signalInputfilename, "r")
for iline, line in enumerate(signalInputfile):
	line = line.rstrip('\n')
       	splitline = line.split(" ")
	#print(splitline)
	modell = splitline[0]
        category = splitline[1]
        M1 = splitline[2]
        L = splitline[3]
        Nsig = splitline[4]
        e_Nsig = splitline[5]
        mu = splitline[6]
        e_mu = splitline[7]
        ssigma = splitline[8]
        e_ssigma = splitline[9]
        a1 = splitline[10]
        e_a1 = splitline[11]
        n1 = splitline[12]
        e_n1 = splitline[13]
        a2 = splitline[14]
        e_a2 = splitline[15]
        n2 = splitline[16]
        e_n2 = splitline[17]
        
        #L=L.replace("p",".")
        print(modell)	
        if not category == cat:
            continue
        
        print("current signal = ", modell, category, M1, L, Nsig, e_Nsig, mu, e_mu, ssigma, e_ssigma, a1, e_a1, n1, e_n1, a2, e_a2, n2, e_n2)

		
        signalString = modell+"_"+category#+"_"+"M"+str(int(float(M1)))+"_R0p"+str(int(float(L)*10))+"_"+category
	print("signal_string   "+signalString)        	


        ## Signal pdf
        mean = RooRealVar("mean_"+signalString,"mean_"+signalString,float(mu)) 
        width = RooRealVar("width_"+signalString,"width_"+signalString,float(ssigma)) 
        alpha1 = RooRealVar("alpha1_"+signalString,"alpha1_"+signalString,float(a1))
        n1 = RooRealVar("n1_"+signalString,"n1_"+signalString,float(n1))
        alpha2 = RooRealVar("alpha2_"+signalString,"alpha2_"+signalString,float(a2))
        n2 = RooRealVar("n2_"+signalString,"n2_"+signalString,float(n2))
        nsig = RooRealVar("ParametricSignalPdf_"+signalString+"_norm","ParametricSignalPdf_"+signalString+"_norm",float(Nsig),0,100000000)

        var_tmp = RooRealVar("var_tmp","var_tmp",500, 5000)
        var_tmp.setRange("maximum_range", 500, 5000)  # create range to integrate over
        signalPdf_tmp = RooDoubleCBFast("CB_tmp", "CB_tmp", var_tmp, mean, width, alpha1, n1, alpha2, n2)            
        intrinsicNorm = signalPdf_tmp.createIntegral(RooArgSet(var_tmp), RooFit.NormSet(RooArgSet(var_tmp)), RooFit.Range("maximum_range")) 
        var_tmp.setRange("fit_range", var_min_set, var_max_set)  # create range to integrate over
        signal_integral = signalPdf_tmp.createIntegral(RooArgSet(var_tmp), RooFit.NormSet(RooArgSet(var_tmp)), RooFit.Range("fit_range"))
        nsig.setVal( float(Nsig)*signal_integral.getVal()/intrinsicNorm.getVal() )
	print(str(nsig.getVal()))

        
        signalPdf = RooDoubleCBFast("signalPdf_"+signalString, "signalPdf_"+signalString, var[icat], mean, width, alpha1, n1, alpha2, n2)

	
        ParametricSignalPdf_ = RooParametricShapeBinPdf("ParametricSignalPdf_"+signalString, "ParametricSignalPdf_"+signalString, signalPdf, var[icat], RooArgList(mean, width, alpha1, n1, alpha2, n2), th1_rebin[icat])



		
