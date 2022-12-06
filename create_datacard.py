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
import create_workspaces_and_datacards_utils as cwd_utils

fitFunction_name = "std_4par"
fitparam = []
#gROOT.LoadMacro(os.path.dirname(os.path.abspath(__file__))+"/../../src/libCpp/RooDoubleCBFast.cc+")
gROOT.SetBatch(True)
gErrorIgnoreLevel = kFatal
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

inputdir = "/data/mcampana/CMS/CMSSW_10_6_28_LQAna/src/RootTreeAnalyzer/Fit_Signal/category"
filenameInput = "test_h1_mmuj_ak4.root"
subDirList = next(os.walk(inputdir))[1]
print(subDirList)

## Output directories  
outputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC"
outputdirdatacards = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards"
weboutputdir = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_plot"
os.system("mkdir -p "+outputdir)
os.system("mkdir -p "+outputdirdatacards)
os.system(u"rm -f "+outputdirdatacards+"/*")
os.system("mkdir -p "+weboutputdir)
scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
os.system("cp "+scriptsPath+"/index.php "+weboutputdir)
filenameworkspace = "workspace_gen.root"
filenameworkspacePath = outputdir+"/"+filenameworkspace
workspaceName = "w"








## Histogram
#if opt.fitData:
#    generateToy = -1
#    histoname = "h1_mmuj_ak4__DATA"
#else:
generateToy = 1
histoname = "htot_bkg_MC_m_muj_ak4"

statMultiplier = 1
varname = "m_muj_ak4"
vartitle = "m_{\muj}_ak4 [GeV]"

## Signal input
signalInputfilename = "/data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output/signals.txt"


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

### ML Fit with combine ###
#Data no signal
#combine -M MultiDimFit datacards/datacard_gg_M3000_R0p1.txt --rMin -5 --rMax 5 -S 0 --algo singles --cl=0.68 --saveFitResult --trackParameters "rgx{.*}" -t 0 --name _santanas_ --keepFailures (--robustHesse 1)

#Toy no signal
#combine -M MultiDimFit datacards/datacard_gg_M3000_R0p1.txt --rMin -5 --rMax 5 -S 0 --algo singles --cl=0.68 --saveFitResult --saveToys --trackParameters "rgx{.*}" -t 10 --name _santanas_ --keepFailures -s 123456 (--robustHesse 1)

#Toy with signal
#combine -M MultiDimFit datacards/datacard_gg_M3000_R0p1.txt --rMin -5 --rMax 5 -S 0 --algo singles --cl=0.68 --saveFitResult --saveToys --trackParameters "rgx{.*}" -t 10 --name _santanas_ --keepFailures --expectSignal 1 -s 123456 (--robustHesse 1)

#Access fit results (higgsCombineTest.MultiDimFit.mH120.123456.root)
#limit->Scan("iToy:quantileExpected:r:covQual:trackedParam_r:trackedParamErr_r","quantileExpected==-1")

### Old ###
#combine -M FitDiagnostics /afs/cern.ch/work/s/santanas/Releases/CMSSW_9_4_0_trijet/src/CMSJET/RootTreeAnalyzer/output/prod1_26_09_18/fitMjetCategories4/datacards/datacard_gg_M3000_R0p1.txt --verbose 3 -S 0 --rMin -20 --rMax 20 (--robustHesse 1)
#Signal fits

##############################################################################################

# Definition of arrays
rootfile = [None] * ncategories

var = [None] * ncategories

th1_fromFile = [None] * ncategories
th1_original = [None] * ncategories
th1_rebin = [None] * ncategories
th1_rebin_bkg = [None] * ncategories
th1_rebin_pull = [None] * ncategories
rooHist = [None] * ncategories
rooHist_bkg = [None] * ncategories
numberOfEvents = [None] * ncategories

nbkg = [None] * ncategories
p1 = [None] * ncategories
p2 = [None] * ncategories
p3 = [None] * ncategories
bkgPdf = [None] * ncategories
bkgExtPdf = [None] * ncategories
bkgExtPdfTF1 = [None] * ncategories
fitResult = [None] * ncategories
ndof = [None] * ncategories
Chi2 = [None] * ncategories
reducedChi2 = [None] * ncategories


ParametricBkgPdf = [None] * ncategories



canvas = [None] * ncategories

## Create Workspace
w = RooWorkspace(workspaceName,workspaceName)
    # output tree with chi2
M1_br = array('f', [0])
R_br  = array('f', [0])
icat_br = array('i', [0])
Nevt  = array('i', [0])

ndof = array('i', [0])
chi2 = array('f', [0.])
rchi2 = array('f', [0.])
pval = array('f', [0.])
orig_pval = array('f', [0.])
pval_ok = array('i', [0])

ndof_allbins = array('i', [0])
chi2_allbins = array('f', [0.])
reducedchi2_allbins = array('f', [0.])

Npos   = array('f', [0.])
Nneg   = array('f', [0.])
WWruns = array('f', [0.])
WWpval = array('f', [0.])

fitrange_L = array('f', [0.])
fitrange_R = array('f', [0.])


## Loop over signal categories  
#var_min_limit = M1val*0.85
for icat, cat in enumerate(subDirList):
	#if ("2Muon" in str(cat)):	
	#	continue
        print("\n")
        print("#######################################################")
        print("Cat: "+str(icat))
        print("#######################################################")
    
        counter = 0
        #category_edges = cat.split("_")
        #P3mjet_low = float(category_edges[1])
        #P3mjet_high = float(category_edges[2])
        #Res2mjet_low = float(category_edges[4])
        #Res2mjet_high = float(category_edges[5])

        var_min_set = 453
        var_max_set = 4000
        
        ## Modify variable binning
        i_min_border = 0
        i_max_border = 0
        varBins = []
        for iborder, border in enumerate(varBins_all):
            if border<var_min_set:
                i_min_border += 1
                continue
            if border>var_max_set:
                continue
            else:
                varBins.append(float(border))
                i_max_border = iborder

        if fitRanges:
            print(var_min_set, var_max_set)
            if var_min_set < varBins[0] and abs(var_min_set - varBins[0])>0.001:
                varBins.insert(0, var_min_set)
            else:
                var_min_set = varBins[0]
            if var_max_set > varBins[-1]:
                varBins.append(var_max_set)
        else:
            var_min_set = varBins[0 ] 
            var_max_set = varBins[-1]
        print("var_min_set: ", str(var_min_set), "varBins[0]: ", str(varBins[0]), "varBins[1]: ", str(varBins[1]))

        NvarBins = len(varBins)-1
        canvas[icat] = TCanvas("canvas_"+cat, "canvas_"+cat, 200, 10, 700, 500 )


        #print(varBins)
        #prova=array('d',varBins)
        #print(prova)	
        ## Get original TH1 histogram from root file
        rootfilename = inputdir+"/"+cat+"/"+filenameInput
        #print rootfilename    
        rootfile[icat] = TFile.Open(rootfilename)
        print("Get "+histoname+" from file "+rootfilename)
        th1_fromFile[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
        th1_fromFile[icat].Draw()
        canvas[icat].SaveAs("tets.png")
        integral = int(th1_fromFile[icat].Integral())
        print("diocaro", integral)

	#outRoot=TFile( outputdir+"/background"+cat+".root","recreate")
	#th1_fromFile[icat].Write()
	#getattr(w,'import')(th1_fromFile[icat])
	#print("dai")
	#print("cazo")
	#clone_histo.Write()

        

      	  
        ## Create RooDataHist in fit range from TH1
        #rooHist[icat] = RooDataHist("RooDataHist_"+cat,"RooDataHist_"+cat,RooArgList(var[icat]),RooFit.Import(th1_rebin[icat]))
        #numberOfEvents[icat] = rooHist[icat].sumEntries()
  
        test=th1_fromFile[icat].ComputeIntegral() 
        print(test)
        ## Generate toy histogram
        gRandom = TRandom()
        if(generateToy==1):
            th1_original[icat] = TH1D("Toy","Toy", th1_fromFile[icat].GetNbinsX(), th1_fromFile[icat].GetXaxis().GetXmin(), th1_fromFile[icat].GetXaxis().GetXmax())
            #th1_original[icat] = TH1D("","", 5000, 0, 5000)
            gRandom.SetSeed(0)
            print("MAREMMA MAIALA")
            th1_original[icat].FillRandom(th1_fromFile[icat],integral)
            print("MAREMMA MAIALA")
            #th1_original[icat].Reset()
            #gRandom.SetSeed(0)
            #th1_original[icat].FillRandom(th1_fromFile[icat],integral)
	    #th1_original[icat].Write()
        else:
            th1_original[icat] = rootfile[icat].Get(histoname) # 1 GeV bin histogram
            # th1_original[icat].Add(th1_injected_2[icat], 1)
            # th1_original[icat].Add(th1_injected_1[icat], -2)

        var[icat] = RooRealVar(varname+"_"+cat,vartitle,varBins_all[0],varBins_all[-1])
        var[icat].Print()

        print("------------------------------------------------------------")
        print("                   FINAL FIT                                ")
        print("------------------------------------------------------------")

        var[icat] = RooRealVar(varname+"_"+cat,vartitle,var_min_set,var_max_set)
        var[icat].Print()

        RooFit.SumW2Error(kTRUE)
        ## Create data histogram with coarser binning
        th1_rebin_bkg[icat] = th1_original[icat].Rebin(NvarBins,"th1_rebin_"+cat,array('d',varBins))
        th1_rebin[icat] = th1_fromFile[icat].Rebin(NvarBins,"th1_rebin_"+cat,array('d',varBins))
        print(varBins)
        ## Create RooDataHist in fit range from TH1
        rooHist[icat] = RooDataHist("RooDataHist_"+cat,"RooDataHist_"+cat,RooArgList(var[icat]),RooFit.Import(th1_rebin[icat]))
        rooHist_bkg[icat] = RooDataHist("RooDataHist_bkg_"+cat,"RooDataHist_bkg_"+cat,RooArgList(var[icat]),RooFit.Import(th1_rebin_bkg[icat]))
        numberOfEvents[icat] = rooHist_bkg[icat].sumEntries()
        nbkg[icat] = rooHist[icat].sumEntries()
        th1_rebin[icat].SetBinErrorOption(TH1.kPoisson)
        th1_rebin_bkg[icat].SetBinErrorOption(TH1.kPoisson)

        ## Main physics observable defined in fit range
    ## Loop over signals
        dirRootFile=os.path.dirname(signalInputfilename)
        signalInputfile = io.open(signalInputfilename, "r")
        for iline, line in enumerate(signalInputfile):

            line = line.rstrip('\n')
            splitline = line.split(" ") 
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
            InSigRoot=TFile(dirRootFile+"/h1_mmuj_ak4__"+modell+"_"+cat+".root","w")
            #print(InSigRoot)
            test= "h1_mmuj_ak4__"+str(modell)
		    #print(test)
            th1Sig=InSigRoot.Get(test)
            sig_integral = int(th1Sig.Integral())
            sig_rebi = th1Sig.Rebin(NvarBins,"th1_rebin_"+cat,array('d',varBins))
        ## Create RooDataHist in fit range from TH1
            rooHistSign = RooDataHist("RooDataHist_"+modell+cat,"RooDataHist_"+modell+cat,RooArgList(var[icat]),RooFit.Import(sig_rebi))
		#print(th1Sig)
		#InSigRoot.Close()
		#outSigRoot=TFile(outputdir+"/"+modell+cat+".root","recreate")
		#print(outSigRoot)
		#th1Sig.SetName("sig")
		#th1Sig.Write()
		#sigclone.SetName("sig")
		#sigclone.Write()
		#clone_histo.Write()
		#data_obs.Write()
            getattr(w,'import')(rooHistSign)

        	## Signal pdf
        	

            os.system("mkdir -p "+outputdirdatacards+"/"+modell)    
            os.system("mkdir -p "+outputdirdatacards+"/"+modell+"/categories")    
        	    ## Create datacard for current signal
            datacardfilename = outputdirdatacards+"/"+modell+"/categories/"+"datacard_gen_"+signalString+".txt"
            outputfile = io.open(datacardfilename,'w')
        	
        	## Create datacard for current signal
        	#datacardfilename = outputdirdatacards+"/"+"datacard_"+signalString+".txt"
        	#outputfile = io.open(datacardfilename,'w')
        	
            outputfile.write( u"## Datacard for "+signalString+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"imax * number of channels"+"\n" )
            outputfile.write( u"jmax * number of backgrounds"+"\n" )
            outputfile.write( u"kmax * number of nuisance parameters (source of systematic uncertainties)"+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
                #outputfile.write( u"shapes * * " +str(outSigRoot.GetName())+ " $PROCESS\n" )
                #outputfile.write( u"shapes sig * " +str(outSigRoot.GetName())+ " $PROCESS\n" )
                #outputfile.write( u"shapes bkg * " +str(outSigRoot.GetName())+ " $PROCESS\n" )
                #outputfile.write( u"shapes data_obs * " +str(outSigRoot.GetName())+ " $PROCESS\n" )
                #outputfile.write( u"shapes bkg "+cat+" " +str(outSigRoot.GetName())+ " "+clone_histo.GetName()+"\n" )
                #outputfile.write( u"shapes data_obs "+cat+" " +str(outSigRoot.GetName())+ " "+data_obs.GetName()+"\n" )
            outputfile.write( u"shapes sig "+cat+" "+filenameworkspacePath+ " "+workspaceName+":"+rooHistSign.GetName()+"\n")
            outputfile.write( u"shapes bkg "+cat+" "+filenameworkspacePath+ " "+workspaceName+":"+rooHist[icat].GetName()+"\n")
            outputfile.write( u"shapes data_obs "+cat+" "+filenameworkspacePath+" "+workspaceName+":"+rooHist_bkg[icat].GetName()+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"bin \t\t"+cat+"\n" )
            outputfile.write( u"observation \t"+str(numberOfEvents[icat])+"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
            outputfile.write( u"bin \t\t"+cat+"\t"+cat+"\n" )
            outputfile.write( u"process \t"+"sig \t\t"+ "bkg"+"\n" )
            outputfile.write( u"process \t"+ "0" +" \t\t"+"1"+"\n" )
            outputfile.write( u"rate \t\t"+"-1"+ "\t\t"+"-1"+"\n" )
            #outputfile.write( u"process \t"+"0 \t\t"+ "1"+"\n" )
            #outputfile.write( u"rate \t\t"+str(int(nsig.getValV()))+" \t\t"+ str(int(nbkg[icat].getValV())) +"\n" )
            outputfile.write( u"----------------------------------------------------------------------------------------------------------------------------------"+"\n" )
        	#outputfile.write( p1[icat].GetName()+u" flatParam"+"\n" )
        	#outputfile.write( p2[icat].GetName()+u" flatParam"+"\n" )
        	#outputfile.write( p3[icat].GetName()+u" flatParam"+"\n" )
        	#outputfile.write( nbkg[icat].GetName()+u" flatParam"+"\n" )
		#outputfile.write( u"lumi lnN 1.018 -\n")
            outputfile.close()
    #workspacePath = workspacePath.split("/")[-1]
		#outSigRoot.Close()
        getattr(w,'import')(rooHist[icat])
        getattr(w,'import')(rooHist_bkg[icat])

        signalInputfile.close()

    ## Fill Workspace for data and background
    # Import data

   ############### Plots ##################

   ## Create canvas
        
  ##
        rootfile[icat].Close()

w.writeToFile(filenameworkspacePath)
w.Delete()
print("Workspace at "+filenameworkspacePath)

#### Create final datacards (by combinging all categories for a given signal) ####

## Loop over signals (find the modells)
listOfSignalModels = []

signalInputfile = io.open(signalInputfilename, "r")
for iline, line in enumerate(signalInputfile):
    
    line = line.rstrip('\n')
    splitline = line.split(" ")
    
    modell = splitline[0]
    #category = splitline[1]
    M1 = splitline[2]
    L = splitline[3]

    signalStringModel = modell#+"_"+"M"+str(int(float(M1)))+"_R0p"+str(int(float(R)*10))

    if signalStringModel not in listOfSignalModels:
        listOfSignalModels.append(signalStringModel)

signalInputfile.close()
print(listOfSignalModels)

datacardList = io.open(outputdirdatacards+"/datacardList_gen.txt", 'w')
datacardList_single = io.open(outputdirdatacards+"/datacardList_gen_single_category.txt", 'w')
## Loop over signals and categories (create combined datacard)
for signal in listOfSignalModels:
    currentPath = os.getcwd()
    command = ""
    os.system("mkdir -p "+outputdirdatacards+"/"+signal)
    ## Loop over event categories    
    for icat, cat in enumerate(subDirList):
        signalStringCat = signal+"_"+cat
        datacardfilenameAll = outputdirdatacards+"/"+signal+"/"+"datacard_gen_"+signal+".txt"
        datacardfilename = outputdirdatacards+"/"+signal+"/categories/"+"datacard_gen_"+signalStringCat+".txt"        

        print(signal)
        if icat==0:
            command = "python "+combineCardScriptPath

        if not os.path.isfile(datacardfilename):
            continue        

        ## Select only mjet categories within -2ssigma and +1ssigma of the Res2 mass
        #resolution = 0.1
        #mjetMin = float(cat.split("_")[1])
        #mjetMax = float(cat.split("_")[2])
        #mjetmean = (mjetMax+mjetMin)/2 
        #mass = float(signal.split("_")[1].split("M")[1])
        #R = float((signal.split("_")[2].split("R0p")[1]))*0.1
        #mass2 = mass*R
        #ssigma_mass2 = mass2*resolution
        #print(mjetMin, mjetMax, mjetmean, mass, R, mass2, ssigma_mass2
        #if ( 
        #    ( (mass2-mjetMax)<=2*ssigma_mass2 and (mass2-mjetMax)>=0 ) 
        #    or ( (mjetMin-mass2)<=1*ssigma_mass2 and (mjetMin-mass2)>=0 )
        #    or ( mass2>=mjetMin and mass2<mjetMax )  
        #    ): 
        command += " "+cat+"="+datacardfilename
        print(signalStringCat)

        splitline = signalStringCat.split("_")
        modell = splitline[0]
        M1val = (splitline[1]).strip("M")
        Lval  = (splitline[2]).strip("L")
        Lval =  Lval.replace("p",".")
        category = splitline[3]         
        datacardList_single.write(modell+" "+M1val+" "+Lval+" "+category+" "+datacardfilename+"\n")

    command += " -S > "+datacardfilenameAll
    print(command)
    os.system("cd "+outputdirdatacards+" ; "+command+" ; "+ "cd .." )
    print("Created final datacard at "+datacardfilenameAll)
   
    splitline = signal.split("_")
    modell = splitline[0]
    M1val = (splitline[1]).strip("M")
    Lval  = (splitline[2]).strip("L")
    Lval =  Lval.replace("p",".")

 
    datacardList.write(modell+" "+M1val+" "+Lval+" "+datacardfilenameAll+"\n")
    print("Created datacard List at "+outputdirdatacards+"/datacardList.txt")

#datacardList.close()


## Make multidimfit with fixed signal strength r=0

