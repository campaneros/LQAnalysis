#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import subprocess
from glob import glob
from collections import defaultdict
from collections import OrderedDict

from ROOT import *

cmsswbase = os.environ['CMSSW_BASE']

usage = "usage: To be run from RootTreeAnalyzer/scripts/fit/trijetAna:  python launch_bias_study.py -f [--fitfunctiondirlist] \"[datacard_f1,datacard_f2,...]\" -o output/bias_study -w weboutputdir --expectSignal xsec_limits.root --limit xsec_limits.root --resubmit"

parser = optparse.OptionParser(usage)

parser.add_option("-o", "--output", dest="outputdir",
                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

#parser.add_option("-w", "--weboutput", dest="weboutputdir",
#                  help="Web directory containing bias pull plots.")

parser.add_option("-f", "--funcdirlist", dest="fitfunctiondirlist",
                  help="fit function datacard directory")

parser.add_option("-t", "--toys", dest="toys", default=100,
                  help="number fo toys")

parser.add_option("--expectSignal", dest="expectSignal", default=None,
                  help="file with limits")

parser.add_option("-l", "--limit", dest="limit", default=None,
                  help="file with limits")

parser.add_option("--resubmit", action="store_true", dest="resubmit", default=False,
                  help="file with limits")

(opt, args) = parser.parse_args()

if  not opt.outputdir:
    parser.error("output directory not provided.")

if  not opt.fitfunctiondirlist:
    parser.error("Fit function datacad directory not provided.")


######################################################################
gROOT.SetBatch

# Get expected limits from output of plotLimits_2d_th2.py (xsec_limits.root)
if opt.expectSignal:
    print opt.expectSignal
    limitfilename = opt.expectSignal
    limitfile = TFile.Open(limitfilename, "r")
    limitth2  = limitfile.Get("h2_limit")

if opt.limit:
    print opt.limit
    limitfilename = opt.limit
    limitfile2 = TFile.Open(limitfilename, "r")
    limitth22  = limitfile2.Get("h2_limit")

for gendir in opt.fitfunctiondirlist.split(","):
    for fitdir in opt.fitfunctiondirlist.split(","):

        # Make list of signal samples to analyze
        print gendir, fitdir
        samplelist   = []
        g_samplelist = []
        f_samplelist = []

        g_samplelist = sorted(next(os.walk(gendir))[1])
        f_samplelist = sorted(next(os.walk(fitdir))[1])
        if g_samplelist == f_samplelist:
            samplelist = g_samplelist
        else:
            print "[ERROR]: gen sample list does not match fit sample list"
            print list(set(g_samplelist)-set(f_samplelist))
            print list(set(f_samplelist)-set(g_samplelist))
            print g_samplelist
            print f_samplelist
        genfuncname = ""
        for string in gendir.split("/"):
            if string.startswith("datacards_"):
                splitstring = string.split("_")
                genfuncname = splitstring[1]+"_"+splitstring[2]
        #print genfuncname

        fitfuncname = ""
        for string in fitdir.split("/"):
            if string.startswith("datacards_"):
                splitstring = string.split("_")
                fitfuncname = splitstring[1]+"_"+splitstring[2]
        #print fitfuncname

        outputdir        = opt.outputdir
        #weboutputdir     = opt.weboutputdir
        combineoutputdir = "%s/gen_%s_fit_%s" % ( outputdir, genfuncname, fitfuncname )
        logoutputdir     = "%s/log"           % ( outputdir ) 
        scriptsdir       = "%s/sh"            % ( outputdir )

        os.system("mkdir -p "+outputdir)
        #os.system("mkdir -p "+weboutputdir)
        scriptsPath = os.path.dirname(os.path.abspath(__file__))+"/../../"
        #os.system("cp "+scriptsPath+"/index.php "+weboutputdir)
        os.system("mkdir -p "+combineoutputdir)
        os.system("mkdir -p "+logoutputdir)
        os.system("mkdir -p "+scriptsdir)

for isample, sample in enumerate(samplelist):
    outputdir        = opt.outputdir
    #weboutputdir     = opt.weboutputdir
    logoutputdir     = "%s/log"           % ( outputdir ) 
    scriptsdir       = "%s/sh"            % ( outputdir )

    tmpsubmitfilename = "%s/bias_%s.sh"  % ( scriptsdir, sample )
    tmpsubmitfile = open(tmpsubmitfilename,'w')
    tmpsubmitfile.write("#!/bin/bash\n")
    tmpsubmitfile.write("cd "+cmsswbase+"/src/CMSJET/RootTreeAnalyzer/scripts/fit/trijetAna\n")
    tmpsubmitfile.write("cmsenv\n")
    tmpsubmitfile.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
    tmpsubmitfile.write("eval `scramv1 runtime -sh`\n")
    tmpsubmitfile.close()

    to_resubmit = 0
    for gendir in opt.fitfunctiondirlist.split(","):
        for fitdir in opt.fitfunctiondirlist.split(","):

            # Make list of signal samples to analyze
            print gendir, fitdir
            genfuncname = ""
            for string in gendir.split("/"):
                if string.startswith("datacards_"):
                    splitstring = string.split("_")
                    genfuncname = splitstring[1]+"_"+splitstring[2]
                    break
            #print genfuncname

            fitfuncname = ""
            for string in fitdir.split("/"):
                if string.startswith("datacards_"):
                    splitstring = string.split("_")
                    fitfuncname = splitstring[1]+"_"+splitstring[2]
                    break
            #print fitfuncname

            combineoutputdir = "%s/gen_%s_fit_%s" % ( outputdir, genfuncname, fitfuncname )
    
            M1sample = float( (sample.split("_")[1]).replace("M"  ,""  ) )
            Rsample  = float( (sample.split("_")[2]).replace("R0p","0.") )
            limit    = 0
            if opt.expectSignal:
                ibin  = limitth2.FindBin( M1sample, Rsample )
                limit = limitth2.GetBinContent(ibin)
            if opt.limit:
                ibin2  = limitth22.FindBin( M1sample, Rsample )
                limit2 = limitth22.GetBinContent(ibin2)
            print M1sample, Rsample, limit
            gen_sampledatacard = "%s/datacard_%s.txt" % ( gendir, sample )
            fit_sampledatacard = "%s/datacard_%s.txt" % ( fitdir, sample )

            if opt.resubmit:
                biasoutputfileName = "%s/datacard_%s_t_%s_syst1_seed1/higgsCombine_toys%s_expectSignal%s.MultiDimFit.mH120.123456.root"  % ( combineoutputdir, sample, str(opt.toys), str(opt.toys), str(limit) )
                biasoutputfile = TFile.Open(biasoutputfileName)
                if biasoutputfile:
                    if biasoutputfile.Get("limit"):
                        continue
                    else:
                        to_resubmit = 1
                else:
                    to_resubmit = 1
                        
            command = "python bias_study.py"
            command+= " -g %s"                % ( gen_sampledatacard )
            command+= " -d %s"                % ( fit_sampledatacard )
            command+= " -o %s"                % ( combineoutputdir ) 
            command+= " -t %s"                % ( str(opt.toys) )
            command+= " -S 1"
            command+= " -l %s"                % ( str(limit2) )
            command+= " -s 12"
            if limit!=0:
                command+= " --expectSignal %s" % ( str(limit) )
            print "\n"+command+"\n"

            plot_command = "python plot_bias.py"
            plot_command+= " -i %s/datacard_%s_t_%s_syst1_seed12345/higgsCombine_toys%s_expectSignal%s.MultiDimFit.mH120.123456.root"  % ( combineoutputdir, sample, str(opt.toys), str(opt.toys), str(limit) )
            plot_command+= " -o %s"                                                                               % ( outputdir )
            plot_command+= " -f gen_%s_fit_%s_%s"                                                                 % ( genfuncname, fitfuncname, sample )
            plot_command+= " -r %s"                                                                               % ( str(limit) )
            print plot_command

            if to_resubmit == 1:
                os.system(command)
                os.system(plot_command)

	    os.system(plot_command)
