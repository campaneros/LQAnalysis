#! /usr/bin/env python

#from ROOT import *
import ROOT as ROOT



function_STD= {
      "category1Muon_BDT_loose_btag" : "std_3par",
      "category1Muon_BDT_loose_nobtag" : "std_3par",
      "category2Muon_BDT_loose_btag" : "std_3par",
      "category2Muon_BDT_loose_nobtag" : "std_3par",
      "category2Muon_BDT_tight_btag" : "std_3par",
      "category2Muon_BDT_tight_nobtag" : "std_3par",
      "category1Muon_BDT_tight_btag" : "std_4par",
      "category1Muon_BDT_tight_nobtag" : "std_4par",
      }


function_UA2 = {
      "category1Muon_BDT_loose_btag" : "UA2_3par",
      "category1Muon_BDT_loose_nobtag" : "UA2_3par",
      "category2Muon_BDT_loose_btag" : "UA2_2par",
      "category2Muon_BDT_loose_nobtag" : "UA2_3par",
      "category2Muon_BDT_tight_btag" : "UA2_2par",
      "category2Muon_BDT_tight_nobtag" : "UA2_3par",
      "category1Muon_BDT_tight_btag" : "UA2_4par",
      "category1Muon_BDT_tight_nobtag" : "UA2_3par",
      }

nested_dict = {"function_STD": function_STD, "function_UA2": function_UA2}



limit_dict_1p0= {
    "700": 0.02, 
    "1000": 0.009,
    "2000": 0.002,
    "3000": 0.001,
    "4000": 0.00009,
    "5000" : 0.0008,
}

nested_dict_limit = {"limit_dict_1p0": limit_dict_1p0}

cross_section_dict_umu_0p1= {
    "700": 3.46008E-03, 
    "1000": 9.45814E-04,
    "2000": 4.85404E-05,
    "3000": 5.28391E-06,
    "4000": 7.33823E-07,
    "5000": 1.07973E-07, 
}

cross_section_dict_umu_1p0= {
    "700": 3.45971E-01, 
    "1000": 9.49567E-02,
    "2000": 4.97741E-03,
    "3000": 5.66914E-04,
    "4000": 8.74072E-05,
    "5000": 1.65938E-05,
} 


cross_section_dict_umu_1p5= {
    "700": 7.77979E-01, 
    "1000": 2.14614E-01,
    "2000": 1.15361E-02,
    "3000": 1.38123E-03,
    "4000": 2.35351E-04,
    "5000": 5.34450E-05, 
}

cross_section_dict_umu_2p0= {
    "700": 1.38027E+00, 
    "1000": 3.83304E-01,
    "2000": 2.12795E-02,
    "3000": 2.70481E-03,
    "4000": 5.11174E-04,
    "5000": 1.34026E-04, 
}

nested_dict_umu_cross_section = {"cross_section_dict_umu_0p1": cross_section_dict_umu_0p1, "cross_section_dict_umu_1p0": cross_section_dict_umu_1p0, "cross_section_dict_umu_1p5": cross_section_dict_umu_1p5, "cross_section_dict_umu_2p0": cross_section_dict_umu_2p0}

cross_section_dict_bmu_0p1= {
    "700": 2.73661E-04, 
    "1000": 5.27037E-05,
    "2000": 1.09773E-06,
    "3000": 5.87234E-08,
    "4000": 4.23152E-09,
    "5000": 2.71429E-10, 
}

cross_section_dict_bmu_1p0= {
    "700": 2.77177E-02, 
    "1000": 5.41827E-03,
    "2000": 1.24580E-04,
    "3000": 8.85890E-06,
    "4000": 1.32583E-06,
    "5000": 3.83615E-07,
} 


cross_section_dict_bmu_1p5= {
    "700": 6.33016E-02, 
    "1000": 1.25938E-02,
    "2000": 3.20714E-04,
    "3000": 2.81736E-05,
    "4000": 5.46531E-06,
    "5000": 1.84223E-06, 
}

cross_section_dict_bmu_2p0= {
    "700": 1.14575E-01, 
    "1000": 2.33070E-02,
    "2000": 6.67795E-04,
    "3000": 7.00758E-05,
    "4000": 1.58225E-05,
    "5000": 5.68700E-06, 
}

nested_dict_bmu_cross_section = {"cross_section_dict_bmu_0p1": cross_section_dict_bmu_0p1, "cross_section_dict_bmu_1p0": cross_section_dict_bmu_1p0, "cross_section_dict_bmu_1p5": cross_section_dict_bmu_1p5, "cross_section_dict_bmu_2p0": cross_section_dict_bmu_2p0}

nested_dict_cross_section = {"nested_dict_umu_cross_section": nested_dict_umu_cross_section, "nested_dict_bmu_cross_section": nested_dict_bmu_cross_section}    

def set_bkg_fit_function(icat, cat, var, fitFunction_name, fitparam, nbkg, bkgPdf, ParametricBkgPdf, bkgExtPdf, th1_rebin, numberOfEvents):

    del fitparam[:]

    if fitFunction_name == "std_2par": # 3 par (including normalization) standard dijet
        fitparam.append( ROOT.RooRealVar("p1_std"+cat,"p1_std"+cat,5,-50,50) )

        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*10)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "pow( (1-@0/13000), @1 )", ROOT.RooArgList(var[icat], fitparam[0]))  # 2 par (excl norm) standard dijet
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricBkgExtPdf_"+cat,"ParametricBkgExtPdf_"+cat,bkgPdf[icat],nbkg[icat]) 

    if fitFunction_name == "std_3par": # 3 par (including normalization) standard dijet
        fitparam.append( ROOT.RooRealVar("p1_std"+cat,"p1_std"+cat,7,-50,50) )
        fitparam.append( ROOT.RooRealVar("p2_std"+cat,"p2_std"+cat,5,  0,50) )
    
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*10)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "pow( (1-@0/13000), @1 )/pow(@0/13000, @2)", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1]))  # 2 par (excl norm) standard dijet
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0], fitparam[1]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricBkgExtPdf_"+cat,"ParametricBkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 

    if fitFunction_name == "std_4par": # 4 par (including normalization) standard dijet
        fitparam.append( ROOT.RooRealVar("p1_std"+cat,"p1_std"+cat,5,-50,50) )
        fitparam.append( ROOT.RooRealVar("p2_std"+cat,"p2_std"+cat,5,-50,50) )
        fitparam.append( ROOT.RooRealVar("p3_std"+cat,"p3_std"+cat,1,-10,10) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*3)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "pow( (1-@0/13000), @1 )/pow(@0/13000, @2+@3*log(@0/13000))", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1], fitparam[2]))  # 3 par (excl norm) standard dijet
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0], fitparam[1], fitparam[2]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricbkgExtPdf_"+cat,"bkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 

    if fitFunction_name == "std_5par": # 4 par (including normalization) standard dijet
        fitparam.append( ROOT.RooRealVar("p1_std"+cat,"p1_std"+cat,5,-50,50) )
        fitparam.append( ROOT.RooRealVar("p2_std"+cat,"p2_std"+cat,5,-50,50) )
        fitparam.append( ROOT.RooRealVar("p3_std"+cat,"p3_std"+cat,1,-10,10) )
        fitparam.append( ROOT.RooRealVar("p4_std"+cat,"p4_std"+cat,1,-10,10) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*3)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "pow( (1-@0/13000), @1 )/pow( @0/13000 , @2 + @3*log(@0/13000) + @4*pow( log( @0/13000) , 2 ) )", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1], fitparam[2], fitparam[3]))  # 4 par (excl norm) standard dijet
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0], fitparam[1], fitparam[2],fitparam[3]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricbkgExtPdf_"+cat,"bkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 
            
    if fitFunction_name == "modExp_2par": # 2 par (including normalization) exponential
        fitparam.append( ROOT.RooRealVar("p1_modExp"+cat,"p1_modExp"+cat,-35,-100,0) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*10)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "exp( @1*@0/13000 )", ROOT.RooArgList(var[icat], fitparam[0]))  # 1 par (excl norm) exponential pdf
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricbkgExtPdf_"+cat,"ParametricbkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 

    if fitFunction_name == "modExp_3par": # 3 par (including normalization) exponential
        fitparam.append( ROOT.RooRealVar("p1_modExp"+cat,"p1_modExp"+cat,  -35, -1000,100) )
        fitparam.append( ROOT.RooRealVar("p2_modExp"+cat,"p2_modExp"+cat, 0.53,   -1, 10) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*3)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "exp( @1*pow(@0/13000, @2) )", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1]))  # 2 par (excl norm) modified exponential pdf
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0], fitparam[1]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricBkgExtPdf_"+cat,"ParametricBkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 

    if fitFunction_name == "modExp_4par": # 4 par (including normalization) exponential
        fitparam.append( ROOT.RooRealVar("p1_modExp"+cat,"p1_modExp"+cat,-35,-100,100) )
        fitparam.append( ROOT.RooRealVar("p2_modExp"+cat,"p2_modExp"+cat,0.53,-1,10) )
        fitparam.append( ROOT.RooRealVar("p3_modExp"+cat,"p3_modExp"+cat,0.03,-1,10) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*10)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "exp( @1*(pow((@0/13000), @2) + pow((1-@0/13000), @3)) )", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1], fitparam[2]))
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0], fitparam[1], fitparam[2]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricbkgExtPdf_"+cat,"ParametricbkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 


    if fitFunction_name == "modExp_5par": # 5 par (including normalization) modified exponential
        fitparam.append( ROOT.RooRealVar("p1_modExp"+cat,"p1_modExp"+cat,-35,-100,100) )
        fitparam.append( ROOT.RooRealVar("p2_modExp"+cat,"p2_modExp"+cat,0.53,-100,100) )
        fitparam.append( ROOT.RooRealVar("p3_modExp"+cat,"p3_modExp"+cat,0.03,-100,100) )
        fitparam.append( ROOT.RooRealVar("p4_modExp"+cat,"p4_modExp"+cat,-10,-100,100) )
        
        nbkg[icat] = ROOT.RooRealVar("bkgPdf_"+cat+"_norm","bkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*10)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "exp( @1*pow((@0/13000), @2) + @3*pow((1-@0/13000), @4) )", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1], fitparam[2], fitparam[3]))
        bkgExtPdf[icat] = ROOT.RooExtendPdf("bkgExtPdf_"+cat,"bkgExtPdf_"+cat,bkgPdf[icat],nbkg[icat]) 

    if fitFunction_name == "UA2_2par": # 2 par (including normalization) UA2/ATLAS family
        fitparam.append( ROOT.RooRealVar("p1_UA2"+cat,"p1_UA2"+cat,5,-50,50) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*3)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "pow( @0/13000, -@1 )", ROOT.RooArgList(var[icat], fitparam[0]))
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricbkgExtPdf_"+cat,"bkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 
            
    if fitFunction_name == "UA2_3par": # 3 par (including normalization) UA2/ATLAS family
        fitparam.append( ROOT.RooRealVar("p1_UA2"+cat,"p1_UA2"+cat,5, 0,50) )
        fitparam.append( ROOT.RooRealVar("p2_UA2"+cat,"p2_UA2"+cat,15,-50,50) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*3)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "exp( -@2*(@0/13000) )/pow( @0/13000, @1)", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1]))
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0], fitparam[1]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricbkgExtPdf_"+cat,"bkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 
            
    if fitFunction_name == "UA2_4par": # 4 par (including normalization) UA2/ATLAS family
        fitparam.append( ROOT.RooRealVar("p1_UA2"+cat,"p1_UA2"+cat,5,-50,50) )
        fitparam.append( ROOT.RooRealVar("p2_UA2"+cat,"p2_UA2"+cat,15,-100,100) )
        fitparam.append( ROOT.RooRealVar("p3_UA2"+cat,"p3_UA2"+cat,35,-100,100) )
        
        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*10)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "exp( -@2*(@0/13000) -@3*pow( @0/13000, 2) )/pow( @0/13000, @1)", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1], fitparam[2]))
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0], fitparam[1]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricbkgExtPdf_"+cat,"bkgExtPdf_"+cat,ParametricBkgPdf[icat],nbkg[icat]) 

    if fitFunction_name == "modStd_3par": # 3 par (including normalization) standard dijet
        fitparam.append( ROOT.RooRealVar("p1_modStd"+cat,"p1_modStd"+cat,10,0,50) )
        fitparam.append( ROOT.RooRealVar("p2_modStd"+cat,"p2_modStd"+cat,0.5,0,5) )

        nbkg[icat] = ROOT.RooRealVar("ParametricBkgPdf_"+cat+"_norm","ParametricBkgPdf_"+cat+"_norm",numberOfEvents[icat],0,numberOfEvents[icat]*3)
        bkgPdf[icat] = ROOT.RooGenericPdf("bkgPdf_"+cat, "(1-@2*@0/13000)/pow(@0/13000, @1)", ROOT.RooArgList(var[icat], fitparam[0], fitparam[1]))  # 2 par (excl norm) standard dijet
        ParametricBkgPdf[icat] = ROOT.RooParametricShapeBinPdf( "ParametricBkgPdf_"+cat, "ParametricBkgPdf_"+cat, bkgPdf[icat], var[icat], ROOT.RooArgList(fitparam[0]), th1_rebin[icat])
        bkgExtPdf[icat] = ROOT.RooExtendPdf("ParametricBkgExtPdf_"+cat,"ParametricBkgExtPdf_"+cat,bkgPdf[icat],nbkg[icat]) 

def bkgRooPdf_to_TF1(icat, cat, fitFunction_name, fitparam, bkgExtPdfTF1, nbkg, var_min_set, var_max_set):
    if fitFunction_name == "std_2par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*pow( 1-x/13000, [1] )",var_min_set,var_max_set) 
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
            
    if fitFunction_name == "std_3par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*pow( 1-x/13000, [1] )/pow(x/13000, [2])",var_min_set,var_max_set) 
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
            
            
    if fitFunction_name == "std_4par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*pow( 1-x/13000, [1] )/pow(x/13000, [2]+[3]*log(x/13000))",var_min_set,var_max_set)
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        bkgExtPdfTF1[icat].SetParameter(3,fitparam[2].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    

    if fitFunction_name == "std_5par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*pow( 1-x/13000, [1] )/pow(x/13000, [2]+[3]*log(x/13000)+[4]*pow(log(x/13000),2))",var_min_set,var_max_set)
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        bkgExtPdfTF1[icat].SetParameter(3,fitparam[2].getValV())
        bkgExtPdfTF1[icat].SetParameter(4,fitparam[3].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
            
    if fitFunction_name == "modExp_2par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*exp( [1] * x/13000 )",var_min_set,var_max_set) # 1 par (excl norm) exponential pdf  
        
        bkgExtPdfTF1[icat].SetParameter(0,1)     
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
            
    if fitFunction_name == "modExp_3par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*exp( [1]*pow(x/13000, [2]) )",var_min_set,var_max_set) # 2 par (excl norm) modified exponential pdf 
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    

    if fitFunction_name == "modExp_4par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*exp( [1]*(pow((x/13000), [2]) + pow((1-x/13000), [3])) )",var_min_set,var_max_set) # 3 (excl norm) par modified exponential pdf
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        bkgExtPdfTF1[icat].SetParameter(3,fitparam[2].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
        
    if fitFunction_name == "modExp_5par": # 4 (excl norm) par modified exponential pdf
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*exp( [1]*pow((x/13000), [2]) + [3]*pow((1-x/13000), [4]) )",var_min_set,var_max_set)
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        bkgExtPdfTF1[icat].SetParameter(3,fitparam[2].getValV())
        bkgExtPdfTF1[icat].SetParameter(4,fitparam[3].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
                
    if fitFunction_name == "UA2_2par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat, "[0]*pow( x/13000, -[1] )", var_min_set, var_max_set)
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
        
    if fitFunction_name == "UA2_3par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat, "[0]*exp( -[2]*(x/13000) )/pow( x/13000, [1])", var_min_set, var_max_set)
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
            
    if fitFunction_name == "UA2_4par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat, "[0]*exp( -[2]*(x/13000) -[3]*pow( x/13000, 2) )/pow( x/13000, [1])", var_min_set, var_max_set)
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        bkgExtPdfTF1[icat].SetParameter(3,fitparam[2].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    
        
    if fitFunction_name == "modStd_3par":
        bkgExtPdfTF1[icat] = ROOT.TF1("bkgExtPdfTF1_"+cat,"[0]*( 1-[2]*x/13000)/pow(x/13000, [1])",var_min_set,var_max_set) 
        
        bkgExtPdfTF1[icat].SetParameter(0,1)    
        bkgExtPdfTF1[icat].SetParameter(1,fitparam[0].getValV())
        bkgExtPdfTF1[icat].SetParameter(2,fitparam[1].getValV())
        integral_TF1 = bkgExtPdfTF1[icat].Integral(var_min_set,var_max_set)
        bkgExtPdfTF1[icat].SetParameter(0,float(nbkg[icat].getValV())/float(integral_TF1))    


def evaluate_chi2(icat, ndof, Chi2, rChi2, ndof_allbins, Chi2_allbins, reducedChi2_allbins, th1_rebin, th1_rebin_pull, Npar):
    ndof[0] = 0
    Chi2[0] = 0. 
    rChi2[0] = 0.

    ndof_allbins[0] = 0 
    Chi2_allbins[0] = 0. 
    reducedChi2_allbins[0] = 0.

    ## Get Chi2 from final fit
    for bin in range(th1_rebin[icat].GetNbinsX()):
        if bin == 0:
            continue

        data = float(th1_rebin[icat].GetBinContent(bin))

        # pull histo
        if data!=0:
            pull = th1_rebin_pull[icat].GetBinContent(bin)

            #Chi2 for all not empty bins
            ndof_allbins[0] += 1
            Chi2_allbins[0] += pull*pull
            
            #chi2
            if(data>10):
                Chi2[0]+=pull*pull
                ndof[0]+=1

    ndof[0] -= Npar
    if ndof[0]>0:
        rChi2[0] = Chi2[0]/float(ndof[0])
    ndof_allbins[0] -= Npar
    if ndof_allbins[0]>Npar:
        reducedChi2_allbins[0] = Chi2_allbins[0]/float(ndof_allbins[0])

def evaluate_chi2_lowstat(icat, ndof, Chi2, rChi2, ndof_allbins, Chi2_allbins, reducedChi2_allbins, th1_rebin, th1_rebin_pull, Npar):
    ndof[0] = 0
    Chi2[0] = 0. 
    rChi2[0] = 0.

    ndof_allbins[0] = 0 
    Chi2_allbins[0] = 0. 
    reducedChi2_allbins[0] = 0.

    ## Get Chi2 from final fit
    for bin in range(th1_rebin[icat].GetNbinsX()):
        if bin == 0:
            continue

        data = float(th1_rebin[icat].GetBinContent(bin))

        # pull histo
        if data!=0:
            pull = th1_rebin_pull[icat].GetBinContent(bin)

            #Chi2 for all not empty bins
            ndof_allbins[0] += 1
            Chi2_allbins[0] += pull*pull
            
            #chi2
            if(data<100):
                Chi2[0]+=pull*pull
                ndof[0]+=1

    ndof[0] -= Npar
    if ndof[0]>0:
        rChi2[0] = Chi2[0]/float(ndof[0])
    ndof_allbins[0] -= Npar
    if ndof_allbins[0]>Npar:
        reducedChi2_allbins[0] = Chi2_allbins[0]/float(ndof_allbins[0])

def parse_syst_file( syst_22cat, syst_9cat, syst_1cat, syst_dir ):
    
    f_22cat = open(syst_dir+"/syst_22cat.txt")
    f_22cat_lines = f_22cat.readlines()

    for line in f_22cat_lines:
        line = line.rstrip('\n')
        splitline  = line.split(" ")
        systlabel  = splitline[0]        
        splitline.pop(0)
        systparam  = splitline        
        syst_22cat[systlabel] = systparam

    f_9cat = open(syst_dir+"/syst_9cat.txt")
    f_9cat_lines = f_9cat.readlines()

    for line in f_9cat_lines:
        line = line.rstrip('\n')
        splitline  = line.split(" ")
        systlabel  = splitline[0]        
        splitline.pop(0)
        systparam  = splitline        
        syst_9cat[systlabel] = systparam

    f_1cat = open(syst_dir+"/syst_1cat.txt")
    f_1cat_lines = f_1cat.readlines()

    for line in f_1cat_lines:
        line = line.rstrip('\n')
        splitline  = line.split(" ")
        systlabel  = splitline[0]        
        splitline.pop(0)
        systparam  = splitline        
        syst_1cat[systlabel] = systparam
        #print systlabel, systparam

