import os
import argparse
import itertools
from multiprocessing import Pool
import create_workspaces_and_datacards_utils as cwd_utils





def run_bias_study(args):
    m, l, exp, toy , category = args
    print(m,l,exp,toy,category)
    fitting_function = ['STD','UA2']
    for function_gen, function_fit in itertools.product(fitting_function, fitting_function):
     fit_fucntion_name_gen= cwd_utils.nested_dict["function_%s"%function_gen]["%s"%category]
     fit_fucntion_name_fit= cwd_utils.nested_dict["function_%s"%function_fit]["%s"%category]
     os.system("python bias_study.py -d /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/new_umu_datacards_finalcat_"+function_gen+"_450/umuLQumu_M"+str(m)+"_L"+str(l)+"/categories/datacard_"+fit_fucntion_name_gen+"_umuLQumu_M"+str(m)+"_L"+str(l)+"_"+category+".txt -g /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/new_umu_datacards_finalcat_"+function_fit+"_450/umuLQumu_M"+str(m)+"_L"+str(l)+"/categories/datacard_"+fit_fucntion_name_fit+"_umuLQumu_M"+str(m)+"_L"+str(l)+"_"+category+".txt -t "+str(toy)+" --expectSignal "+str(exp)+" -l 0.001 --gen 1 -o /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/gen_"+function_gen+"_fit_"+function_fit+"")
     os.system("python plot_bias.py -i ../Fit_Signal/output_MC/gen_"+function_gen+"_fit_"+function_fit+"/datacard_"+fit_fucntion_name_gen+"_umuLQumu_M"+str(m)+"_L"+str(l)+"_"+category+"_t_"+str(toy)+"_syst0_seed123456/higgsCombine_toys"+str(toy)+"_expectSignal"+str(exp)+"_std_4par.MultiDimFit.mH120.123456.root -o ../Fit_Signal_BDT_tests/output_MC/test_STD_2/ -f fit_"+function_fit+"_gen_"+function_gen+"_toys"+str(toy)+"_"+category+"_expect"+str(exp)+"_"+str(m)+"_"+str(l)+" -r "+str(exp)+"")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--toys', '-t', type=int, default=1000)
    parser.add_argument('--number', '-n', type=int, default=1)
    parser.add_argument('--datacard', '-d', action='store_true')
    parser.add_argument('--fit', '-f', action='store_true')
    args = parser.parse_args()

    #Mass = [4000, 5000]
    Mass = [700, 1000, 2000, 3000, 4000, 5000]
    #Mass = [700]
    L = ['0p1','1p0','1p5','2p0']
    expect_signal = [0]
    #expect_signal = [0, 0.01, 0.001, 0.0001]
    #categories = ['category2Muon_BDT_loose_btag']
    categories = ['category2Muon_BDT_loose_btag','category2Muon_BDT_loose_nobtag','category2Muon_BDT_tight_nobtag','category2Muon_BDT_tight_btag','category1Muon_BDT_loose_btag','category1Muon_BDT_loose_nobtag','category1Muon_BDT_tight_nobtag','category1Muon_BDT_tight_btag']


    if args.datacard:
        os.system("python create_datacard.py")

    if args.fit:
        os.system("python create_categories.py -F std_4par")

    pool = Pool(processes=15)
    tasks = list(itertools.product(Mass, L, expect_signal, [args.toys],categories))
    pool.map(run_bias_study, tasks)


if __name__ == '__main__':
    main()
