import os
import argparse
import itertools
from multiprocessing import Pool
import create_workspaces_and_datacards_utils as cwd_utils

def run_bias_study(args):
    m, l, toy = args
    fitting_function = ['STD','UA2']
    exp = 3*float(cwd_utils.nested_dict_limit["limit_dict_%s"%l]["%s"%m])
    for function_gen, function_fit in itertools.product(fitting_function, fitting_function):
     if function_gen == "STD":
        fit_fucntion_name_gen = "std_3par"
     else:
        fit_fucntion_name_gen = "UA2_3par"
     if function_fit == "STD":
        fit_fucntion_name_fit = "std_3par"
     else:
        fit_fucntion_name_fit = "UA2_3par"
     #fit_fucntion_name_gen= cwd_utils.nested_dict["function_%s"%function_gen]["%s"%category]
     #fit_fucntion_name_fit= cwd_utils.nested_dict["function_%s"%function_fit]["%s"%category]
     os.system("python bias_study.py -d /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/new_umu_datacards_finalcat_"+function_gen+"_450/umuLQumu_M"+str(m)+"_L"+str(l)+"/datacard_"+fit_fucntion_name_gen+"_umuLQumu_M"+str(m)+"_L"+str(l)+".txt -g /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/new_umu_datacards_finalcat_"+function_fit+"_450/umuLQumu_M"+str(m)+"_L"+str(l)+"/datacard_"+fit_fucntion_name_fit+"_umuLQumu_M"+str(m)+"_L"+str(l)+".txt -t "+str(toy)+" --expectSignal "+str(round(exp,4))+" -l 0.001 --gen 1 -o /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/umu_gen_"+function_gen+"_fit_"+function_fit+"")
     os.system("python plot_bias.py -i ../Fit_Signal/output_MC/umu_gen_"+function_gen+"_fit_"+function_fit+"/datacard_"+fit_fucntion_name_gen+"_umuLQumu_M"+str(m)+"_L"+str(l)+"_t_"+str(toy)+"_syst0_seed123456/higgsCombine_toys"+str(toy)+"_expectSignal"+str(round(exp,4))+"_std_4par.MultiDimFit.mH120.123456.root -o ../Fit_Signal_BDT_tests/output_MC/test_STD_2/ -f umu_fit_"+function_fit+"_gen_"+function_gen+"_toys"+str(toy)+"_expect"+str(round(exp,4))+"_"+str(m)+"_"+str(l)+" -r "+str(round(exp,4))+" -a 1")
    

#def run_bias_study(args):
#    m, l, exp, toy = args
#    os.system("python bias_study.py -d /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/new_umu_datacards_finalcatSTD450_std_family_real/umuLQumu_M{}_L{}/datacard_std_3par_umuLQumu_M{}_L{}.txt -g /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/new_umu_datacards_finalcatSTD450_std_family_real/umuLQumu_M{}_L{}/datacard_std_3par_umuLQumu_M{}_L{}.txt -t {} --expectSignal {} -l 0.01 --gen 1".format(m, l, m, l, m, l, m, l, toy, exp))
#    os.system("python plot_bias.py -i ../Fit_Signal/output_MC/datacard_std_3par_umuLQumu_M{}_L{}_t_{}_syst0_seed123456/higgsCombine_toys{}_expectSignal{}_std_4par.MultiDimFit.mH120.123456.root -o ../Fit_Signal_BDT_tests/output_MC/test_STD/ -f toys_{}_expect{}_{}_{} -r {}".format(m, l, toy, toy, exp, toy, exp, m, l, exp))
#

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--toys', '-t', type=int, default=1000)
    parser.add_argument('--number', '-n', type=int, default=1)
    parser.add_argument('--datacard', '-d', action='store_true')
    parser.add_argument('--fit', '-f', action='store_true')
    args = parser.parse_args()

    Mass = [700, 1000, 2000, 3000, 4000, 5000]
    L = ['0p1','1p0','1p5','2p0']
    #expect_signal = [0]
    #expect_signal = [0, 0.01, 0.001, 0.0001]

    if args.datacard:
        os.system("python create_datacard.py")

    if args.fit:
        os.system("python create_categories.py -F std_4par")

    pool = Pool()
    tasks = list(itertools.product(Mass, L, [args.toys]))
    pool.map(run_bias_study, tasks)


if __name__ == '__main__':
    main()
