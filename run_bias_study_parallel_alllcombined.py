import os
import argparse
import itertools
from multiprocessing import Pool


def run_bias_study(args):
    m, l, exp, toy = args
    os.system("python bias_study.py -d /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/newsample_datacards_finalcatSTD450_std_family_real/umuLQumu_M{}_L{}/datacard_std_3par_umuLQumu_M{}_L{}.txt -g /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal_BDT_data_all/newsample_datacards_finalcatSTD450_std_family_real/umuLQumu_M{}_L{}/datacard_std_3par_umuLQumu_M{}_L{}.txt -t {} --expectSignal {} -l 0.01 --gen 1".format(m, l, m, l, m, l, m, l, toy, exp))
    os.system("python plot_bias.py -i ../Fit_Signal/output_MC/datacard_std_3par_umuLQumu_M{}_L{}_t_{}_syst0_seed123456/higgsCombine_toys{}_expectSignal{}_std_4par.MultiDimFit.mH120.123456.root -o ../Fit_Signal_BDT_tests/output_MC/test_STD/ -f toys_{}_expect{}_{}_{} -r {}".format(m, l, toy, toy, exp, toy, exp, m, l, exp))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--toys', '-t', type=int, default=1000)
    parser.add_argument('--number', '-n', type=int, default=1)
    parser.add_argument('--datacard', '-d', action='store_true')
    parser.add_argument('--fit', '-f', action='store_true')
    args = parser.parse_args()

    Mass = [700, 2000, 3000, 4000, 5000]
    L = ['1p0']
    expect_signal = [0]
    #expect_signal = [0, 0.01, 0.001, 0.0001]

    if args.datacard:
        os.system("python create_datacard.py")

    if args.fit:
        os.system("python create_categories.py -F std_4par")

    pool = Pool()
    tasks = list(itertools.product(Mass, L, expect_signal, [args.toys]))
    pool.map(run_bias_study, tasks)


if __name__ == '__main__':
    main()
