#!/usr/bin/bash
PARSED_OPTIONS=$(getopt -n "$0" -o t:n:df --long "toys,number:datacard,fit"  -- "$@")


eval set -- "$PARSED_OPTIONS"

dodatac=0
dofit=0
toy=1
n=1

while true;
do
    case "$1" in
        -t|--toys)
            if [ -n "$2" ]; then
                toy=$2
                echo "Generating #${toy} toys"
            fi
            shift 2;;
        -n|--number)
            if [ -n "$2" ]; then
                n=$2
                echo "Running  chi2 on #${n} toys"
            fi
            shift 2;;

        -d|--datacard)
            dodatac=1
            echo "creating datacard"
            shift;;
        -f|--fit)
            dofit=1
            echo "fitting bkg"
            shift;;
        --)
        shift
        break;;
    esac
done





Mass=(1000 2000 3000)
L=(0p1 1p0)
fit=(std_4par)
expect_signal=(0 0.01 0.001 0.0001)


source ../../../setup_cmssw.sh
cmsenv


if [ $dodatac -eq 1 ]; then 
    python create_datacard.py
fi

if [ $toy -eq 1 ]; then
    nplot=${toy}
else
    nplot=$((${toy}/100))    
    #echo ${nplot}
fi

for f in "${fit[@]}"
do
    if [ $dofit -eq 1 ]; then
        python create_categories.py -F ${f}
    fi 


    for m in "${Mass[@]}"
    do 
        #echo ${m}
        for l in "${L[@]}"
        do
        for exp in "${expect_signal[@]}"
            do
                python bias_study.py -g /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M${m}_L${l}/datacard_gen_LQumu_M${m}_L${l}.txt -d /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M${m}_L${l}/datacard_${f}_LQumu_M${m}_L${l}.txt -t ${toy} --expectSignal ${exp} -l 0.0001 --gen 1
                python plot_bias.py -i  /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_${f}_LQumu_M${m}_L${l}_t_${toy}_syst0_seed123456/higgsCombine_toys${toy}_expectSignal${exp}_${f}.MultiDimFit.mH120.123456.root -o /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_${f}_LQumu_M${m}_L${l}_t_${toy}_syst0_seed123456 -f expect${exp} -r ${exp}
                python plotSimultaneousFit_toy.py -t /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_${f}_LQumu_M${m}_L${l}_t_${toy}_syst0_seed123456/higgsCombine_toys${toy}_expectSignal${exp}.GenerateOnly.mH120.123456.root -f /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_${f}_LQumu_M${m}_L${l}_t_${toy}_syst0_seed123456/higgsCombine_toys${toy}_expectSignal${exp}_${f}.MultiDimFit.mH120.123456.root -c /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M${m}_L${l}/categories/ -w /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/workspace_${f}.root -n ${nplot}  -o /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/plotSimFit_${f}/LQumu_M${m}_L${l} -F ${f} -b expect_signal${exp}_${f}
                python FitToy_Wrapper.py -t /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_${f}_LQumu_M${m}_L${l}_t_${toy}_syst0_seed123456/higgsCombine_toys${toy}_expectSignal${exp}.GenerateOnly.mH120.123456.root -f /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/datacard_${f}_LQumu_M${m}_L${l}_t_${toy}_syst0_seed123456/higgsCombine_toys${toy}_expectSignal${exp}_${f}.MultiDimFit.mH120.123456.root -c /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/datacards/LQumu_M${m}_L${l}/categories/ -w /data/mcampana/CMS/CMSSW_8_1_0_LQ/src/Fit_Signal/output_MC/workspace_${f}.root -n ${n} -r ${exp}
            done
        done
    done
done