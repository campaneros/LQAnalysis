#!/usr/bin/bash

L=( 0p1 1p0 1p5 2p0 )
#type=( limits limits_bmu )
#type=( limits limits_bmu )
type=( limits )


for t in ${type[@]}; do
    for l in ${L[@]}; do
        echo "Running on ${t} ${l}"
        python make_table.py  limit_${t}_cut_${l}.txt limit_8cat_${t}_${l}.txt limit_8cat_${t}_${l}.txt
    done
done
