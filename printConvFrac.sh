#!/bin/bash

stat=64

while [ "${stat}" -le 128 ]
do
    for bin in {0..5}
    do
	# echo "=== bin ${bin} ==="
	for i in {0..2}
	do
	    arr[$((3*${bin}+${i}))]=$(grep "statCode - ff ${i}" logs_parSub_fracTest/simfit_recoMC_fullAngularMass_*_${bin}_${stat}.out | wc -l)
	done
    done
    echo ${arr[@]}
    stat=$((${stat}*2))
done
