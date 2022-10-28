#!/bin/bash

par=1

multi=0
nsam=0

plot=1
save=2

vers=${1}

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults ]; then mkdir -p simFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
if make simfit_recoMC_fullAngular; then
 
    # while read -a line; do
	# bin=${line[0]}

    bin=4
	
    for year in {2016..2018}; do
	
	./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${vers} 0 ${plot} ${save} ${year} \
	    &>logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${multi}_${nsam}_${vers}_${year}.out &
	
    done

    ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${vers} 0 ${plot} ${save} 2016 2017 2018 \
	&>logs_simFit/simfit_recoMC_fullAngular_XGBv4_${bin}_${par}_${multi}_${nsam}_${vers}_2016_2017_2018.out &

    # done < ../confSF/KDE_SF.list

fi
