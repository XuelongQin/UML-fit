#!/bin/bash

par=0

## full stat
multi=0
nsam=0
#### data like one job
# multi=0
# nsam=${1}

plot=1
save=1

outdir="/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/"

# Create directories for fit logs, results and plots
if [ ! -d "${outdir}logs_simFit4d" ]; then mkdir ${outdir}logs_simFit4d; fi
if [ ! -d "${outdir}simFitResults4d" ]; then mkdir ${outdir}simFitResults4d; fi
if [ ! -d "${outdir}plotSimFit4d_d" ]; then mkdir ${outdir}plotSimFit4d_d; fi

# Compile dictionary and macro
# make AngDict
if make simfit_recoMC_fullAngularMass; then

    while read -a line; do
	bin=${line[0]}

	# for year in {2016..2018}; do
	
	#     ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 0 ${plot} ${save} ${year} \
	# 	&>logs_simFit/simfit_recoMC_fullAngularMass_${bin}_${par}_${multi}_${nsam}_${year}.out &
	
	# done

	if [ "$bin" -lt 8 ]; then
	    nohup ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 0 ${plot} ${save} 2016 2017 2018 \
		&>${outdir}logs_simFit4d/simfit_recoMC_fullAngularMass_${bin}_${par}_${multi}_${nsam}_2016_2017_2018.out &
	else
	    nohup ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 0 ${plot} ${save} 2016 \
		&>${outdir}logs_simFit4d/simfit_recoMC_fullAngularMass_${bin}_${par}_${multi}_${nsam}_2016.out &
	fi

    done < ../confSF/KDE_SF.list

fi
