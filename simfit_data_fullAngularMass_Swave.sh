#!/bin/bash

par=1

multi=0
nsam=0
q2stat=0

plot=0
save=2

bin=${1}
yearConf=${2}

eosdir="/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave"

# Create directories for fit logs, results and plots
if [ ! -d "${eosdir}/logs_simFit4d_freeSigmas" ]; then mkdir ${eosdir}/logs_simFit4d_freeSigmas; fi
if [ ! -d simFitResults4d ]; then mkdir simFitResults4d; fi
if [ ! -d plotSimFit4d_d ]; then mkdir plotSimFit4d_d; fi



# Compile dictionary and macro
# make AngDict
if make simfit_data_fullAngularMass_Swave; then

    case "$yearConf" in

	0)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2016 2017 2018
	    nohup ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2016 2017 2018 \
		&> ${eosdir}/logs_simFit4d_freeSigmas/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_0_${plot}_${save}_2016_2017_2018.log &
	    ;;

	1)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2016
	    ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2016 \
		&> ${eosdir}/logs_simFit4d_freeSigmas/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_0_${plot}_${save}_2016.log &
	    ;;

	2)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2017
	    ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2017 \
		&> ${eosdir}/logs_simFit4d_freeSigmas/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_0_${plot}_${save}_2017.log &
	    ;;

	3)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2018
	    ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} 0 ${plot} ${save} 2018 \
		&> ${eosdir}/logs_simFit4d_freeSigmas/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_0_${plot}_${save}_2018.log &
	    ;;

    esac

fi
