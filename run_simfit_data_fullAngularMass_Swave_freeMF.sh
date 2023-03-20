#!/bin/bash

par=1

multi=0
nsam=${1}
q2stat=${4}

plot=1
save=${5}

bin=${2}
# ibin=${2}

yearConf=${3}

XGBv=8

export HOME=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/sys/UML-fit
export CMSSWDIR=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/CMSSW_10_4_0
export SAMPLEDIR=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/sys/UML-fit/dataset
export EFFDIR=/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/eff-XGBv8

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

# nbin=0
# while read -a line; do
#     abin[$nbin]=${line[0]}
#     nbin=$((nbin+1))
# done < $HOME/../confSF/KDE_SF.list
# bin=${abin[$ibin]}

echo 'now submitting for bin ' ${bin}

if [ "${par}" == 0 ]; then
    parstr="ev"
else
    parstr="od"
fi

for iy in {2016..2018}
do
    [ "$yearConf" -gt 0 ] && [ "$((${yearConf}+2015))" != "$iy" ] && continue
    dataname="${SAMPLEDIR}/recoDATADataset_b${bin}_${iy}.root"
    # dataname="${SAMPLEDIR}/${iy}/lmnr/newphi/recoDATADataset_b${bin}_${iy}.root"
    effname="${EFFDIR}/KDEeff_b${bin}_${parstr}_${iy}.root"
    if [ ! -r "${dataname}" ]; then
	echo "${dataname}" not found
	exit 1
    fi
    if [ ! -r "${effname}" ]; then
	echo "${effname}" not found
	exit 1
    fi
    cp "${dataname}" .
    cp "${effname}" .
done

if [ ! -r $HOME/simfit_data_fullAngularMass_Swave_freeMF ]; then
    echo $HOME/simfit_data_fullAngularMass_Swave_freeMF not found
    exit 1
fi
cp $HOME/simfit_data_fullAngularMass_Swave_freeMF .
cp $HOME/*.pcm .

mkdir simFitResults4d
mkdir plotSimFit4d_d

case "$yearConf" in

    0)
	echo ./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2016 2017 2018
	./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2016 2017 2018
	;;

    1)
	echo ./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2016
	./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2016
	;;

    2)
	echo ./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2017
	./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2017
	;;

    3)
	echo ./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2018
	./simfit_data_fullAngularMass_Swave_freeMF ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} 1 ${plot} ${save} 2018
	;;

esac

if [ ! -d $HOME/simFitResults4d ]; then
    mkdir $HOME/simFitResults4d
fi
if [ ! -d $HOME/plotSimFit4d_d ]; then
    mkdir $HOME/plotSimFit4d_d
fi
cp plotSimFit4d_d/* $HOME/plotSimFit4d_d/
cp simFitResults4d/* $HOME/simFitResults4d/
