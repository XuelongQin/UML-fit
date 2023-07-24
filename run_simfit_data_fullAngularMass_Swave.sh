#!/bin/bash

par=1

multi=0
nsam=0

XGBv=8
localFile=0
fitopt=0
unbl=4

plot=1
save=2

xgbv=8

bin=4
unbl=4

yearConf=0

wscale=${1}

if [ "${USER}" == "aboletti" ]; then
    export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-JpsiFit
    export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
elif [ "${USER}" == "fiorendi" ]; then
    export HOME=/afs/cern.ch/work/f/fiorendi/private/effKDE/UML-fit
    export CMSSWDIR=/afs/cern.ch/work/f/fiorendi/private/effKDE/CMSSW_10_4_0/src
else
    echo no user found
    exit 1
fi

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

if [ ! -r $HOME/simfit_data_fullAngularMass_Swave ]; then
    echo $HOME/simfit_data_fullAngularMass_Swave not found
    exit 1
fi
cp $HOME/simfit_data_fullAngularMass_Swave .
cp $HOME/*.pcm .

mkdir simFitResults4d
mkdir plotSimFit4d_d

case "$yearConf" in
# void simfit_data_fullAngularMass_SwaveBin(int q2Bin, int parity, bool multiSample, uint nSample, uint q2stat, int fitOption, int XGBv, int unblind, bool localFiles, bool plot, int save, std::vector<int> years)

    0)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016 2017 2018
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016 2017 2018
	;;

    1)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016
	;;

    2)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2017
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2017
	;;

    3)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2018
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${wscale} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2018
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
