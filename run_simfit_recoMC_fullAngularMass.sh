#!/bin/bash

par=0

multi=0
nsam=${1}

plot=0
save=1

ibin=${2}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-JpsiFit
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
export EFFDIR=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/files
export OUTDIR=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

nbin=0
while read -a line; do
    abin[$nbin]=${line[0]}
    nbin=$((nbin+1))
done < $HOME/../confSF/KDE_SF.list
bin=${abin[$ibin]}

echo 'now submitting for bin ' ${bin}

if [ "${par}" == 0 ]; then
    parstr="ev"
else
    parstr="od"
fi

for iy in {2016..2018}
do	      
    dataname="${SAMPLEDIR}/${iy}/lmnr/newphi/recoMCDataset_b${bin}_${iy}.root"
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

if [ ! -r $HOME/simfit_recoMC_fullAngularMass ]; then
    echo $HOME/simfit_recoMC_fullAngularMass not found
    exit 1
fi
cp $HOME/simfit_recoMC_fullAngularMass .

mkdir simFitResults4d
mkdir plotSimFit4d_d

if [ "$bin" -lt 8 ]; then
    echo ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018
    ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018
else
    echo ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016
    ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016
fi


if [ ! -d $OUTDIR/simFitResults4d ]; then
    mkdir $OUTDIR/simFitResults4d
fi
if [ ! -d $OUTDIR/plotSimFit4d_d ]; then
    mkdir $OUTDIR/plotSimFit4d_d
fi
cp plotSimFit4d_d/* $OUTDIR/plotSimFit4d_d/
cp simFitResults4d/* $OUTDIR/simFitResults4d/
