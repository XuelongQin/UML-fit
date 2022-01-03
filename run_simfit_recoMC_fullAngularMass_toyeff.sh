#!/bin/bash

par=0

ntoys=10
nsam=1

plot=0
save=1

ibin=${1}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-JpsiFit
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
export OUTDIR=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave

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

for iy in {2016..2018}
do	      
    dataname="${SAMPLEDIR}/${iy}/lmnr/newphi/recoMCDataset_b${bin}_${iy}.root"
    if [ ! -r "${dataname}" ]; then
	echo "${dataname}" not found
	exit 1
    fi
    cp "${dataname}" .
done

if [ ! -r $HOME/simfit_recoMC_fullAngularMass_toyeff ]; then
    echo $HOME/simfit_recoMC_fullAngularMass_toyeff not found
    exit 1
fi
cp $HOME/simfit_recoMC_fullAngularMass_toyeff .

mkdir simFitResults4d

if [ "$bin" -lt 8 ]; then
    echo ./simfit_recoMC_fullAngularMass_toyeff ${bin} ${par} ${ntoys} ${nsam} 1 ${plot} ${save} 2016 2017 2018
    ./simfit_recoMC_fullAngularMass_toyeff ${bin} ${par} ${ntoys} ${nsam} 1 ${plot} ${save} 2016 2017 2018
else
    echo ./simfit_recoMC_fullAngularMass_toyeff ${bin} ${par} ${ntoys} ${nsam} 1 ${plot} ${save} 2016
    ./simfit_recoMC_fullAngularMass_toyeff ${bin} ${par} ${ntoys} ${nsam} 1 ${plot} ${save} 2016
fi

outfiledir="$OUTDIR/simFitResults4d"
if [ ! -d "${outfiledir}" ]; then
    mkdir "${outfiledir}"
fi
cp simFitResults4d/* "${outfiledir}/"

rm -rf simFitResults4d

rm simfit_recoMC_fullAngularMass_toyeff
rm recoMCDataset_b*
