#!/bin/bash

par=0

neff=${1}
nsam=0

xgb=8

plot=0
save=1

ibin=${2}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-JpsiFit
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
export OUTDIR=/eos/user/a/aboletti/BdToKstarMuMu/MCstat

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

nbin=0
while read -a line; do
    [ ${line[0]} == "4" ] && continue
    [ ${line[0]} == "6" ] && continue
    abin[$nbin]=${line[0]}
    nbin=$((nbin+1))
done < $HOME/../confSF/KDE_SF.list
bin=${abin[$ibin]}

echo 'now submitting for bin ' ${bin}

# for iy in {2016..2018}
# do	      
#     dataname="${SAMPLEDIR}/${iy}/lmnr/newphi/recoMCDataset_b${bin}_${iy}.root"
#     if [ ! -r "${dataname}" ]; then
# 	echo "${dataname}" not found
# 	exit 1
#     fi
#     cp "${dataname}" .
# done

if [ ! -r $HOME/simfit_recoMC_fullAngularMass_toybkg_toyeff ]; then
    echo $HOME/simfit_recoMC_fullAngularMass_toybkg_toyeff not found
    exit 1
fi
cp $HOME/simfit_recoMC_fullAngularMass_toybkg_toyeff .
cp $HOME/*.pcm .

mkdir simFitResults4d
mkdir plotSimFit4d_d

if [ "$bin" -lt 8 ]; then
    echo ./simfit_recoMC_fullAngularMass_toybkg_toyeff ${bin} ${par} ${neff} ${nsam} ${xgb} 0 ${plot} ${save} 2016 2017 2018
    ./simfit_recoMC_fullAngularMass_toybkg_toyeff ${bin} ${par} ${neff} ${nsam} ${xgb} 0 ${plot} ${save} 2016 2017 2018
else
    echo ./simfit_recoMC_fullAngularMass_toybkg_toyeff ${bin} ${par} ${neff} ${nsam} ${xgb} 0 ${plot} ${save} 2016
    ./simfit_recoMC_fullAngularMass_toybkg_toyeff ${bin} ${par} ${neff} ${nsam} ${xgb} 0 ${plot} ${save} 2016
fi

outfiledir="$OUTDIR/simFitResults4d/xgbv8"
if [ ! -d "${outfiledir}" ]; then
    mkdir -p "${outfiledir}"
fi
cp simFitResults4d/* "${outfiledir}/"

outplotdir="$OUTDIR/plotSimFit4d_d/xgbv8"
if [ ! -d "${outplotdir}" ]; then
    mkdir -p "${outplotdir}"
fi
cp plotSimFit4d_d/* "${outplotdir}/"
