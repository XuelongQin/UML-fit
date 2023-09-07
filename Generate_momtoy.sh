#!/bin/bash
statbin=(0 1 2 3 5 7)
istat=${1}
q2stat=${statbin[${istat}]}

Fsinlist=(0.00 0.04 0.08)
iFsin=${2}
Fsin=${Fsinlist[${iFsin}]}

echo q2stat is ${q2stat}
echo Fsin is ${Fsin}


export CMSSWDIR=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/CMSSW_10_4_0/src
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
export WORKDIR=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/fit/latest/UML-fit
cd $WORKDIR
if [ ! -r $WORKDIR/Generate_momtoy ]; then
    echo $WORKDIR/Generate_momtoy not found
    exit 1
fi
echo ./Generate_momtoy 4 1 1 100 ${q2stat} 8 0 ${Fsin} 2016 2017 2018
./Generate_momtoy 4 1 1 100 ${q2stat} 8 0 ${Fsin} 2016 2017 2018
