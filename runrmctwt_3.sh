#!/bin/bash

export CMSSWDIR=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/fit/CMSSW_10_4_0/src
export WORKDIR=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/fit/UML-fit

cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd $WORKDIR

for q2Bin in 4 6
do
    echo ./simfit_data_fullAngularMass_Swave_rmctwt ${q2Bin} 1 0 0 0 1 1 2016 2017 2018
    ./simfit_data_fullAngularMass_Swave_rmctwt ${q2Bin} 1 0 0 0 1 1 2016 2017 2018
done
