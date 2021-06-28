#!/bin/bash

par=1

multi=0
nsam=${1}

plot=0
save=1

ibin=${2}

yearConf=${3}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-JpsiFit
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE

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

if [ ! -r $SAMPLEDIR/2016/lmnr/newphi/recoDATADataset_b${bin}_2016.root ]; then
    echo $SAMPLEDIR/2016/lmnr/newphi/recoDATADataset_b${bin}_2016.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2017/lmnr/newphi/recoDATADataset_b${bin}_2017.root ]; then
    echo $SAMPLEDIR/2017/lmnr/newphi/recoDATADataset_b${bin}_2017.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2018/lmnr/newphi/recoDATADataset_b${bin}_2018.root ]; then
    echo $SAMPLEDIR/2018/lmnr/newphi/recoDATADataset_b${bin}_2018.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2016/lmnr/newphi/KDEeff_b${bin}_od_2016.root ]; then
    echo $SAMPLEDIR/2016/lmnr/newphi/KDEeff_b${bin}_od_2016.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2017/lmnr/newphi/KDEeff_b${bin}_od_2017.root ]; then
    echo $SAMPLEDIR/2017/lmnr/newphi/KDEeff_b${bin}_od_2017.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2018/lmnr/newphi/KDEeff_b${bin}_od_2018.root ]; then
    echo $SAMPLEDIR/2018/lmnr/newphi/KDEeff_b${bin}_od_2018.root not found
    exit 1
fi
if [ ! -r $HOME/simfit_data_fullAngularMass_Swave ]; then
    echo $HOME/simfit_data_fullAngularMass_Swave not found
    exit 1
fi

cp $SAMPLEDIR/2016/lmnr/newphi/recoDATADataset_b${bin}_2016.root .
cp $SAMPLEDIR/2017/lmnr/newphi/recoDATADataset_b${bin}_2017.root .
cp $SAMPLEDIR/2018/lmnr/newphi/recoDATADataset_b${bin}_2018.root .
cp $SAMPLEDIR/2016/lmnr/newphi/KDEeff_b${bin}_od_2016.root .
cp $SAMPLEDIR/2017/lmnr/newphi/KDEeff_b${bin}_od_2017.root .
cp $SAMPLEDIR/2018/lmnr/newphi/KDEeff_b${bin}_od_2018.root .
cp $HOME/simfit_data_fullAngularMass_Swave .

mkdir simFitResults4d
mkdir plotSimFit_d

case "$yearConf" in

    0)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018
	;;

    1)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016
	;;

    2)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2017
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2017
	;;

    3)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2018
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2018
	;;

esac

if [ ! -d $HOME/simFitResults4d ]; then
    mkdir $HOME/simFitResults4d
fi
if [ ! -d $HOME/plotSimFit_d ]; then
    mkdir $HOME/plotSimFit_d
fi
cp plotSimFit_d/* $HOME/plotSimFit_d/
cp simFitResults4d/* $HOME/simFitResults4d/
# for file in simFitResults/* ; do cp $file $HOME/${file//.root/_${multi}s${nsam}.root}; done

rm -rf plotSimFit_d
rm -rf simFitResults4d

rm simfit_data_fullAngularMass_Swave
rm recoDATADataset_b*
rm KDEeff_b*
