# UML-fit

# How to run

**Table of contents**  
[Setup working area](#setup)  
[Create datasets](#createDatasets)  
[Fit to signal MC](#performFit)  
* [Fit to angular variables](#angular)  
* [Fit to angular + mass variables](#angMass)  
[Plotting macros](#plotMacros)  

<a name="setup"/>

## Setup working area
Make sure to be in an area with ROOT v6.12 (or later) available. The code could work with previous versions, but it is not guaranteed.
If you are on CERN lxplus machines, you can achieve it by running:
```sh
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src/ && cmsenv && cd ../..
```
Clone this branch in the working directory:
```sh
git clone -b master git@github.com:CMSKStarMuMu/UML-fit.git
cd UML-fit
```
<a name="createDatasets"/>

## Create datasets
If needed, change the [location of the ntuples](createDataset.cc#L55-L62), which need to be produced with the code in the [B0KstMuMuNtuple repository](https://github.com/CMSKStarMuMu/B0KstMuMuNtuple).  
Then, you can produce the files which will contain the needed datasets, that is, correctly and wrongly tagged reco'ed events from the MC. The datasets will include the following variables: ctK, ctL, phi, mass, rand.
rand is a random variable uniformly generated between 0 and 1, to be used to select a given statistics from the MC sample.  
Example on how to run for 2017 ntuples:
```sh
root -q -b 'createDataset.cc(7)'
```


<a name="performFit"/>

## Perform fits to MC dataset with PDF*eff

<a name="angular"/>
### Fit to the angular variables only
The fit is performed by the simfit_recoMC_singleComponent code.
It requires in input:
* a root file containing the workspace with the roodatasets to be fitted (produced in the step above)
* the files containing the KDE modeling of the efficiency and, if possible, the corresponding partial integrals. These are produced by the workflow described in the eff-KDE repository. 

Compile and run with:
```sh
source simfit_recoMC_singleComponent.sh
```
where you have to set the datasets to be considered (set year = 0 to not include the dataset). 
The variable "datalike" sets the statistics to be considered:  
* datalike = 0 means consider the full MC stat (half of it actually)  
* datalike = 1 means consider a data-like statistics.  

The code will produce a root file `simFitResults/simfitResult_recoMC_singleComponentXXXX.root` containing the RooFitResult objects, where XXXX describes the considered datasets.
Corresponding fit projection plots are created in `plotSimFit_d/simfitResult_recoMC_singleComponent_*.pdf`.

#### Plot and compare fit results
```sh
root -b -q 'plotSimFitResults.cc(1)' # for fit with odd efficiency on even dataset and full MC stat
```
while use
```
root -b -q 'plotSimFitResults.cc(1,-1,true, false, false,true)' 
```
to plot the results for the data-like stat fits.

<a name="angMass"/>

### Fit to the angular variables + mass variable
The fit is performed by the simfit4d_recoMC_singleComponent code.
It requires in input:
* a root file containing the workspace with the roodatasets to be fitted (produced in the step above)
* the files containing the KDE modeling of the efficiency and, if possible, the corresponding partial integrals. These are produced by the workflow described in the eff-KDE repository. 
* a root file containing the workspace with the RooFitResults to the MC mass distributions (as produced by [this code](https://github.com/CMSKStarMuMu/selection_and_fits/blob/master/perBin_massFit.py))

Compile and run with:
```sh
source simfit4d_recoMC_singleComponent.sh
```
where you have to set the datasets to be considered (set year = 0 to not include the dataset). 
The variable "datalike" sets the statistics to be considered (datalike = 0 -> full MC stat, datalike = 1 -> data-like statistics).  

The code will produce a root file `simFit4dResults/simfit4dResult_recoMC_singleComponentXXXX.root` containing the RooFitResult objects, where XXXX describes the considered datasets.
Corresponding fit projection plots are created in `plotSimFit4d_d/simfitResult_recoMC_singleComponent_*.pdf`.

<a name="plotMacros"/>

## Plotting macros

Here a list of the macros used to plot the fit results:
- `plot_simfit_*.cc`: plot the projections of a single fit result, reading the workspace containing PDF and dataset. These plots can also be produced during the fit by activating the plot flag
- `plotSimFitResults.cc`: crate a graph comparing the angular parameters resulting from one or more fits to RECO MC sample and compare them to the GEN-fit result, as a function of the q2 bin. It also prints out the formatted table for the "efficiency mismodelling" systematics. It reads the workspace with GEN fit result and the TTrees with the RECO fit results
- `plotMultiFit.cc`: plot distributions of the resulting angular parameters and their uncertainties from fits to a set of data-like samples, and compare them with the full-sample fit results. It also prints out the formatted table for the "fit bias" systematics. It reads the TTrees with the fit results
- `plotMultiDataFit.cc`: same as above for set of fits to control-region subsamples, comparing them with the full control-region result
- `plotSimFitComparison.cc`: create a graph with comparison of data simultaneous fit results, with single-year results and the mass constraints. It reads the workspace with resulting PDFs
- `plotSimFitComparison_manual.cc`: same as above, but read the simultaneous fit result from the log file (needed for Jpsi control region)

### Outdated macros

To use the following macros with the latest version of the fit results, some changes are necessary:
- `plotMultiResults.cc`: plot 1D and 2D distribution of set of MC subsample fits, together with boundary slice and slice of the log-likelihood from full-MC fit. It allows also to plot 2D slices of individual requirements composing the boundary
- `plotPulls.cc`: plot the pull plot of a set of toy fit results and compute the coverage of the confidence interval (not recently tested, it may still work with the latest version of [simfit_toy_fullAngular.cc](simfit_toy_fullAngular.cc)
